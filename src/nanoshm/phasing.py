"""Separate a minor haplotype / paralogue / contaminant read population.

The problem this solves
-----------------------
Competitive alignment removes reads whose *best* hit is a contaminant, but not a
population that still maps best to the target: a second haplotype, a paralogous
locus, or a near-identical contaminant. Those reads carry a coherent block of
differences from the reference, and pooling them into one pileup makes every
position in that block look like a low-frequency variant.

Measured on the JH4 data (docs/calibration/FINDINGS.md): 14% of Aicda-/- reads
carried *all sixteen* of the worst-offending positions while 83% carried none.
Removing that population accounts for roughly a sevenfold drop in the apparent
background.

Why co-occurrence, not frequency
--------------------------------
The obvious approach, "flag positions with a suspiciously high alt fraction",
does not work. The same ~60-site block appears in all six JH4 samples but its
abundance ranges from ~2% of reads in one mouse to ~17% in another, so any
frequency cutoff either misses it or swallows real signal.

Co-occurrence is the invariant. Somatic hypermutation does not co-occur: each
read comes from a different B cell and carries its own mutations, so two SHM
positions are statistically independent. A haplotype block is the opposite, the
same reads carry every site. So we seed on one variable position and keep only
the sites that vary *in the same reads*.
Why not drop every block (changed October 2026)
-----------------------------------------------
Clonal expansion breaks the independence argument above. Reads from one
expanded germinal-centre clone share that clone's mutations, so a clone looks
exactly like a haplotype block. On the bundled wild-type JH4 files, de novo
phasing removed 19% of Wt3 reads on a 4-site block and 6% of Wt4 on a 10-site
block, at different positions in each mouse and absent from every Aicda-/-
mouse: real SHM. A haplotype, paralogue or contaminant is the opposite: it
does not need AID, so it also appears in mutation-free controls, and it recurs
across animals.

So by default (phase_mode="auto") a block is only *dropped* when it matches a
block already seen in control reads or recurring in another sample
(`known_blocks`). Any other block is *flagged* in the stats and left in.
phase_mode="drop" restores the old drop-everything behaviour.
"""
from __future__ import annotations

from collections import defaultdict

from .config import Params

BASES = ("A", "C", "G", "T")


def _pileup(aligned, lo: int, hi: int):
    counts: dict = defaultdict(lambda: defaultdict(int))
    depth: dict = defaultdict(int)
    for r in aligned:
        for rpos, base, _q in r.pairs:
            if lo <= rpos <= hi and base in BASES:
                counts[rpos][base] += 1
                depth[rpos] += 1
    return counts, depth


def candidate_sites(aligned, ref: str, p: Params, lo: int, hi: int) -> list[tuple[int, str]]:
    """Positions variable enough to be worth testing for co-occurrence.

    The floor is deliberately low (phase_min_freq, default 1%) because a minor
    haplotype can sit at a few percent. Specificity comes from the co-occurrence
    test in `partition`, not from this threshold.
    """
    counts, depth = _pileup(aligned, lo, hi)
    sites = []
    for pos in range(lo, hi + 1):
        d = depth.get(pos, 0)
        if d == 0:
            continue
        refbase = ref[pos].upper()
        if refbase not in BASES:
            continue
        for alt in BASES:
            if alt == refbase:
                continue
            frac = counts[pos].get(alt, 0) / d
            if p.phase_min_freq <= frac < p.germline_freq_threshold:
                sites.append((pos, alt, frac))
    sites.sort(key=lambda x: -x[2])
    return [(pos, alt) for pos, alt, _ in sites]


def _index(aligned, positions: set) -> list[dict]:
    """One pass: for each read, its observed base at each candidate position.

    Built once and reused, otherwise block detection is O(reads x sites x
    alignment length) and takes minutes on a few thousand reads.
    """
    out = []
    for r in aligned:
        out.append({rpos: b for rpos, b, _q in r.pairs if rpos in positions})
    return out


def _block_from_seed(idx: list[dict], sites, seed) -> list:
    """Given a seed site, find the sites that vary in the same reads.

    A site joins the block if most seed-carrying reads have it and few others
    do. That asymmetry separates a haplotype from sequencing error (rare in both
    groups) and from real SHM (sporadic in both).
    """
    spos, salt = seed
    inside = [a for a in idx if a.get(spos) == salt]
    outside = [a for a in idx if spos in a and a.get(spos) != salt]
    if len(inside) < 5:
        return []

    def frac(group, pos, alt):
        n = c = 0
        for a in group:
            b = a.get(pos)
            if b is not None:
                n += 1
                c += (b == alt)
        return (c / n) if n else 0.0

    return [(pos, alt) for pos, alt in sites
            if frac(inside, pos, alt) >= 0.5 and frac(outside, pos, alt) <= 0.05]


def find_blocks(aligned, ref: str, p: Params, lo: int, hi: int,
                idx: list[dict] | None = None, sites: list | None = None,
                max_seeds: int = 10) -> list[list[tuple[int, str]]]:
    """All distinct co-occurring blocks, largest first.

    Seeds on the most variable sites in turn; a seed that already belongs to a
    found block is skipped, so a clone and a haplotype in the same sample come
    back as two blocks rather than the larger one hiding the other.
    """
    if sites is None:
        sites = candidate_sites(aligned, ref, p, lo, hi)
    if len(sites) < p.phase_min_sites:
        return []
    if idx is None:
        idx = _index(aligned, {pos for pos, _ in sites})
    blocks: list[list] = []
    claimed: set = set()
    for seed in sites[:max_seeds]:
        if seed in claimed:
            continue
        block = _block_from_seed(idx, sites, seed)
        if len(block) >= p.phase_min_sites:
            blocks.append(block)
            claimed.update(block)
    blocks.sort(key=len, reverse=True)
    return blocks


def matches_known(block, known_blocks, min_shared: int) -> bool:
    """A block matches if it shares most of its sites with a known block."""
    bs = set(map(tuple, block))
    for kb in known_blocks or []:
        ks = set(map(tuple, kb))
        shared = len(bs & ks)
        if shared >= min(min_shared, len(bs), len(ks)) and shared >= 0.5 * min(len(bs), len(ks)):
            return True
    return False


def _carriers(idx: list[dict], block, p: Params) -> list[int]:
    out = []
    for i, seen in enumerate(idx):
        covered = sum(1 for pos, _alt in block if pos in seen)
        if covered < p.phase_min_cover:
            continue
        carried = sum(1 for pos, alt in block if seen.get(pos) == alt)
        if carried / covered >= p.phase_score_threshold:
            out.append(i)
    return out


def _mode(p: Params) -> str:
    return "off" if not p.phase_reads else p.phase_mode


def partition(aligned, ref: str, p: Params, lo: int, hi: int,
              known_blocks: list | None = None):
    """Split reads into (major, minor) populations. Returns (major, stats).

    Deliberately conservative: blocks that cannot be shown to be non-SHM are
    reported, not removed; and if removing the blocks that can would discard
    more than phase_max_drop of the reads, nothing is dropped and the reason is
    reported. A tool that silently discarded half the data would be worse than
    one that does nothing.
    """
    aligned = list(aligned)
    mode = _mode(p)
    stats = {"phasing_applied": False, "phase_mode": mode, "haplotype_sites": 0,
             "reads_dropped_minor": 0, "blocks": [], "reason": ""}
    if mode == "off":
        stats["reason"] = "disabled"
        return aligned, stats
    if not aligned:
        stats["reason"] = "no aligned reads"
        return aligned, stats

    sites = candidate_sites(aligned, ref, p, lo, hi)
    if len(sites) < p.phase_min_sites:
        stats["reason"] = f"only {len(sites)} variable sites (< {p.phase_min_sites})"
        return aligned, stats
    idx = _index(aligned, {pos for pos, _ in sites})
    blocks = find_blocks(aligned, ref, p, lo, hi, idx=idx, sites=sites)
    if not blocks:
        stats["reason"] = f"no co-occurring block of >= {p.phase_min_sites} sites"
        return aligned, stats

    drop: set = set()
    dropped_sites = 0
    for block in blocks:
        carriers = _carriers(idx, block, p)
        if mode == "drop":
            action, why = "drop", "phase_mode=drop"
        elif mode == "auto" and matches_known(block, known_blocks, p.phase_min_sites):
            action, why = "drop", "matches a block seen in controls or another sample"
        else:
            action, why = "flag", ("not seen in controls or other samples; may be a clone, "
                                   "kept")
        # how many of the block's sites each read carries, as a histogram. A
        # haplotype or paralogue is all-or-nothing (reads carry every site or
        # none); hypermutation and even a clone are graded. This is the
        # evidence behind the action taken, and viz_project draws it.
        carried: dict = {}
        for seen in idx:
            covered = sum(1 for pos, _alt in block if pos in seen)
            if covered < p.phase_min_cover:
                continue
            k = sum(1 for pos, alt in block if seen.get(pos) == alt)
            carried[k] = carried.get(k, 0) + 1
        stats["blocks"].append({
            "sites": [[pos + 1, alt] for pos, alt in sorted(block)],   # 1-based for humans
            "n_sites": len(block), "reads": len(carriers), "action": action, "why": why,
            "carried_hist": dict(sorted(carried.items()))})
        if action == "drop" and carriers:
            drop.update(carriers)
            dropped_sites += len(block)

    if not drop:
        n = len(blocks)
        stats["reason"] = (f"{n} co-occurring block(s) flagged, none dropped "
                           f"(phase_mode={mode})")
        return aligned, stats
    if len(drop) / len(aligned) > p.phase_max_drop:
        stats["reason"] = (f"minor population is {100 * len(drop) / len(aligned):.0f}% of reads "
                           f"(> {100 * p.phase_max_drop:.0f}%); refusing to phase")
        return aligned, stats

    major = [r for i, r in enumerate(aligned) if i not in drop]
    stats.update(phasing_applied=True, reads_dropped_minor=len(drop),
                 haplotype_sites=dropped_sites,
                 reason=f"removed {len(drop)} reads sharing a "
                        f"{dropped_sites}-site block")
    return major, stats
