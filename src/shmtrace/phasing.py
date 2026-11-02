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

Which B cells carry a block (changed November 2026)
---------------------------------------------------
Where the amplicon reads into the V(D)J junction (junction_ref_start), each read
names the B-cell clone it came from, and that settles the question the rule above
can only guess at. A clone's shared mutations sit in reads of one junction (or a
few, where sequencing error splits it). A block carried by reads of many unrelated
junctions cannot be hypermutation: independent B cells do not acquire the same four
or more mutations. So in "auto" mode a block whose carriers are polyclonal
(`polyclonal`) is dropped even when it is seen in no other sample, which a second
haplotype or a contaminant present in one mouse only would otherwise survive; a
block held by one clone is kept and said to be a clone. On the bundled JH4
libraries every block in a wild-type mouse is held by one or two junctions and is
kept, as before.
"""
from __future__ import annotations

from collections import Counter, defaultdict

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
            # with a declared ancestral reference a selected variant may reach fixation and still be real
            hi = 2.0 if p.ancestral_reference else p.germline_freq_threshold
            if p.phase_min_freq <= frac < hi:
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
    tried = 0
    # max_seeds counts the seeds tried, not the sites passed over: a site already in a
    # found block is skipped free, or a block of ten sites would hide every rarer one
    for seed in sites:
        if tried >= max_seeds:
            break
        if seed in claimed:
            continue
        tried += 1
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


# a block is polyclonal when at least POLY_MIN_READS of its carriers have a readable
# junction, they fall in at least POLY_MIN_CLONES junction groups, and the two largest
# hold under POLY_TOP2 of them: a clone, even one split by sequencing error, is not that
POLY_MIN_READS, POLY_MIN_CLONES, POLY_TOP2 = 8, 4, 0.5


def junctions(aligned, carriers: list[int], p: Params) -> dict | None:
    """How the carriers of a block spread over V(D)J junctions (clones.py), or None
    where the amplicon has no junction (junction_ref_start unset)."""
    if not p.junction_ref_start:
        return None
    from .clones import cluster, junction
    j0 = p.junction_ref_start - 1
    keys = {}
    for i in carriers:
        k = junction(aligned[i], j0, p.junction_length)
        if k:
            keys[aligned[i].qname] = k
    groups = sorted(Counter(cluster(keys, p.junction_max_edits).values()).values(), reverse=True)
    n = len(keys)
    return {"with_junction": n, "junctions": len(groups), "largest": groups[0] if groups else 0,
            "top_two_share": round(sum(groups[:2]) / n, 3) if n else None}


def polyclonal(jx: dict | None) -> bool:
    return bool(jx) and jx["with_junction"] >= POLY_MIN_READS and jx["junctions"] >= POLY_MIN_CLONES \
        and jx["top_two_share"] < POLY_TOP2


def _mode(p: Params) -> str:
    return "off" if not p.phase_reads else p.phase_mode


def detect(aligned, ref: str, p: Params, lo: int, hi: int):
    """The co-occurring blocks in `aligned`: a list of (block, carrier read indices,
    histogram of sites carried), or a string giving why there are none. Depends only
    on the reads and settings, not on the phase mode or the known blocks."""
    sites = candidate_sites(aligned, ref, p, lo, hi)
    if len(sites) < p.phase_min_sites:
        return f"only {len(sites)} variable sites (< {p.phase_min_sites})"
    idx = _index(aligned, {pos for pos, _ in sites})
    blocks = find_blocks(aligned, ref, p, lo, hi, idx=idx, sites=sites)
    if not blocks:
        return f"no co-occurring block of >= {p.phase_min_sites} sites"
    out = []
    for block in blocks:
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
        out.append((block, _carriers(idx, block, p), carried))
    return out


def partition(aligned, ref: str, p: Params, lo: int, hi: int,
              known_blocks: list | None = None, detected=None, keep: dict | None = None):
    """Split reads into (major, minor) populations. Returns (major, stats).

    Deliberately conservative: blocks that cannot be shown to be non-SHM are
    reported, not removed; and no more than phase_max_drop of the reads is ever
    removed. A tool that silently discarded half the data would be worse than one
    that does nothing. The cap is applied block by block, the block with the most
    sites (the surest not to be hypermutation) first: each is dropped while the
    reads removed so far, with its own, stay within the cap, and one that would
    take them over is kept and says so. Until November 2026 the cap was applied to
    all the blocks together, so one large block could keep another, however clear,
    in the library; a control carrying a 58-site paralogue-like block in a quarter
    of its reads beside a common 6-site one kept both, and its background hid every
    call in the run.

    `detected` is detect()'s result for these same reads, in this order, when it is
    already known (batch.py: a worker detects, the main process decides); `keep`, a
    dict, receives detect()'s result under "detected".
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

    found = detect(aligned, ref, p, lo, hi) if detected is None else detected
    if keep is not None:
        keep["detected"] = found
    if isinstance(found, str):
        stats["reason"] = found
        return aligned, stats
    blocks = [b for b, _c, _h in found]

    decided = []
    for block, carriers, carried in found:
        jx = junctions(aligned, carriers, p)
        if mode == "drop":
            action, why = "drop", "phase_mode=drop"
        elif mode == "auto" and matches_known(block, known_blocks, p.phase_min_sites):
            action, why = "drop", "matches a block seen in controls or another sample"
        elif mode == "auto" and polyclonal(jx):
            action, why = "drop", (f"carried by reads of {jx['junctions']} different V(D)J junctions; "
                                   "hypermutation does not repeat across B-cell clones, so a "
                                   "haplotype, paralogue or contaminant")
        elif jx and jx["with_junction"] >= POLY_MIN_READS:
            action, why = "flag", (f"its {jx['with_junction']} carriers with a readable junction fall in "
                                   f"{jx['junctions']} junction(s), the largest holding {jx['largest']}: "
                                   "an expanded clone, kept")
        else:
            action, why = "flag", ("not seen in controls or other samples; may be a clone, "
                                   "kept")
        decided.append({"block": block, "carriers": carriers, "action": action, "why": why,
                        "row": {"sites": [[pos + 1, alt] for pos, alt in sorted(block)],   # 1-based for humans
                                "n_sites": len(block), "reads": len(carriers),
                                "carried_hist": dict(sorted(carried.items())), "junctions": jx}})

    drop: set = set()
    dropped_sites = 0
    capped = []
    for d in sorted((d for d in decided if d["action"] == "drop" and d["carriers"]),
                    key=lambda d: (-len(d["block"]), -len(d["carriers"]))):
        more = drop | set(d["carriers"])
        if len(more) / len(aligned) > p.phase_max_drop:
            d["action"] = "kept (cap)"
            d["why"] += (f"; removing it as well would take more than {100 * p.phase_max_drop:.0f}% "
                         f"of the reads (it is {100 * len(d['carriers']) / len(aligned):.0f}%), "
                         "so it was kept")
            capped.append(d)
            continue
        drop, dropped_sites = more, dropped_sites + len(d["block"])
        d["kept_block"] = True
    # a read with part of a dropped block (a PCR chimera, or one too short for most of
    # its sites) can fall under phase_score_threshold and still carry ten or more of
    # its differences, more than any real read does (contaminant_min_alleles)
    if p.contaminant_min_alleles > 0:
        big = [{(pos, alt) for pos, alt in d["block"]} for d in decided
               if d.get("kept_block") and len(d["block"]) >= p.contaminant_min_alleles]
        if big:
            idx = _index(aligned, {pos for b in big for pos, _ in b})
            for i, seen in enumerate(idx):
                if i not in drop and any(sum(seen.get(pos) == alt for pos, alt in b)
                                         >= p.contaminant_min_alleles for b in big):
                    drop.add(i)
    for d in decided:
        stats["blocks"].append({**d["row"], "action": d["action"], "why": d["why"]})

    if not drop:
        if capped:
            share = len(set().union(*(d["carriers"] for d in capped))) / len(aligned)
            stats["reason"] = (f"minor population is {100 * share:.0f}% of reads "
                               f"(> {100 * p.phase_max_drop:.0f}%); refusing to phase")
        else:
            stats["reason"] = (f"{len(blocks)} co-occurring block(s) flagged, none dropped "
                               f"(phase_mode={mode})")
        return aligned, stats

    major = [r for i, r in enumerate(aligned) if i not in drop]
    stats.update(phasing_applied=True, reads_dropped_minor=len(drop),
                 haplotype_sites=dropped_sites,
                 reason=f"removed {len(drop)} reads sharing a "
                        f"{dropped_sites}-site block"
                        + (f"; {len(capped)} more block(s) kept, over the {100 * p.phase_max_drop:.0f}% cap"
                           if capped else ""))
    return major, stats
