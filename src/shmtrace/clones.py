"""Count B-cell clones, not reads, using the V(D)J junction as a barcode.

Without UMIs every read counts once, so a clone that expanded in the germinal
centre, or one template that PCR copied many times, counts as often as it was
read. The review of this project called that its main unresolved problem.

An amplicon that starts in the V region carries a natural barcode: the bases
just 5' of the J segment are the V(D)J junction, made by recombination and
unique to the clone. Reads are grouped by it, each group's shared mutations are
taken as the clone's, and the clone counts once.

This gives counts per clone, not per molecule: two cells of one clone, or two
PCR copies of one template, look the same. It is the unit that matters for the
expansion problem, and it needs no change to the protocol.

What tells a clone apart from a coincidence (validated on J꜀4,
docs/calibration/CLONES.md): reads sharing a junction share their mutations far
more than reads that do not. A junction group whose reads share nothing is
either an unmutated clone or not a clone at all, and the figure shows which.

The junction is read from the unaligned 5' sequence: for J꜀4 the germline
reference upstream of the J segment is not what a rearranged read carries, so
the aligner starts the alignment at or just after the J segment, and the bases
before that are the junction. A read whose alignment starts far from the J
segment (an unrearranged allele, or a read that does not reach it) has no
junction and is left out of the clone count.
"""
from __future__ import annotations

from collections import Counter, defaultdict

SLACK = 18      # bp either side of the J segment start that an alignment may begin


def _edit(a: str, b: str, cap: int) -> int:
    """Edit distance with the 5' edge free, capped at `cap` + 1.

    The window's 5' edge is an arbitrary cut 20 bases from J, so an indel
    inside it slides one extra base in at that edge. Plain Levenshtein charges
    that twice (the indel and the slid base); leading gaps cost nothing here,
    so one indel costs one. On simulated clones with ten times the real indel
    rate this lifts recovery from adjusted Rand 0.96 to 0.99, and it changes
    nothing at the real rate (docs/calibration/LINEAGE_SIM.md).
    """
    prev = [0] * (len(b) + 1)
    for ca in a:
        cur = [0] * (len(b) + 1)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return min(prev[-1], cap + 1)


def junction(read, j_start0: int, length: int) -> str | None:
    """The `length` bases 5' of where the read's alignment begins, if that is at J."""
    if not (j_start0 - SLACK <= read.r_st <= j_start0 + SLACK):
        return None
    q = getattr(read, "q_st", 0)
    if q < length:
        return None
    return read.seq[q - length:q]


def cluster(keys: dict, max_edits: int) -> dict:
    """read name -> clone id. Greedy by abundance against clone centroids.

    The commonest junction founds the first clone; each junction joins the
    first centroid within `max_edits` or founds its own. Centroids, not
    chains, so a run of near neighbours cannot merge two different clones.
    """
    counts = Counter(keys.values())
    centroids: list[str] = []
    to_clone: dict = {}
    for k, _n in counts.most_common():
        hit = next((i for i, c in enumerate(centroids) if _edit(k, c, max_edits) <= max_edits),
                   None)
        if hit is None:
            centroids.append(k)
            hit = len(centroids) - 1
        to_clone[k] = hit
    return {name: to_clone[k] for name, k in keys.items()}


def call_clones(aligned, variants, p, lo: int, hi: int) -> dict | None:
    """Group a sample's reads into clones and count each clone's mutations once.

    A clone's mutations are the called mutations carried by at least half of
    its reads that cover the position. A one-read clone keeps all of its own.
    """
    if not p.junction_ref_start:
        return None
    j0 = p.junction_ref_start - 1
    keys = {}
    for r in aligned:
        k = junction(r, j0, p.junction_length)
        if k:
            keys[r.qname] = k
    if not keys:
        return None
    clone_of = cluster(keys, p.junction_max_edits)
    called = {(v.pos - 1, v.alt) for v in variants}
    reads = {r.qname: r for r in aligned if r.qname in clone_of}
    members = defaultdict(list)
    for name, c in clone_of.items():
        members[c].append(name)
    bp = hi - lo + 1
    rows = []
    for c, names in members.items():
        cover = Counter()
        carry = Counter()
        for n in names:
            for pos, base, _q in reads[n].pairs:
                if lo <= pos <= hi:
                    cover[pos] += 1
                    if (pos, base) in called:
                        carry[(pos, base)] += 1
        need = 1 if len(names) == 1 else None
        clonal = [m for m, k in carry.items()
                  if k >= (need or max(1, (cover[m[0]] + 1) // 2))]
        rows.append({"clone": c, "reads": len(names),
                     "junction": keys[names[0]],
                     "mutations": len(clonal),
                     "positions": sorted(pos + 1 for pos, _a in clonal)})
    rows.sort(key=lambda r: -r["reads"])
    n_reads = sum(r["reads"] for r in rows)
    weighted = sum(r["mutations"] for r in rows) / (len(rows) * bp)
    return {
        "reads_with_junction": n_reads,
        "reads_analysed": len(aligned),
        "clones": len(rows),
        "largest_clone_pct": round(100 * rows[0]["reads"] / n_reads, 1),
        # every clone weighted once: robust to one clone's expansion
        "clone_weighted_frequency_x1e3": round(1e3 * weighted, 4),
        # distinct mutations over clones, a mutation counting once per clone
        "clonal_mutation_events": sum(r["mutations"] for r in rows),
        "sharing": sharing(reads, clone_of, called),
        "rows": rows,
        # read name -> clone; the caller keeps it with the reads (airr.py), not in the summary
        "clone_of": clone_of,
    }


def sharing(reads: dict, clone_of: dict, called: set, max_pairs: int = 40000) -> dict:
    """Mutation sharing (Jaccard) between reads of one clone and of different clones.

    The check that the junction is a barcode at all. Within a real clone, reads
    inherit the clone's mutations and share them; across clones they should
    share almost none. Pairs where neither read carries a mutation are skipped.
    """
    import itertools
    import random
    muts = {n: {(p, b) for p, b, _q in r.pairs if (p, b) in called} for n, r in reads.items()}
    names = list(reads)
    pairs = list(itertools.combinations(names, 2))
    random.Random(0).shuffle(pairs)
    within, between = [], []
    for a, b in pairs[:max_pairs]:
        A, B = muts[a], muts[b]
        if not A and not B:
            continue
        j = len(A & B) / len(A | B)
        (within if clone_of[a] == clone_of[b] else between).append(j)
    mean = lambda xs: round(sum(xs) / len(xs), 4) if xs else None   # noqa: E731
    return {"within": mean(within), "between": mean(between),
            "pairs_within": len(within), "pairs_between": len(between)}
