"""UMI consensus: one read per starting molecule.

A unique molecular identifier (UMI) is a random tag added to each template molecule
before PCR, so every read sharing a UMI is a copy of one molecule. Building one
consensus per UMI does two things for SHM:

  * sequencing errors, which differ between copies, are outvoted, so the background
    falls and single-molecule mutations become callable;
  * PCR duplicates and jackpots collapse to one molecule, so frequencies count
    molecules rather than amplification.

Where the UMI is (Params.umi):

  "flank"   next to an anchor at a read end. The anchor is a known sequence just
            inside the UMI, by default the first or last `umi_anchor_length` bases
            of the reference (the primer); the UMI is the `umi_length` bases just
            outside it, at `umi_offset` bases' distance (for a spacer). Found with
            up to `umi_anchor_edits` edits. `umi_ends`: "left", "right" or "both"
            (dual UMIs, joined into one key).
  "header"  in the read name, as captured by `umi_header_regex` (group 1), e.g.
            umi_tools' "readname_ACGTACGT" or "...;UMI=ACGT...".

Reads are grouped by UMI with errors allowed (`umi_max_edits` per UMI end), greedily
by abundance against group centroids, as clones.cluster does; a UMI seen once or twice
may also join a group three or more times its size within two more edits (UMI-tools'
"directional" idea: it is far likelier an error copy than a molecule of its own). Groups with fewer than
`umi_min_reads` reads are dropped: a consensus of one or two reads corrects little.

The consensus is built column by column on the reference: at each position, the
base (or deletion) with the highest summed log-likelihood over the group's reads,
each read's base weighted by its own quality. Its quality is the posterior against
the runner-up, capped at Q60. Where reads disagree the quality is low, so a
disputed base is not trusted as a mutation downstream.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

from .align import AlignedRead

QMAX = 60


IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def _semiglobal(text: str, pat: str) -> tuple[int, int, int]:
    """Best match of all of `pat` within `text`: (edits, start, end), end exclusive.

    `pat` may carry IUPAC codes (degenerate primers): R matches A or G, and so on."""
    n, m = len(text), len(pat)
    if not n or not m:
        return m, 0, 0
    # rows: pattern, columns: text; free start and end in text
    prev = [0] * (n + 1)
    starts_prev = list(range(n + 1))
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        starts = [0] * (n + 1)
        ok = IUPAC.get(pat[i - 1], pat[i - 1])
        for j in range(1, n + 1):
            d = prev[j - 1] + (text[j - 1] not in ok)
            u = prev[j] + 1
            left = cur[j - 1] + 1
            if d <= u and d <= left:
                cur[j], starts[j] = d, starts_prev[j - 1]
            elif u <= left:
                cur[j], starts[j] = u, starts_prev[j]
            else:
                cur[j], starts[j] = left, starts[j - 1]
        prev, starts_prev = cur, starts
    best = min(range(1, n + 1), key=lambda j: (prev[j], -j))
    return prev[best], starts_prev[best], best


def extract(read: AlignedRead, ref: str, p) -> str | None:
    """The read's UMI key, or None if it cannot be read."""
    if p.umi == "header":
        m = re.search(p.umi_header_regex, read.qname)
        return m.group(1).upper() if m else None
    L, A = p.umi_length, p.umi_anchor_length
    seq = read.seq
    parts = []
    if p.umi_ends in ("left", "both"):
        anchor = (p.umi_left_anchor or ref[:A]).upper()
        # the anchor sits near where the reference starts in the read
        approx = read.q_st - read.r_st
        lo = max(0, approx - 40)
        hi = min(len(seq), approx + len(anchor) + 40)
        ed, s, _e = _semiglobal(seq[lo:hi], anchor)
        s += lo
        end = s - p.umi_offset
        start = end - L
        # a UMI at the very end of a read may have lost a base or two: keep what is
        # there if it is short by no more than the grouping tolerates
        if ed > p.umi_anchor_edits or end - max(start, 0) < L - p.umi_max_edits:
            return None
        parts.append(seq[max(start, 0):end])
    if p.umi_ends in ("right", "both"):
        anchor = (p.umi_right_anchor or ref[-A:]).upper()
        q_en = read.q_en if read.q_en >= 0 else read.q_st + len(read.pairs)
        approx = q_en + (len(ref) - read.r_en) - len(anchor)
        lo = max(0, approx - 40)
        hi = min(len(seq), approx + len(anchor) + 40)
        ed, _s, e = _semiglobal(seq[lo:hi], anchor)
        e += lo
        start = e + p.umi_offset
        if ed > p.umi_anchor_edits or min(start + L, len(seq)) - start < L - p.umi_max_edits:
            return None
        parts.append(seq[start:start + L])
    return "".join(parts).upper() if parts else None


def _edit(a: str, b: str, cap: int) -> int:
    """Levenshtein distance, stopping once it exceeds `cap` (returns cap + 1)."""
    if abs(len(a) - len(b)) > cap:
        return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        lo = cap + 1
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
            lo = min(lo, cur[j])
        if lo > cap:
            return cap + 1
        prev = cur
    return prev[-1]


def group(keys: dict[str, str], max_edits: int) -> dict[str, int]:
    """read name -> group id: greedy by abundance against centroids.

    Candidates are found through exact matches of UMI segments: with at most
    `max_edits` edits, one of `max_edits + 1` segments of the UMI is untouched
    (pigeonhole), so only centroids sharing a segment need an edit distance.
    """
    counts = Counter(keys.values())
    nseg = max_edits + 1
    centroids: list[str] = []
    index: dict = defaultdict(set)
    to_group: dict[str, int] = {}

    def segments(u: str):
        k = max(1, len(u) // nseg)
        return [(i, u[i * k:(i + 1) * k if i < nseg - 1 else len(u)]) for i in range(nseg)]

    size: list[int] = []
    wide = max_edits + 2
    for u, n in counts.most_common():
        cand = set()
        for i, s in segments(u):
            cand |= index.get((i, s), set())
        hit = None
        for c in sorted(cand):
            if _edit(u, centroids[c], max_edits) <= max_edits:
                hit = c
                break
        if hit is None and n <= 2:
            # directional, as in UMI-tools: a UMI seen once or twice, a little beyond
            # the tolerance of a group at least three times its size, is far likelier
            # an error copy of that molecule than a molecule of its own
            for c in range(len(centroids)):
                if size[c] >= 3 * n and _edit(u, centroids[c], wide) <= wide:
                    hit = c
                    break
        if hit is None:
            hit = len(centroids)
            centroids.append(u)
            size.append(0)
            for i, s in segments(u):
                index[(i, s)].add(hit)
        size[hit] += n
        to_group[u] = hit
    return {name: to_group[u] for name, u in keys.items()}


def _consensus(members: list[AlignedRead], ref: str, name: str) -> AlignedRead | None:
    cols: dict[int, list] = defaultdict(list)
    for r in members:
        got = {rp: (b, q) for rp, b, q in r.pairs}
        for rp in range(r.r_st, r.r_en):
            cols[rp].append(got.get(rp, ("-", 20)))      # covered but not aligned: deleted
    pairs, deleted = [], []
    for rp in sorted(cols):
        obs = cols[rp]
        if len(obs) < max(2, len(members) // 2):
            continue                                      # too few reads reach here
        ll = {}
        for cand in "ACGT-":
            s = 0.0
            for b, q in obs:
                e = min(0.75, 10 ** (-max(q, 1) / 10))
                s += math.log(1 - e) if b == cand else math.log(e / 4)
            ll[cand] = s
        order = sorted(ll, key=ll.get, reverse=True)
        best, second = order[0], order[1]
        diff = ll[best] - ll[second]                      # natural log odds
        q = QMAX if diff > QMAX / 4.343 else int(round(10 * diff / math.log(10)))
        q = max(0, min(QMAX, q))
        if best == "-":
            deleted.append(rp)
            continue
        pairs.append((rp, best, q))
    if not pairs:
        return None
    rep = max(members, key=lambda r: sum(q for _r, _b, q in r.pairs) / max(1, len(r.pairs)))
    strand = Counter(r.strand for r in members).most_common(1)[0][0]
    indels = [(rp, "D", 1, 30) for rp in deleted]
    return AlignedRead(name, rep.ctg, pairs[0][0], pairs[-1][0] + 1, pairs, rep.seq, strand,
                       indels, rep.q_st, rep.q_en)


def consensus_reads(aligned: list[AlignedRead], ref: str, p,
                    fates: dict | None = None) -> tuple[list[AlignedRead], dict]:
    """Group aligned reads by UMI and return one substitution consensus per molecule.

    Insertion voting and contiguous deletion reconstruction are not implemented;
    the resulting reads must not be used for quantitative indel analysis.
    """
    if not p.exclude_indels:
        raise ValueError("Indel analysis after UMI consensus is not supported. "
                         "Set exclude_indels=True to analyze substitutions.")
    keys: dict[str, str] = {}
    no_umi = 0
    for r in aligned:
        u = extract(r, ref, p)
        if u and "N" not in u:
            keys[r.qname] = u
        else:
            no_umi += 1
            if fates is not None:
                fates[r.qname] = {"fate": "no UMI"}
    ends = 2 if (p.umi == "flank" and p.umi_ends == "both") else 1
    gid = group(keys, p.umi_max_edits * ends)
    by_group: dict[int, list] = defaultdict(list)
    byname = {r.qname: r for r in aligned}
    for name, g in gid.items():
        by_group[g].append(byname[name])
    out, small = [], 0
    sizes = []
    for g, members in sorted(by_group.items()):
        sizes.append(len(members))
        if len(members) < p.umi_min_reads:
            small += 1
            if fates is not None:
                for r in members:
                    fates[r.qname] = {"fate": "UMI group too small", "reads": len(members)}
            continue
        c = _consensus(members, ref, f"umi{g}_{len(members)}reads")
        if c is not None:
            c.members = [r.qname for r in members]     # for scoring simulations
            out.append(c)
    sizes.sort()
    stats = {"umi_reads_in": len(aligned), "umi_reads_without_umi": no_umi,
             "umi_groups": len(by_group), "umi_groups_too_small": small,
             "umi_molecules": len(out),
             "umi_median_reads_per_group": sizes[len(sizes) // 2] if sizes else 0}
    return out, stats
