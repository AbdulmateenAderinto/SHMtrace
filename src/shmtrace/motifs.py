"""SHM hot-/cold-spot motif scoring (SHMTool feature).

Motifs are defined on the reference around the mutated base (underlined). We
score (a) how many mutations fall in each motif and (b) how many motif *sites*
exist in the analysed region — the denominator for a base-composition-corrected
motif mutation frequency. Strand bias falls out of comparing a top-strand motif
(e.g. WRC, C mutated) with its bottom-strand complement (GYW, G mutated).

Degenerate IUPAC codes: W=A/T R=A/G Y=C/T S=G/C D=A/G/T K=G/T M=A/C H=A/C/T N=any
"""
from __future__ import annotations

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "W": "AT", "R": "AG", "Y": "CT", "S": "GC",
    "K": "GT", "M": "AC", "D": "AGT", "H": "ACT",
    "B": "CGT", "V": "ACG", "N": "ACGT",
}

# (pattern, index or indices of the mutated base). Grouped by biological category.
MOTIFS: dict[str, tuple[str, int | tuple[int, ...]]] = {
    # AID hotspots (canonical WRC / complement GYW; extended RGYW / WRCY and WRCH / DGYW)
    "WRC": ("WRC", 2), "GYW": ("GYW", 0),
    "RGYW": ("RGYW", 1), "WRCY": ("WRCY", 2),
    "DGYW": ("DGYW", 1), "WRCH": ("WRCH", 2),
    # WRC and GYW overlapping, so AID can deaminate both strands at once: both the G and the C
    "WGCW": ("WGCW", (1, 2)),
    # Pol-eta / A:T hotspots (WA top strand, TW complement), and the two overlapping
    "WA": ("WA", 1), "TW": ("TW", 0),
    "WATW": ("WATW", (1, 2)),
    # AID coldspots
    "SYC": ("SYC", 2), "GRS": ("GRS", 0),
}

CATEGORIES = {
    "AID_hotspot": ["WRC", "GYW", "RGYW", "WRCY"],
    "polymerase_eta_hotspot": ["WA", "TW"],
    "coldspot": ["SYC", "GRS"],
}
STRAND_PAIR = ("WRC", "GYW")   # top vs bottom strand AID hotspot -> strand bias

# The motif families reported for each call and each sample, as (key, label, bases, motifs, what
# it is). A family is a motif and its other-strand form; a call can sit in several, since WRCH is
# a WRC and WGCW is both a WRC and a GYW. The percentage is of the sample's calls at those bases.
FAMILIES: list[tuple[str, str, str, tuple[str, ...], str]] = [
    ("WRC_GYW", "WRC/GYW", "GC", ("WRC", "GYW"), "AID hotspot, the canonical form"),
    ("WRCH_DGYW", "WRCH/DGYW", "GC", ("WRCH", "DGYW"), "AID hotspot, the extended form that predicts mutability best"),
    ("WRCY_RGYW", "WRCY/RGYW", "GC", ("WRCY", "RGYW"), "AID hotspot, the classic form"),
    ("WGCW", "WGCW", "GC", ("WGCW",), "overlapping WRC/GYW: an AID hotspot on both strands"),
    ("SYC_GRS", "SYC/GRS", "GC", ("SYC", "GRS"), "AID coldspot"),
    ("WA_TW", "WA/TW", "AT", ("WA", "TW"), "Pol η hotspot: A:T mutations made in mismatch repair"),
    ("WATW", "WATW", "AT", ("WATW",), "overlapping WA/TW: a Pol η hotspot on both strands"),
]


def _match(pattern: str, ref: str, start: int) -> bool:
    if start < 0 or start + len(pattern) > len(ref):
        return False
    for k, code in enumerate(pattern):
        if ref[start + k].upper() not in IUPAC[code]:
            return False
    return True


def is_site(motif: str, ref: str, pos0: int) -> bool:
    """Is reference position pos0 (0-based) the mutated base of `motif`?"""
    pattern, idx = MOTIFS[motif]
    return any(_match(pattern, ref, pos0 - i) for i in ((idx,) if isinstance(idx, int) else idx))


def families_at(ref: str, pos0: int) -> list[str]:
    """The labels of the families whose site is reference position pos0 (0-based)."""
    return [label for _k, label, _b, ms, _w in FAMILIES if any(is_site(m, ref, pos0) for m in ms)]


def family_share(ref: str, positions0, key: str) -> tuple[int, int]:
    """(sites in the family, positions at its bases) among positions0: hits and their denominator."""
    _k, _l, bases, ms, _w = next(f for f in FAMILIES if f[0] == key)
    at = [p for p in positions0 if 0 <= p < len(ref) and ref[p].upper() in bases]
    return sum(1 for p in at if any(is_site(m, ref, p) for m in ms)), len(at)


def background(ref: str, lo: int, hi: int) -> dict[str, float | None]:
    """Each family's share of the window's bases it can sit on: what a call would hit by chance."""
    out = {}
    for key, *_ in FAMILIES:
        hits, n = family_share(ref, range(lo, hi + 1), key)
        out[key] = round(100 * hits / n, 1) if n else None
    return out


def site_counts(ref: str, positions: list[int]) -> dict[str, int]:
    """Number of positions in `positions` (0-based) that are each motif's site."""
    counts = {m: 0 for m in MOTIFS}
    for p in positions:
        for m in MOTIFS:
            if is_site(m, ref, p):
                counts[m] += 1
    return counts


def score_variant_positions(ref: str, positions0: list[int]) -> dict[str, int]:
    """How many of the given mutated positions fall in each motif."""
    hits = {m: 0 for m in MOTIFS}
    for p in positions0:
        for m in MOTIFS:
            if is_site(m, ref, p):
                hits[m] += 1
    return hits
