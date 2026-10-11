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

# (pattern, index-of-mutated-base). Grouped by biological category.
MOTIFS: dict[str, tuple[str, int]] = {
    # AID hotspots (canonical WRC / complement GYW; extended RGYW / WRCY)
    "WRC": ("WRC", 2), "GYW": ("GYW", 0),
    "RGYW": ("RGYW", 2), "WRCY": ("WRCY", 2),
    "DGYW": ("DGYW", 2), "WRCH": ("WRCH", 2),
    # Pol-eta / A:T hotspots (WA top strand, TW complement)
    "WA": ("WA", 1), "TW": ("TW", 0),
    # AID coldspots
    "SYC": ("SYC", 2), "GRS": ("GRS", 0),
}

CATEGORIES = {
    "AID_hotspot": ["WRC", "GYW", "RGYW", "WRCY"],
    "polymerase_eta_hotspot": ["WA", "TW"],
    "coldspot": ["SYC", "GRS"],
}
STRAND_PAIR = ("WRC", "GYW")   # top vs bottom strand AID hotspot -> strand bias


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
    return _match(pattern, ref, pos0 - idx)


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
