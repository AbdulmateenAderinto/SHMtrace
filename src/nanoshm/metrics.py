"""Downstream SHM summaries: unique-variant density, spectrum, call clustering."""
from __future__ import annotations

from bisect import bisect_left
from math import exp

from .config import Params
from .variants import BASES, Variant


def unique_variant_density(variants: list[Variant], analysed_bp: int,
                           mean_depth: float, p: Params) -> float:
    """Unique variants / (analysed_bp x depth) x scale  (matches the paper)."""
    denom = analysed_bp * max(1.0, mean_depth)
    if denom == 0:
        return 0.0
    return len(variants) / denom * p.density_scale


def spectrum_matrix(variants: list[Variant], weighted: bool = False) -> dict:
    """12-type substitution matrix. weighted=False counts each unique variant
    once; weighted=True counts each by its supporting-read number."""
    mat = {r: {a: 0 for a in BASES} for r in BASES}
    for v in variants:
        if v.ref in BASES and v.alt in BASES:
            mat[v.ref][v.alt] += v.supporting_reads if weighted else 1
    return mat


def gc_at_summary(variants: list[Variant], weighted: bool = False) -> dict:
    """Fraction of substitutions occurring at G:C vs A:T reference bases —
    the hallmark AID readout."""
    gc = at = 0
    for v in variants:
        w = v.supporting_reads if weighted else 1
        if v.ref in ("G", "C"):
            gc += w
        elif v.ref in ("A", "T"):
            at += w
    total = gc + at
    return {"at_GC": gc, "at_AT": at,
            "fraction_GC": round(gc / total, 4) if total else 0.0}


def positional_clusters(variants: list[Variant], analysed_bp: int, window: int = 10,
                        alpha: float = 1e-3) -> list[dict]:
    """Windows holding far more calls than a random scatter would give.

    A localised PCR or sequencing artefact puts several calls into a few bases
    of one library. SHM does not: it is spread over the amplicon, and a real
    co-occurring block (a second haplotype or an expanded clone) is caught by
    read phasing instead, because its calls share reads.

    For each window of `window` bases holding k calls, the Poisson tail
    P(X >= k) at lambda = (calls x window / analysed_bp), multiplied by the
    number of windows, must be below `alpha`. Only the strongest overlapping
    window of a run is reported.
    """
    if analysed_bp <= 0 or len(variants) < 3:
        return []
    # every call, not every position: two alternative bases called at one
    # position is itself a sign of an artefact rather than one mutation
    pos = sorted(v.pos for v in variants)
    n_tests = max(1, analysed_bp - window + 1)
    lam = len(variants) * window / analysed_bp
    hits = []
    for i, start in enumerate(pos):
        if i and pos[i - 1] == start:
            continue
        j = bisect_left(pos, start + window, i)     # pos is sorted: no linear scan
        k = j - i
        if k < 3:
            continue
        members = pos[i:j]
        # P(X >= k) for X ~ Poisson(lam)
        term, cdf = exp(-lam), exp(-lam)
        for j in range(1, k):
            term *= lam / j
            cdf += term
        p_adj = min(1.0, max(0.0, 1.0 - cdf) * n_tests)
        if p_adj < alpha:
            hits.append({"start": members[0], "end": members[-1], "calls": k,
                         "p_adjusted": p_adj})
    out: list[dict] = []
    for h in sorted(hits, key=lambda h: (h["p_adjusted"], -h["calls"])):
        if all(h["end"] < o["start"] or h["start"] > o["end"] for o in out):
            out.append(h)
    return sorted(out, key=lambda h: h["start"])


def spectrum_note(variants: list[Variant], expected: tuple[str, ...]) -> dict:
    """Share of calls in the substitution classes the experiment can produce.

    `expected` is a tuple of "C>T"-style classes. In a Ung-/-Msh2-/- background
    AID can give only C>T and G>A, so a low share is evidence the calls are not
    AID's, whatever their frequency (docs/calibration/ECTOPIC_AID.md).
    """
    total = sum(1 for v in variants if v.ref in BASES and v.alt in BASES)
    if not total:
        return {"calls": 0, "in_expected": 0, "fraction": None}
    k = sum(1 for v in variants
            if v.ref in BASES and v.alt in BASES and f"{v.ref}>{v.alt}" in expected)
    return {"calls": total, "in_expected": k, "fraction": round(k / total, 4)}


__all__ = ["unique_variant_density", "spectrum_matrix", "gc_at_summary",
           "positional_clusters", "spectrum_note"]
