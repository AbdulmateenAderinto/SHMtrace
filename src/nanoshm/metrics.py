"""Downstream SHM summaries: unique-variant density and substitution spectrum."""
from __future__ import annotations

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
