"""Which repair routes made the mutations? Spectrum indices and decomposition.

AID deaminates C to U. What happens next leaves a signature in the spectrum:

  replication over U                    G:C transitions (C>T, G>A)
  UNG excision, then TLS over the       G:C transversions (REV1 inserts C: C>G, G>C)
    abasic site                           and some transitions
  MSH2/MSH6-EXO1 patch, Pol-eta         A:T mutations, more at WA than at TW
    (non-canonical mismatch repair)       (Pol-eta errs on one strand)

In Ung-/- mice almost every G:C mutation is a transition; in Msh2-/- mice A:T
mutations nearly vanish; in Ung-/- Msh2-/- mice only G:C transitions remain (Rada et
al. 2002, 2004). The indices below read these routes off a spectrum, per animal, so
groups can be compared with the animal as the unit (replicates.permutation_p and Welch).

Indices (unique mutations; per animal and per group):
  pct_AT              A:T mutations / all                    MSH2-Pol-eta
  pct_GC_Ts_of_GC     G:C transitions / G:C mutations         replication over U
  pct_GC_Tv_of_GC     G:C transversions / G:C mutations       UNG-dependent TLS
  pct_CG_of_GC_Tv     C>G + G>C / G:C transversions           REV1 dC insertion
  A_over_T            A mutations / T mutations (top strand)  Pol-eta strand bias
  WA_over_TW          mutations at WA / at TW
  WRC_over_GYW        mutations at WRC / at GYW               AID strand

Decomposition: a group's 12-type spectrum as a non-negative mixture of reference
spectra (groups of the same project, e.g. Ung-/-, Msh2-/-, double knockouts, or any
spectra given), fitted by non-negative least squares on proportions, with 95%
intervals from resampling animals and their mutations.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
from scipy import optimize, stats

from . import motifs
from .replicates import permutation_p

SUBS12 = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
TS = {"C>T", "G>A", "A>G", "T>C"}


def indices(calls: list[tuple[int, str, str]], ref: str) -> dict:
    """calls: (1-based position, ref base, alt base)."""
    c = Counter(f"{r}>{a}" for _p, r, a in calls)
    n = sum(c.values())
    gc = sum(v for k, v in c.items() if k[0] in "GC")
    at = n - gc
    gc_ts = c["C>T"] + c["G>A"]
    gc_tv = gc - gc_ts
    a_m = sum(v for k, v in c.items() if k[0] == "A")
    t_m = sum(v for k, v in c.items() if k[0] == "T")
    wa = sum(1 for p, r, _a in calls if r == "A" and motifs.is_site("WA", ref, p - 1))
    tw = sum(1 for p, r, _a in calls if r == "T" and motifs.is_site("TW", ref, p - 1))
    wrc = sum(1 for p, r, _a in calls if r == "C" and motifs.is_site("WRC", ref, p - 1))
    gyw = sum(1 for p, r, _a in calls if r == "G" and motifs.is_site("GYW", ref, p - 1))

    def pct(a, b):
        return round(100 * a / b, 2) if b else None

    def ratio(a, b):
        return round(a / b, 3) if b else None
    return {"mutations": n, "pct_AT": pct(at, n), "pct_GC_Ts_of_GC": pct(gc_ts, gc),
            "pct_GC_Tv_of_GC": pct(gc_tv, gc), "pct_CG_of_GC_Tv": pct(c["C>G"] + c["G>C"], gc_tv),
            "A_over_T": ratio(a_m, t_m), "WA_over_TW": ratio(wa, tw),
            "WRC_over_GYW": ratio(wrc, gyw),
            "_counts": {"at": at, "gc": gc, "gc_ts": gc_ts, "gc_tv": gc_tv}}


def proportions(spec: dict) -> np.ndarray:
    v = np.array([spec.get(k, 0) for k in SUBS12], float)
    return v / v.sum() if v.sum() else v


def decompose(spec: dict, references: dict[str, dict]) -> dict:
    """Non-negative mixture weights (summing to 1) of reference spectra for `spec`."""
    names = list(references)
    A = np.column_stack([proportions(references[n]) for n in names])
    b = proportions(spec)
    # NNLS with the weights summing to 1: append a heavily weighted row of ones
    w_row = 100.0
    A2 = np.vstack([A, w_row * np.ones(len(names))])
    b2 = np.concatenate([b, [w_row]])
    w, _res = optimize.nnls(A2, b2)
    w = w / w.sum() if w.sum() else w
    fitted = A @ w
    ss_res = float(np.sum((b - fitted) ** 2))
    ss_tot = float(np.sum((b - b.mean()) ** 2)) or 1.0
    return {"weights": dict(zip(names, (float(x) for x in w), strict=True)),
            "r2": 1 - ss_res / ss_tot, "residual": dict(zip(SUBS12, (float(x) for x in b - fitted),
                                                            strict=True))}


def decompose_animals(animals: list[dict], references: dict[str, list[dict]], boot: int = 999,
                      seed: int = 1) -> dict:
    """Decompose the pooled spectrum of `animals` (each a 12-type count dict), with 95%
    intervals from resampling animals, and mutations within each animal, in both the
    sample and the references (whose spectra are themselves estimates)."""
    rng = np.random.default_rng(seed)

    def pooled(xs):
        out = Counter()
        for x in xs:
            out.update(x)
        return dict(out)

    def resample(xs):
        pick = [xs[i] for i in rng.integers(0, len(xs), len(xs))]
        out = Counter()
        for x in pick:
            v = np.array([x.get(k, 0) for k in SUBS12], float)
            if v.sum():
                draw = rng.multinomial(int(v.sum()), v / v.sum())
                out.update(dict(zip(SUBS12, draw.tolist(), strict=True)))
        return dict(out)
    point = decompose(pooled(animals), {k: pooled(v) for k, v in references.items()})
    reps = {k: [] for k in references}
    for _ in range(boot):
        d = decompose(resample(animals), {k: resample(v) for k, v in references.items()})
        for k, v in d["weights"].items():
            reps[k].append(v)
    point["ci"] = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
                   for k, v in reps.items()}
    return point


def compare_indices(res: dict, case: str, control: str) -> list[dict]:
    """Each index per animal of two groups of a project result, compared animal-wise."""
    ref = res["reference"]
    by_animal: dict = {}
    for v in res["variants"]:
        by_animal.setdefault(v["sample"], []).append((v["pos"], v["ref"], v["alt"]))
    groups = {s["sample"]: s["group"] for s in res["samples"]}
    idx = {a: indices(calls, ref) for a, calls in by_animal.items()}
    rows = []
    for key in ("pct_AT", "pct_GC_Ts_of_GC", "pct_GC_Tv_of_GC", "pct_CG_of_GC_Tv", "A_over_T",
                "WA_over_TW", "WRC_over_GYW"):
        a = [idx[s][key] for s, g in groups.items() if g == case and s in idx
             and idx[s][key] is not None]
        b = [idx[s][key] for s, g in groups.items() if g == control and s in idx
             and idx[s][key] is not None]
        row = {"index": key, "group": case, "control": control,
               "mean_group": round(float(np.mean(a)), 3) if a else None,
               "mean_control": round(float(np.mean(b)), 3) if b else None,
               "animals": f"{len(a)} vs {len(b)}"}
        if len(a) >= 2 and len(b) >= 2:
            row["p_welch"] = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
            row["p_permutation"] = permutation_p(a, b)
        rows.append(row)
    # pooled counts for the three route shares, Fisher, to show pseudoreplication
    def pooled(g):
        tot = Counter()
        for s, gg in groups.items():
            if gg == g and s in idx:
                tot.update(idx[s]["_counts"])
        return tot
    pa, pb = pooled(case), pooled(control)
    for key, (num, den) in {"pct_AT": ("at", None), "pct_GC_Tv_of_GC": ("gc_tv", "gc")}.items():
        na, nb = pa[num], pb[num]
        da = (pa["at"] + pa["gc"]) if den is None else pa[den]
        db = (pb["at"] + pb["gc"]) if den is None else pb[den]
        if da and db:
            p = stats.fisher_exact([[na, da - na], [nb, db - nb]])[1]
            next(r for r in rows if r["index"] == key)["p_pooled_mutations"] = float(p)
    return rows



def group_spectra(res: dict, group: str) -> list[dict]:
    """Per-animal 12-type spectra (unique mutations) for a group of a project result."""
    members = {s["sample"] for s in res["samples"] if s["group"] == group}
    out: dict = {m: Counter() for m in members}
    for v in res["variants"]:
        if v["sample"] in members:
            out[v["sample"]][f"{v['ref']}>{v['alt']}"] += 1
    return [dict(c) for c in out.values()]


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="shmtrace-repair", description=__doc__.split("\n\n")[0])
    ap.add_argument("result", help="project.json")
    ap.add_argument("--case", required=True, help="group to describe")
    ap.add_argument("--control", help="group to compare it with, animal by animal")
    ap.add_argument("--reference", action="append", default=[],
                    help="NAME=GROUP: a group whose spectrum is a reference route "
                         "(e.g. ung_msh2=UngMsh2KO); repeat for each")
    a = ap.parse_args(argv)
    res = json.load(open(a.result))
    if a.control:
        for r in compare_indices(res, a.case, a.control):
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    if a.reference:
        refs = {n: group_spectra(res, g) for n, g in (x.split("=", 1) for x in a.reference)}
        d = decompose_animals(group_spectra(res, a.case), refs)
        print(json.dumps({k: d[k] for k in ("weights", "ci", "r2")}, indent=1))


if __name__ == "__main__":
    main()
