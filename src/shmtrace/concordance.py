"""Do two measurements of the same animals agree? (e.g. Sanger clones and ONT reads)

The comparison the manuscript makes between platforms, as a tested function.
Each animal is measured twice (A and B); animals are the replicates, so every
paired test uses the animal as the matching factor.

  spectrum     the twelve substitution types: pooled proportions with 95%
               Wilson intervals; Pearson r (Fisher 95% CI, over the twelve
               types) and the least-squares slope of B on A; mean and largest
               absolute difference in percentage points.
  paired       with two or more animals: the platform-by-type interaction of a
               two-way repeated-measures ANOVA, paired t per type with a Sidak
               adjustment over the twelve, and the same for G:C transitions
               and transversions (adjusted over the two).
  positional   share of calls along the window in fixed bins, A against B.
  hotspots     share of G:C calls at WRC/GYW and of A:T calls at WA/TW, per
               animal, against the share of such sites in the window.
  density      each animal's unique-mutation density on both platforms.

Caveats travel with the result (`caveats`): the twelve proportions are a
closed composition, so the correlation describes agreement rather than testing
it, and with three animals a non-significant paired test does not show the
platforms are equal.
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np
from scipy import stats

from . import motifs

SUBS12 = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
GC_TS = ("G>A", "C>T")
GC_TV = ("G>T", "G>C", "C>A", "C>G")

CAVEATS = [
    "The twelve proportions are a closed composition and not independent; Pearson r over them "
    "describes agreement, and its interval is over types, not animals.",
    "With few animals, non-significant paired tests do not show that the two measurements are "
    "equal.",
]


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * max(0.0, c - h), 100 * min(1.0, c + h))


def _props(c: dict) -> np.ndarray:
    tot = sum(c.get(k, 0) for k in SUBS12)
    return np.array([100 * c.get(k, 0) / tot if tot else np.nan for k in SUBS12])


def _sidak(p: float, m: int) -> float:
    return 1 - (1 - p) ** m


def rm_interaction(a: np.ndarray, b: np.ndarray) -> tuple[float, int, int, float]:
    """Platform x type interaction, two-way repeated measures (subjects x 2 x types)."""
    y = np.stack([a, b], axis=1)                  # subject, platform, type
    s, _p, t = y.shape
    g = y.mean()
    m_pt, m_sp, m_st = y.mean(0), y.mean(2), y.mean(1)
    m_p, m_t, m_s = y.mean((0, 2)), y.mean((0, 1)), y.mean((1, 2))
    ss_pt = s * ((m_pt - m_p[:, None] - m_t[None, :] + g) ** 2).sum()
    res = (y - m_pt[None] - m_sp[:, :, None] - m_st[:, None, :]
           + m_p[None, :, None] + m_t[None, None, :] + m_s[:, None, None] - g)
    df1, df2 = t - 1, (s - 1) * (t - 1)
    ss_err = (res ** 2).sum()
    if df2 == 0 or ss_err == 0:
        return float("nan"), df1, df2, float("nan")
    f = (ss_pt / df1) / (ss_err / df2)
    return float(f), df1, df2, float(stats.f.sf(f, df1, df2))


def spectrum(spec_a: dict[str, dict], spec_b: dict[str, dict]) -> dict:
    """spec_x: animal -> {"A>G": count, ...}; the same animals in both."""
    animals = [k for k in spec_a if k in spec_b]
    pool_a, pool_b = Counter(), Counter()
    for k in animals:
        pool_a.update(spec_a[k])
        pool_b.update(spec_b[k])
    na = sum(pool_a.get(k, 0) for k in SUBS12)
    nb = sum(pool_b.get(k, 0) for k in SUBS12)
    pa, pb = _props(pool_a), _props(pool_b)
    r, p = stats.pearsonr(pa, pb)
    z, se = math.atanh(min(r, 0.999999)), 1 / math.sqrt(len(SUBS12) - 3)
    out = {
        "animals": animals, "calls_a": na, "calls_b": nb,
        "types": [{"type": k, "pct_a": round(float(x), 3), "pct_b": round(float(y), 3),
                   "ci_a": [round(v, 3) for v in wilson(pool_a.get(k, 0), na)],
                   "ci_b": [round(v, 3) for v in wilson(pool_b.get(k, 0), nb)]}
                  for k, x, y in zip(SUBS12, pa, pb, strict=True)],
        "pearson_r": float(r), "pearson_p": float(p),
        "pearson_ci": [math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)],
        "slope": float(np.polyfit(pa, pb, 1)[0]),
        "mean_abs_diff_pp": float(np.mean(np.abs(pa - pb))),
        "max_abs_diff_pp": float(np.max(np.abs(pa - pb))),
    }
    usable = [k for k in animals if sum(spec_a[k].get(t, 0) for t in SUBS12)
              and sum(spec_b[k].get(t, 0) for t in SUBS12)]
    if len(usable) >= 2:
        A = np.array([_props(spec_a[k]) for k in usable])
        B = np.array([_props(spec_b[k]) for k in usable])
        f, d1, d2, pf = rm_interaction(A, B)
        per = []
        for i, k in enumerate(SUBS12):
            d = A[:, i] - B[:, i]
            pt = float(stats.ttest_rel(A[:, i], B[:, i]).pvalue) if np.any(d != d[0]) else 1.0
            per.append({"type": k, "mean_diff_pp": float(d.mean()), "p": pt,
                        "p_sidak": _sidak(pt, len(SUBS12))})

        def share(M, keys):
            return np.array([sum(M[j, SUBS12.index(x)] for x in keys) for j in range(len(M))])
        tstv = {}
        for lab, keys in (("gc_transitions", GC_TS), ("gc_transversions", GC_TV)):
            sa, sb = share(A, keys), share(B, keys)
            pt = float(stats.ttest_rel(sa, sb).pvalue) if np.any(sa - sb != (sa - sb)[0]) else 1.0
            tstv[lab] = {"pct_a": [float(x) for x in sa], "pct_b": [float(x) for x in sb],
                         "p": pt, "p_sidak": _sidak(pt, 2)}
        out["paired"] = {"animals": usable, "interaction_F": f, "interaction_df": [d1, d2],
                         "interaction_p": pf, "per_type": per,
                         "min_p_sidak": min(x["p_sidak"] for x in per), "ts_tv": tstv}
    return out


def positional(pos_a: list[int], pos_b: list[int], lo: int, hi: int, bin_bp: int = 25) -> dict:
    """pos_x: 1-based positions of unique calls; lo, hi: 1-based window."""
    edges = np.arange(lo, hi + 2, bin_bp)
    if edges[-1] < hi + 1:
        edges = np.append(edges, hi + 1)
    ha = np.histogram(pos_a, edges)[0]
    hb = np.histogram(pos_b, edges)[0]
    fa = 100 * ha / max(ha.sum(), 1)
    fb = 100 * hb / max(hb.sum(), 1)
    r, p = stats.pearsonr(fa, fb) if ha.sum() and hb.sum() else (float("nan"), float("nan"))
    return {"bin_bp": bin_bp, "starts": [int(x) for x in edges[:-1]],
            "pct_a": [round(float(x), 3) for x in fa], "pct_b": [round(float(x), 3) for x in fb],
            "pearson_r": float(r), "pearson_p": float(p)}


def _hot(calls: list[tuple[int, str]], ref: str) -> tuple[float, float]:
    gc = [p for p, _ in calls if ref[p - 1].upper() in "GC"]
    at = [p for p, _ in calls if ref[p - 1].upper() in "AT"]
    hgc = sum(motifs.is_site("WRC", ref, p - 1) or motifs.is_site("GYW", ref, p - 1) for p in gc)
    hat = sum(motifs.is_site("WA", ref, p - 1) or motifs.is_site("TW", ref, p - 1) for p in at)
    return (100 * hgc / len(gc) if gc else float("nan"),
            100 * hat / len(at) if at else float("nan"))


def hotspots(calls_a: dict, calls_b: dict, ref: str, lo: int, hi: int) -> dict:
    """calls_x: animal -> [(pos1, alt)]."""
    gc = [i for i in range(lo - 1, hi) if ref[i].upper() in "GC"]
    at = [i for i in range(lo - 1, hi) if ref[i].upper() in "AT"]
    exp = (100 * sum(motifs.is_site("WRC", ref, i) or motifs.is_site("GYW", ref, i) for i in gc)
           / max(len(gc), 1),
           100 * sum(motifs.is_site("WA", ref, i) or motifs.is_site("TW", ref, i) for i in at)
           / max(len(at), 1))
    animals = [k for k in calls_a if k in calls_b]
    per = {k: {"a": _hot(calls_a[k], ref), "b": _hot(calls_b[k], ref)} for k in animals}

    def mean(side, j):
        v = [per[k][side][j] for k in animals if not math.isnan(per[k][side][j])]
        return float(np.mean(v)) if v else float("nan")
    return {"per_animal": {k: {"gc_wrc_gyw_a": v["a"][0], "gc_wrc_gyw_b": v["b"][0],
                               "at_wa_tw_a": v["a"][1], "at_wa_tw_b": v["b"][1]}
                           for k, v in per.items()},
            "mean": {"gc_wrc_gyw_a": mean("a", 0), "gc_wrc_gyw_b": mean("b", 0),
                     "at_wa_tw_a": mean("a", 1), "at_wa_tw_b": mean("b", 1)},
            "expected": {"gc_wrc_gyw": exp[0], "at_wa_tw": exp[1]}}


def from_results(res_a: dict, res_b: dict, pairs: dict[str, str] | None = None,
                 label_a: str = "A", label_b: str = "B", bin_bp: int = 25,
                 window: tuple[int, int] | None = None) -> dict:
    """Compare two project results (project.run_project or .public).

    `pairs` maps a sample of A to the same animal's sample in B (default: the
    samples with the same name in both). Controls can be paired too; the
    spectrum, positional and hotspot comparisons use only non-control animals,
    which carry the mutations.
    """
    sa = {s["sample"]: s for s in res_a["samples"]}
    sb = {s["sample"]: s for s in res_b["samples"]}
    pairs = pairs or {k: k for k in sa if k in sb}
    missing = [k for k in pairs if k not in sa] + [v for v in pairs.values() if v not in sb]
    if missing:
        raise ValueError(f"not in the results: {', '.join(missing)}")
    if not pairs:
        raise ValueError("no sample is in both results; give the pairs")
    wa = window or tuple(res_a["params"]["analysis_window"])
    wb = tuple(res_b["params"]["analysis_window"])
    lo, hi = max(wa[0], wb[0]), min(wa[1], wb[1])
    ref = res_a["reference"]
    if ref != res_b["reference"]:
        raise ValueError("the two results use different references; positions cannot be paired")
    mutated = [k for k in pairs if not sa[k].get("control")]

    def calls(res, name):
        return [(v["pos"], v["alt"]) for v in res["variants"]
                if v["sample"] == name and lo <= v["pos"] <= hi]

    def spec(res, name):
        return dict(Counter(f"{v['ref']}>{v['alt']}" for v in res["variants"]
                            if v["sample"] == name and lo <= v["pos"] <= hi))
    ca = {k: calls(res_a, k) for k in mutated}
    cb = {k: calls(res_b, pairs[k]) for k in mutated}
    out = {
        "labels": [label_a, label_b], "pairs": pairs, "window": [lo, hi],
        "spectrum": spectrum({k: spec(res_a, k) for k in mutated},
                             {k: spec(res_b, pairs[k]) for k in mutated}),
        "positional": positional([p for v in ca.values() for p, _ in v],
                                 [p for v in cb.values() for p, _ in v], lo, hi, bin_bp),
        "hotspots": hotspots(ca, cb, ref, lo, hi),
        "density": [{"animal": k, "control": bool(sa[k].get("control")),
                     "density_a": sa[k].get("unique_density_x1e5"),
                     "density_b": sb[pairs[k]].get("unique_density_x1e5")} for k in pairs],
        "caveats": CAVEATS,
    }
    return out


def report_markdown(c: dict) -> str:
    a, b = c["labels"]
    s = c["spectrum"]
    lines = [f"# {a} vs {b}", "",
             f"Window {c['window'][0]}-{c['window'][1]}; animals: {', '.join(s['animals'])}; "
             f"{s['calls_a']} {a} and {s['calls_b']} {b} unique calls.", "",
             f"Spectrum: Pearson r = {s['pearson_r']:.2f} [{s['pearson_ci'][0]:.2f}-"
             f"{s['pearson_ci'][1]:.2f}], slope = {s['slope']:.2f}; mean |difference| "
             f"{s['mean_abs_diff_pp']:.2f} pp, largest {s['max_abs_diff_pp']:.2f} pp."]
    if "paired" in s:
        pr = s["paired"]
        t = pr["ts_tv"]
        lines += [f"Paired: interaction F({pr['interaction_df'][0]},{pr['interaction_df'][1]}) = "
                  f"{pr['interaction_F']:.2f}, p = {pr['interaction_p']:.3g}; smallest Sidak-adjusted "
                  f"per-type p = {pr['min_p_sidak']:.3g}; G:C Ts / Tv adjusted p = "
                  f"{t['gc_transitions']['p_sidak']:.3g} / {t['gc_transversions']['p_sidak']:.3g}."]
    p = c["positional"]
    h = c["hotspots"]["mean"]
    lines += [f"Along the window ({p['bin_bp']}-bp bins): r = {p['pearson_r']:.2f}.",
              f"Hotspots: G:C at WRC/GYW {h['gc_wrc_gyw_a']:.0f}% vs {h['gc_wrc_gyw_b']:.0f}%; "
              f"A:T at WA/TW {h['at_wa_tw_a']:.0f}% vs {h['at_wa_tw_b']:.0f}%.", "",
              "| type | " + a + " % [95% CI] | " + b + " % [95% CI] |", "|---|---|---|"]
    for t in s["types"]:
        lines.append(f"| {t['type']} | {t['pct_a']:.1f} [{t['ci_a'][0]:.1f}-{t['ci_a'][1]:.1f}] | "
                     f"{t['pct_b']:.1f} [{t['ci_b'][0]:.1f}-{t['ci_b'][1]:.1f}] |")
    lines += ["", "| animal | " + a + " density | " + b + " density |", "|---|---|---|"]
    for d in c["density"]:
        lines.append(f"| {d['animal']}{' (control)' if d['control'] else ''} | {d['density_a']} | "
                     f"{d['density_b']} |")
    lines += [""] + [f"- {x}" for x in c["caveats"]]
    return "\n".join(lines) + "\n"


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="shmtrace-concordance", description=__doc__.split("\n\n")[0])
    ap.add_argument("a", help="project.json of the first measurement (e.g. Sanger)")
    ap.add_argument("b", help="project.json of the second (e.g. ONT)")
    ap.add_argument("--pair", action="append", default=[],
                    help="A_SAMPLE=B_SAMPLE, once per animal (default: matching names)")
    ap.add_argument("--labels", default="A,B")
    ap.add_argument("-o", "--out", help="write the result as JSON here")
    a = ap.parse_args(argv)
    ra, rb = json.load(open(a.a)), json.load(open(a.b))
    pairs = dict(x.split("=", 1) for x in a.pair) or None
    la, lb = a.labels.split(",", 1)
    c = from_results(ra, rb, pairs, la, lb)
    if a.out:
        json.dump(c, open(a.out, "w"), indent=1)
    print(report_markdown(c))


if __name__ == "__main__":
    main()
