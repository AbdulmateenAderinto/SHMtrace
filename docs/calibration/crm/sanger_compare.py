"""Table S4: SHMtrace against the two published Sanger analyses, metric by metric.

Every cell is "wild type -> knockout", with an arrow for the direction and stars for a two-tailed
test across mice (Welch on log10 for frequency, Student t on per-mouse percentages otherwise):
*** p < 0.001, ** p < 0.01, * p < 0.05, ns otherwise. Published values are computed from the papers'
per-mouse tables (published/), SHMtrace's from out/sanger_<study>.json. N.R., not reported.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
TS, TV = ("C>T", "G>A"), ("C>A", "C>G", "G>C", "G>T")
AT = ("A>T", "A>C", "A>G", "T>A", "T>C", "T>G")
STUDIES = {
    "cohort2": {"label": "Maea^{fl/fl}Mb1^{cre}", "ko": "Maea_KO", "pub": "published/barbulescu2024_jh4.json",
                "bases": "GCAT", "rows_are": "From"},
    "mkln": {"label": "Mkln1^{−/−}", "ko": "Mkln1_KO", "pub": "published/barbulescu2025_jh4_mkln.json",
             "bases": "ATCG", "rows_are": "To"},
}


def _stars(p):
    return "ns" if p is None or p >= 0.05 else ("*" if p >= 0.01 else ("**" if p >= 0.001 else "***"))


def _cell(wt, ko, p, fmt):
    a, b = np.mean(wt), np.mean(ko)
    arrow = "" if p is None or p >= 0.05 else ("↓" if b < a else "↑")
    return f"{fmt.format(a)} → {fmt.format(b)} {arrow}{_stars(p)}".replace("  ", " ")


def _share(sp, keys):
    tot = sum(sp.values())
    return 100 * sum(sp.get(k, 0) for k in keys) / tot if tot else float("nan")


def _pub_spec(pub, m, bases, rows_are):
    mat = pub["mice"][m]["matrix"]
    return {(f"{bases[i]}>{bases[j]}" if rows_are == "From" else f"{bases[j]}>{bases[i]}"): mat[i][j]
            for i in range(4) for j in range(4) if i != j}


def _pub_freq(pub, m):
    x = pub["mice"][m]
    return 1e4 * x["f_mut"] if "f_mut" in x else 1e4 * x["unique_mutations"] / x["nt"]


def _ttest(a, b, log=False):
    if len(a) < 2 or len(b) < 2:
        return None
    if log:
        a, b = [math.log10(x) for x in a], [math.log10(x) for x in b]
        return float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    return float(stats.ttest_ind(a, b).pvalue)


def study_columns(key):
    """{metric: (shmtrace cell, published cell)} for one study."""
    s = STUDIES[key]
    res = json.load(open(HERE / f"out/sanger_{key}.json"))
    pub = json.load(open(HERE / s["pub"]))
    per = res["per_mouse"]
    samp = {x["sample"]: x for x in res["project"]["samples"]}
    wt = [m for m in per if per[m]["group"] == "WT"]
    ko = [m for m in per if per[m]["group"] == s["ko"]]
    pw = [m for m in pub["mice"] if pub["mice"][m]["group"] == "WT"]
    pk = [m for m in pub["mice"] if pub["mice"][m]["group"] == s["ko"]]
    sp = {m: per[m]["spectrum"] for m in per}
    ps = {m: _pub_spec(pub, m, s["bases"], s["rows_are"]) for m in pub["mice"]}
    out = {}
    fw, fk = [per[m]["freq_x1e4"] for m in wt], [per[m]["freq_x1e4"] for m in ko]
    gw, gk = [_pub_freq(pub, m) for m in pw], [_pub_freq(pub, m) for m in pk]
    out["Mutation frequency (×10^{−4})"] = (_cell(fw, fk, _ttest(fw, fk, log=True), "{:.0f}"),
                                             _cell(gw, gk, _ttest(gw, gk, log=True), "{:.0f}"))
    for name, keys in (("A:T mutations (% of all)", AT), ("G:C transitions (% of all)", TS),
                       ("G:C transversions (% of all)", TV)):
        a, b = [_share(sp[m], keys) for m in wt], [_share(sp[m], keys) for m in ko]
        c, d = [_share(ps[m], keys) for m in pw], [_share(ps[m], keys) for m in pk]
        out[name] = (_cell(a, b, _ttest(a, b), "{:.0f}%"), _cell(c, d, _ttest(c, d), "{:.0f}%"))
    for name, k in (("G:C mutations at WRC/GYW (%)", "pct_GC_in_WRC_GYW"),
                    ("A:T mutations at WA/TW (%)", "pct_AT_in_WA_TW")):
        a = [samp[m][k] for m in wt if samp[m].get(k) is not None]
        b = [samp[m][k] for m in ko if samp[m].get(k) is not None]
        out[name] = (_cell(a, b, _ttest(a, b), "{:.0f}%"), "N.R.")
    a = [per[m]["expected_false_pct"] for m in wt if per[m]["expected_false_pct"] is not None]
    b = [per[m]["expected_false_pct"] for m in ko if per[m]["expected_false_pct"] is not None]
    out["Calls the background explains (%)"] = (f"{np.mean(a):.0f}% → {np.mean(b):.0f}%", "N.R.")
    out["Mice (Wt / knockout)"] = (f"{len(wt)} / {len(ko)}", f"{len(pw)} / {len(pk)}")
    out["Clones (Wt / knockout)"] = (f"{sum(per[m]['clones'] for m in wt)} / {sum(per[m]['clones'] for m in ko)}",
                                     f"{sum(pub['mice'][m]['sequences'] for m in pw)} / "
                                     f"{sum(pub['mice'][m]['sequences'] for m in pk)}")
    return out


def table_rows():
    a, b = study_columns("cohort2"), study_columns("mkln")
    return [[metric, *a[metric], *b[metric]] for metric in a]


HEADER = ["Metric", "Maea: SHMtrace", "Maea: published", "Mkln1: SHMtrace", "Mkln1: published"]

if __name__ == "__main__":
    for r in table_rows():
        print(" | ".join(f"{c:28}" for c in r))
