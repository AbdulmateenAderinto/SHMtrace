"""Data for the SciCanvas Figure S8 (figS8_cohort2.js).

    python docs/calibration/crm/scicanvas/cohort2_data.py

Reads what figs/figS8_sanger_cohort2.py reads (out/sanger_cohort2.json and the transcribed published
tables), reduces each panel to the numbers it plots, and writes them into figS8_cohort2.js between the
DATA markers, so the builder stays one file that can also be pasted into SciCanvas.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "figs"))
from figS8_sanger_cohort2 import (  # noqa: E402
    MICE,
    PUB,
    RES,
    TYPES,
    _bg_share,
    gc,
    one_tailed,
    pub_freq,
    pub_spec,
)

BUILDER = HERE / "figS8_cohort2.js"
GROUPS = ("WT", "Maea_KO", "Aicda_KO")


def main():
    per = RES["per_mouse"]
    mine = {g: [m for m in per if per[m]["group"] == g] for g in GROUPS}
    pub_g = {g: [m for m in MICE if PUB["mice"][m]["group"] == g] for g in GROUPS}

    # A: frequency per mouse, SHMtrace against published
    A = {g: {"sh": [round(per[m]["freq_x1e4"], 2) for m in mine[g]],
             "pub": [round(pub_freq(m), 2) for m in pub_g[g]]} for g in GROUPS}

    # B: per mouse, SHMtrace against published; mice are matched by the sheet's names, which are the paper's
    xs = [pub_freq(m) for m in MICE]
    ys = [per[m]["freq_x1e4"] for m in MICE]
    B = {"points": [[round(x, 3), round(y, 3), PUB["mice"][m]["group"]] for m, x, y in zip(MICE, xs, ys, strict=True)],
         "rho": round(float(stats.spearmanr(xs, ys)[0]), 2),
         "ratio": round(float(10 ** np.mean(np.log10(np.array(ys) / np.array(xs)))), 2)}

    # C: mutations per clone, binned as in the published Fig. 2g
    C = {"bins": PUB["fig2g_bins"], "rows": []}
    for g in ("WT", "Maea_KO"):
        cnt = np.sum([per[m]["load_bins"] for m in mine[g]], axis=0)
        C["rows"].append({"group": g, "who": "published", "n": PUB["fig2g_n"][g],
                          "pct": [round(float(v), 2) for v in PUB["fig2g_percent"][g]]})
        C["rows"].append({"group": g, "who": "SHMtrace", "n": int(cnt.sum()),
                          "pct": [round(float(v), 2) for v in 100 * cnt / cnt.sum()]})

    # D: G:C transitions and transversions per mouse
    D = []
    for g in ("WT", "Maea_KO"):
        for kind, i in (("Ts", 0), ("Tv", 1)):
            D.append({"group": g, "kind": kind, "sh": [round(gc(per[m]["spectrum"])[i], 2) for m in mine[g]],
                      "pub": [round(gc(pub_spec(m))[i], 2) for m in pub_g[g]]})

    # E: pooled 12-type spectrum
    E = {}
    allx, ally = [], []
    for g in ("WT", "Maea_KO"):
        sh = {t: sum(per[m]["spectrum"].get(t, 0) for m in mine[g]) for t in TYPES}
        pb = {t: sum(pub_spec(m)[t] for m in pub_g[g]) for t in TYPES}
        x = [100 * pb[t] / sum(pb.values()) for t in TYPES]
        y = [100 * sh[t] / sum(sh.values()) for t in TYPES]
        allx += x
        ally += y
        E[g] = [[round(a, 3), round(b, 3), t] for a, b, t in zip(x, y, TYPES, strict=True)]
    r = float(stats.pearsonr(allx, ally)[0])

    # F: the published comparisons recovered
    fw, fk = [per[m]["freq_x1e4"] for m in mine["WT"]], [per[m]["freq_x1e4"] for m in mine["Maea_KO"]]
    pw, pk = [pub_freq(m) for m in pub_g["WT"]], [pub_freq(m) for m in pub_g["Maea_KO"]]
    g_of = lambda i, g: [gc(per[m]["spectrum"])[i] for m in mine[g]]                 # noqa: E731
    lineages = [x for x in RES["project"]["qc"]["duplicate_groups"] if x["n_samples"] > 1]
    F = [["Wt / Maea frequency", f"{np.mean(pw) / np.mean(pk):.1f}×", f"{np.mean(fw) / np.mean(fk):.1f}×", 0],
         ["P (one-tailed t)", f"{PUB['fig2h_p']['WT_vs_Maea_KO']:.3f}", f"{one_tailed(fw, fk):.4f}", 1],
         ["G:C Ts P", f"{PUB['fig2i_p']['Ts']:.4f}", f"{stats.ttest_ind(g_of(0, 'WT'), g_of(0, 'Maea_KO')).pvalue:.4f}", 0],
         ["G:C Tv P", f"{PUB['fig2i_p']['Tv']:.4f}", f"{stats.ttest_ind(g_of(1, 'WT'), g_of(1, 'Maea_KO')).pvalue:.3f}", 0],
         ["Mice, Wt / Maea / null", "4 / 4 / 3", f"{len(mine['WT'])} / {len(mine['Maea_KO'])} / {len(mine['Aicda_KO'])}", 0],
         ["Clones, Wt and Maea", str(sum(PUB["mice"][m]["sequences"] for m in MICE)),
          str(sum(per[m]["clones"] for m in mine["WT"] + mine["Maea_KO"])), 0],
         ["Background explains, Wt", "–", _bg_share("WT"), 0],
         ["...                 Maea", "–", _bg_share("Maea_KO"), 1],
         ["Lineages under 2+ mice", "–", f"{len(lineages)} ({sum(x['n_clones'] for x in lineages)} clones)", 0],
         ["Analysis", "hand-edited", f"{RES['seconds']:.0f} s", 0]]
    data = {"A": A, "B": B, "C": C, "D": D, "E": E, "r": round(r, 2), "F": F}
    src = BUILDER.read_text()
    new = re.sub(r"/\*DATA\*/.*?/\*END\*/", lambda _m: "/*DATA*/" + json.dumps(data, separators=(",", ":")) + "/*END*/",
                 src, flags=re.S)
    BUILDER.write_text(new)
    print(f"{BUILDER.name}: data {len(new) - len(src):+d} bytes; rho {B['rho']}, ratio {B['ratio']}, r {r:.2f}")


if __name__ == "__main__":
    main()
