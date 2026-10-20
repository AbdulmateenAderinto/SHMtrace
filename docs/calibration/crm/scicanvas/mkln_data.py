"""Data for the SciCanvas Figure S9 (figS14_mkln.js).

    python docs/calibration/crm/scicanvas/mkln_data.py

Reads what figs/figS14_sanger_mkln.py reads (out/sanger_mkln.json and the transcribed published tables),
reduces each panel to the numbers it plots, and writes them into figS14_mkln.js between the DATA markers,
so the builder stays one file that can also be pasted into SciCanvas.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "figs"))
sys.path.insert(0, str(HERE))
from figS14_sanger_mkln import PUB, RES, TYPES, _bg_share, gc, pub_spec  # noqa: E402
from inject import inject  # noqa: E402

BUILDER = HERE / "figS14_mkln.js"
GROUPS = ("WT", "Mkln1_KO", "Aicda_KO")


def main():
    per = RES["per_mouse"]
    pub_mice = list(PUB["mice"])
    pub_g = {m: PUB["mice"][m]["group"] for m in pub_mice}
    mine = {g: [m for m in per if per[m]["group"] == g] for g in GROUPS}
    shf = {g: [per[m]["freq_x1e4"] for m in mine[g]] for g in GROUPS}
    pbf = {g: [1e4 * PUB["mice"][m]["f_mut"] for m in pub_mice if pub_g[m] == g] for g in GROUPS}

    # A: frequency per mouse, SHMtrace against published
    panelA = {g: {"sh": [round(v, 2) for v in shf[g]], "pub": [round(v, 2) for v in pbf[g]]} for g in GROUPS}

    # B: G:C transitions and transversions per mouse
    panelB = []
    for g in ("WT", "Mkln1_KO"):
        for kind, i in (("Ts", 0), ("Tv", 1)):
            panelB.append({"group": g, "kind": kind,
                           "sh": [round(gc(per[m]["spectrum"])[i], 2) for m in mine[g]],
                           "pub": [round(gc(pub_spec(m))[i], 2) for m in pub_mice if pub_g[m] == g]})

    # C: pooled 12-type spectrum, as percentages
    panelC, rg = {}, {}
    allx, ally = [], []
    for g in ("WT", "Mkln1_KO"):
        sh = {t: sum(per[m]["spectrum"].get(t, 0) for m in mine[g]) for t in TYPES}
        pb = {t: sum(pub_spec(m)[t] for m in pub_mice if pub_g[m] == g) for t in TYPES}
        x = [100 * pb[t] / sum(pb.values()) for t in TYPES]
        y = [100 * sh[t] / sum(sh.values()) for t in TYPES]
        allx += x
        ally += y
        panelC[g] = [[round(a, 3), round(b, 3), t] for a, b, t in zip(x, y, TYPES, strict=True)]
        rg[g] = round(float(stats.pearsonr(x, y)[0]), 2)
    r = float(stats.pearsonr(allx, ally)[0])

    # D: the published comparisons, published against SHMtrace
    t2 = lambda a, b, **k: float(stats.ttest_ind(a, b, **k).pvalue)                 # noqa: E731
    g_of = lambda i, g: [gc(per[m]["spectrum"])[i] for m in mine[g]]                # noqa: E731
    p_of = lambda i, g: [gc(pub_spec(m))[i] for m in pub_mice if pub_g[m] == g]     # noqa: E731
    seqs = {g: sum(PUB["mice"][m]["sequences"] for m in pub_mice if pub_g[m] == g) for g in ("WT", "Mkln1_KO")}
    s2c = PUB["figS2c_p"]["Wt_vs_Mkln1"]
    panelD = [
        ["Wt / Mkln1KO", f"{np.mean(pbf['WT']) / np.mean(pbf['Mkln1_KO']):.1f}×",
         f"{np.mean(shf['WT']) / np.mean(shf['Mkln1_KO']):.1f}×", 0],
        ["P (one-tailed t)ᵇ", f"{t2(pbf['WT'], pbf['Mkln1_KO'], alternative='greater'):.4f}",
         f"{t2(shf['WT'], shf['Mkln1_KO'], alternative='greater'):.4f}", 1],
        ["G:C Ts P", f"{t2(p_of(0, 'WT'), p_of(0, 'Mkln1_KO')):.2f}", f"{t2(g_of(0, 'WT'), g_of(0, 'Mkln1_KO')):.2f}", 0],
        ["G:C Tv P", f"{t2(p_of(1, 'WT'), p_of(1, 'Mkln1_KO')):.3f}", f"{t2(g_of(1, 'WT'), g_of(1, 'Mkln1_KO')):.3f}", 0],
        ["... S2C barsᵃ", f"{s2c['Ts']:.3f} / {s2c['Tv']:.2f}", "–", 1],
        ["Clones, Wt / KO", f"{seqs['WT']} / {seqs['Mkln1_KO']}",
         f"{sum(per[m]['clones'] for m in mine['WT'])} / {sum(per[m]['clones'] for m in mine['Mkln1_KO'])}", 0],
        ["Paralog clones out", "–", f"{RES['excluded'].get('contaminant', 0)}", 0],
        ["Background explains, Wt", "–", _bg_share("WT"), 0],
        ["... Mkln1KO", "–", _bg_share("Mkln1_KO"), 1],
        ["Analysis", "Benchling", f"{RES['seconds']:.0f} s", 0],
    ]
    data = {"A": panelA, "B": panelB, "C": panelC, "r": round(r, 2), "r_group": rg, "D": panelD,
            "traces": RES["traces"], "clones": RES["clones_analysed"]}
    n = inject(BUILDER, data)
    print(f"{BUILDER.name}: data {n:+d} bytes; r = {r:.2f}; {len(panelD)} rows in D")


if __name__ == "__main__":
    main()
