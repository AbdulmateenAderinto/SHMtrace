"""Figure S14, the matplotlib draft: a second published Sanger study re-analysed automatically
(Barbulescu et al., 2025, Mkln1).

The published Figure S14 is built in SciCanvas (scicanvas/figS15_mkln.js), from the numbers that
scicanvas/mkln_data.py takes from this module; running this file draws the draft into out/.

A    mutation frequency per mouse: SHMtrace from the chromatograms against the published values
B    G:C transitions and transversions per mouse, as in the published Supplementary Figure S2C
C    pooled 12-type substitution spectrum, SHMtrace against published, per genotype
D    the published comparisons recovered, and what the read filters removed
Published values are replotted from the paper's per-mouse tables (Supplementary Figure S2A);
SHMtrace's come from out/sanger_mkln.json (sanger_cohort2.py mkln).
"""
import json
from pathlib import Path

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy import stats
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, setup

PUB = json.load(open(Path(__file__).resolve().parents[1] / "published/barbulescu2025_jh4_mkln.json"))
RES = json.load(open(DATA / "sanger_mkln.json"))
SH_C, PUB_C = "#1f6f8b", "#a3621b"
KO, KOF = "#7a4f9a", "#bda3d4"
GCOL = {"WT": (CASE, CASEF), "Mkln1_KO": (KO, KOF), "Aicda_KO": (CTRL, CTRLF)}
GLAB = {"WT": "Wt", "Mkln1_KO": "Mkln1$^{-/-}$", "Aicda_KO": "Aicda$^{-/-}$"}
B = "ATCG"                                    # the published tables are ordered A, T, C, G
TYPES = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
TS, TV = ("C>T", "G>A"), ("C>A", "C>G", "G>C", "G>T")


def pub_spec(m):
    """{'A>T': n, ...} for a published mouse; its matrix is rows = To, columns = From, as printed."""
    mat = PUB["mice"][m]["matrix"]
    return {f"{B[j]}>{B[i]}": mat[i][j] for i in range(4) for j in range(4) if i != j}


def gc(sp):
    tot = sum(sp.values())
    return (100 * sum(sp.get(k, 0) for k in TS) / tot if tot else None,
            100 * sum(sp.get(k, 0) for k in TV) / tot if tot else None)


def _bg_share(group):
    """The share of a group's calls that the AID-null clones' own error rate accounts for."""
    v = [RES["per_mouse"][m]["expected_false_pct"] for m in RES["per_mouse"]
         if RES["per_mouse"][m]["group"] == group and RES["per_mouse"][m]["expected_false_pct"] is not None]
    return f"{min(v):.0f}-{max(v):.0f}%" if v else "-"


def main():
    setup()
    per = RES["per_mouse"]
    pub_mice = list(PUB["mice"])
    pub_g = {m: PUB["mice"][m]["group"] for m in pub_mice}
    fig = plt.figure(figsize=(FULL, 62 * MM))
    gs = gridspec.GridSpec(1, 4, figure=fig, left=0.065, right=0.99, top=0.84, bottom=0.2, wspace=0.62,
                           width_ratios=[1.05, 1.0, 0.95, 1.55])

    # A: frequency per mouse
    ax = fig.add_subplot(gs[0, 0])
    for i, g in enumerate(("WT", "Mkln1_KO", "Aicda_KO")):
        e, f = GCOL[g]
        sh = [per[m]["freq_x1e4"] for m in per if per[m]["group"] == g]
        ax.scatter(i - 0.17 + np.linspace(-0.05, 0.05, len(sh)), sh, s=13, facecolor=f, edgecolor=e, lw=0.6, zorder=3)
        ax.plot([i - 0.3, i - 0.04], [np.mean(sh)] * 2, color=SH_C, lw=1.3)
        pb = [1e4 * PUB["mice"][m]["f_mut"] for m in pub_mice if pub_g[m] == g]
        if pb:
            ax.scatter(i + 0.17 + np.linspace(-0.05, 0.05, len(pb)), pb, s=13, facecolor="white", edgecolor=e,
                       lw=0.6, marker="s", zorder=3)
            ax.plot([i + 0.04, i + 0.3], [np.mean(pb)] * 2, color=PUB_C, lw=1.3)
    ax.set_xticks(range(3), [GLAB[g] for g in ("WT", "Mkln1_KO", "Aicda_KO")])
    ax.set_ylabel("Mutations per 10$^4$ bp")
    ax.set_ylim(0, 215)
    ax.set_title("Frequency per mouse", fontsize=7)
    ax.scatter([], [], s=13, facecolor=SOFT, edgecolor=INK, lw=0.6, label="SHMtrace")
    ax.scatter([], [], s=13, facecolor="white", edgecolor=INK, lw=0.6, marker="s", label="Published")
    ax.legend(loc="upper right", fontsize=6, handletextpad=0.2, borderaxespad=0.1)
    letter(ax, "A", dx=-0.24)

    # B: G:C transitions and transversions
    ax = fig.add_subplot(gs[0, 1])
    ticks = []
    for pos, (g, kind, idx) in enumerate([(g, k, i) for g in ("WT", "Mkln1_KO") for k, i in (("Ts", 0), ("Tv", 1))]):
        e, f = GCOL[g]
        sh = [gc(per[m]["spectrum"])[idx] for m in per if per[m]["group"] == g]
        pb = [gc(pub_spec(m))[idx] for m in pub_mice if pub_g[m] == g]
        for dx, vals, mk, fc, col in ((-0.17, sh, "o", f, SH_C), (0.17, pb, "s", "white", PUB_C)):
            ax.scatter(pos + dx + np.linspace(-0.04, 0.04, len(vals)), vals, s=11, marker=mk, facecolor=fc,
                       edgecolor=e, lw=0.6, zorder=3)
            ax.plot([pos + dx - 0.12, pos + dx + 0.12], [np.mean(vals)] * 2, color=col, lw=1.3)
        ticks.append(f"{kind}\n{GLAB[g].split('$')[0]}")
    ax.set_xticks(range(4), ticks)
    ax.set_ylim(0, 40)
    ax.set_ylabel("% of all mutations")
    ax.set_title("G:C Ts and Tv", fontsize=7)
    letter(ax, "B", dx=-0.3)

    # C: pooled spectrum
    ax = fig.add_subplot(gs[0, 2])
    allx, ally = [], []
    for g in ("WT", "Mkln1_KO"):
        e, f = GCOL[g]
        sh = {t: sum(per[m]["spectrum"].get(t, 0) for m in per if per[m]["group"] == g) for t in TYPES}
        pb = {t: sum(pub_spec(m)[t] for m in pub_mice if pub_g[m] == g) for t in TYPES}
        x = [100 * pb[t] / sum(pb.values()) for t in TYPES]
        y = [100 * sh[t] / sum(sh.values()) for t in TYPES]
        allx += x
        ally += y
        ax.scatter(x, y, s=13, facecolor=f, edgecolor=e, lw=0.6, zorder=3, label=GLAB[g])
    ax.plot([0, 28], [0, 28], color=SOFT, lw=0.6, ls="--")
    ax.set_xlim(0, 28)
    ax.set_ylim(0, 28)
    r = stats.pearsonr(allx, ally)[0]
    ax.text(0.05, 0.95, f"r = {r:.2f}", transform=ax.transAxes, va="top", fontsize=6, color=MID)
    ax.set_xlabel("Published (%)")
    ax.set_ylabel("SHMtrace (%)")
    ax.set_title("Substitution spectrum", fontsize=7)
    ax.legend(loc="lower right", fontsize=6, handletextpad=0.2)
    letter(ax, "C", dx=-0.26)

    # D: what was recovered
    ax = fig.add_subplot(gs[0, 3])
    ax.axis("off")
    wt = [m for m in per if per[m]["group"] == "WT"]
    mk = [m for m in per if per[m]["group"] == "Mkln1_KO"]
    shf = lambda ms: [per[m]["freq_x1e4"] for m in ms]                                      # noqa: E731
    pbf = lambda g: [1e4 * PUB["mice"][m]["f_mut"] for m in pub_mice if pub_g[m] == g]      # noqa: E731
    p_sh = stats.ttest_ind(shf(wt), shf(mk), alternative="greater").pvalue
    ts_sh = stats.ttest_ind([gc(per[m]["spectrum"])[0] for m in wt], [gc(per[m]["spectrum"])[0] for m in mk]).pvalue
    tv_sh = stats.ttest_ind([gc(per[m]["spectrum"])[1] for m in wt], [gc(per[m]["spectrum"])[1] for m in mk]).pvalue
    pw_ts = [gc(pub_spec(m))[0] for m in pub_mice if pub_g[m] == "WT"]
    pk_ts = [gc(pub_spec(m))[0] for m in pub_mice if pub_g[m] == "Mkln1_KO"]
    pw_tv = [gc(pub_spec(m))[1] for m in pub_mice if pub_g[m] == "WT"]
    pk_tv = [gc(pub_spec(m))[1] for m in pub_mice if pub_g[m] == "Mkln1_KO"]
    p_pub = stats.ttest_ind(pbf("WT"), pbf("Mkln1_KO"), alternative="greater").pvalue
    lines = [("", "Published", "SHMtrace"),
             ("Wt / Mkln1$^{-/-}$", f"{np.mean(pbf('WT')) / np.mean(pbf('Mkln1_KO')):.1f}×",
              f"{np.mean(shf(wt)) / np.mean(shf(mk)):.1f}×"),
             ("  P (one-tailed t)", f"{p_pub:.4f}", f"{p_sh:.4f}"),
             ("G:C Ts P", f"{stats.ttest_ind(pw_ts, pk_ts).pvalue:.2f}", f"{ts_sh:.2f}"),
             ("G:C Tv P", f"{stats.ttest_ind(pw_tv, pk_tv).pvalue:.3f}", f"{tv_sh:.3f}"),
             ("  ... S2C bars", f"{PUB['figS2c_p']['Wt_vs_Mkln1']['Ts']} / "
              f"{PUB['figS2c_p']['Wt_vs_Mkln1']['Tv']}$^{{a}}$", "–"),
             ("Clones, Wt / KO", f"{sum(PUB['mice'][m]['sequences'] for m in pub_mice if pub_g[m] == 'WT')}"
              f" / {sum(PUB['mice'][m]['sequences'] for m in pub_mice if pub_g[m] == 'Mkln1_KO')}",
              f"{sum(per[m]['clones'] for m in wt)} / {sum(per[m]['clones'] for m in mk)}"),
             ("Paralog clones removed", "–", f"{RES['excluded'].get('contaminant', 0)}"),
             ("Background explains, Wt", "–", _bg_share("WT")),
             ("  ...   Mkln1$^{-/-}$", "–", _bg_share("Mkln1_KO")),
             ("Analysis", "Benchling", f"{RES['seconds']:.0f} s")]
    for k, (a, b, c) in enumerate(lines):
        y = 0.96 - k * 0.088
        w = "bold" if k == 0 else "normal"
        kw = dict(fontsize=6, fontweight=w, va="top", transform=ax.transAxes, clip_on=False)
        ax.text(0.0, y, a, color=INK, **kw)
        ax.text(0.62, y, b, color=PUB_C if k else INK, ha="center", **kw)
        ax.text(0.88, y, c, color=SH_C if k else INK, ha="center", **kw)
    ax.plot([0, 0.98], [0.895, 0.895], color=SOFT, lw=0.5, transform=ax.transAxes)
    ax.set_title("Published tests", fontsize=7)
    ax.text(0.0, -0.06, "$^{a}$ Ts / Tv; fit the tables only transposed", transform=ax.transAxes,
            fontsize=6, color=MID, va="top")
    letter(ax, "D", dx=-0.04)
    fig.savefig(DATA / "FigureS15_sanger_mkln_matplotlib.png", dpi=300)
    return {"spec_r": float(r), "fold_sh": float(np.mean(shf(wt)) / np.mean(shf(mk))),
            "fold_pub": float(np.mean(pbf("WT")) / np.mean(pbf("Mkln1_KO"))), "p_sh": float(p_sh),
            "ts_p": float(ts_sh), "tv_p": float(tv_sh)}


if __name__ == "__main__":
    print(main())
