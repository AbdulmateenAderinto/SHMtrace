"""Figure 5: nanopore and Sanger agree in the same mice.

A    unique-variant density per mouse, Sanger clones and nanopore reads (529-bp window)
B    pooled 12-type substitution spectrum, Sanger against nanopore (three Wt mice)
C    where calls fall along the window (25-bp bins)
D    hotspot use: G:C calls at WRC/GYW and A:T calls at WA/TW, against the window's share of such sites
E    the same comparison mouse by mouse: four substitution classes, Sanger against nanopore
     (the bookkeeping checks that were panel E are Figure S7, page 2: figS12_bookkeeping.py)
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy import stats
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]
SV = CAL / "figures/sanger_vs_nanopore/recommended_manuscript/stats.json"
SANGER = CAL / "manuscript_corrected/sanger_calls/corrected_sanger.json"
ONT = CAL / "figures/jh4_shmtrace/project.json"
TYPES = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
SANG_C, ONT_C = "#a3621b", "#1f6f8b"


def main():
    setup()
    st = json.load(open(SV))
    fig = plt.figure(figsize=(FULL, 112 * MM))
    gs = gridspec.GridSpec(2, 4, figure=fig, left=0.07, right=0.98, top=0.93, bottom=0.12, wspace=0.6, hspace=0.8,
                           width_ratios=[1.25, 1, 1, 1])
    # A: densities
    ax = fig.add_subplot(gs[0, 0])
    for j, (plat, d, _col) in enumerate((("Sanger", st["density_sanger"], SANG_C), ("Nanopore", st["density_ont"], ONT_C))):
        for k, (grp, f, e) in enumerate((("wt", CASEF, CASE), ("ko", CTRLF, CTRL))):
            vals = list(d[grp].values())
            x = j * 2.4 + k
            ax.scatter(x + np.linspace(-0.12, 0.12, len(vals)), vals, s=14, facecolor=f, edgecolor=e, lw=0.6, zorder=3)
            g = 10 ** np.mean(np.log10(vals))
            ax.plot([x - 0.25, x + 0.25], [g, g], color=e, lw=1.4)
        ax.text(j * 2.4 + 0.5, 3000, f"{d['fold']:.0f}×\np = {d['p_welch_log']:.4f}" if plat == "Nanopore"
                else f"{d['fold']:.0f}×\np = {d['p_welch_log']:.4f}", ha="center", fontsize=6, color=MID)
    ax.set_yscale("log")
    ax.set_ylim(0.2, 8000)
    ax.set_xticks([0, 1, 2.4, 3.4], ["Wt", "Aicda$^{-/-}$", "Wt", "Aicda$^{-/-}$"])
    ax.text(0.5, -0.27, "Sanger clones", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, color=SANG_C)
    ax.text(2.9, -0.27, "Nanopore reads", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, color=ONT_C)
    ax.set_ylabel("Unique variants per 10$^5$ bases")
    # the two densities are not on a common scale (STAR Methods), so the panel is read within a platform
    ax.set_title("Within-platform contrasts", fontsize=7)
    letter(ax, "A", dx=-0.3)
    # B: spectrum scatter
    ax = fig.add_subplot(gs[0, 1])
    sp = st["spectrum"]["pooled_props"]
    xs, ys = [sp["sanger"][t] for t in TYPES], [sp["ont"][t] for t in TYPES]
    gc = [t[0] in "GC" for t in TYPES]
    ax.scatter(xs, ys, s=14, c=["#2c5d8a" if g else "#d9952b" for g in gc], edgecolor="white", lw=0.5, zorder=3)
    for t, x, y in zip(TYPES, xs, ys, strict=False):
        if abs(x - y) > 1.5 or x > 12:
            # at 6 pt two labels would land on neighbouring points: those go above-left and below
            dx, dy, ha = {"A>C": (-3, 3, "right"), "C>T": (3, -8, "left")}.get(t, (3, -2, "left"))
            ax.annotate(t.replace(">", "→"), (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=6,
                        color=MID, ha=ha)
    lim = [0, max(xs + ys) * 1.12]
    ax.plot(lim, lim, color=SOFT, ls="--", lw=0.6)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    r = st["spectrum"]["pearson"]
    # a Fisher interval over 12 compositional categories is not a population confidence interval, so the
    # panel gives the correlation itself, as Results now does
    ax.text(0.04, 0.9, f"r = {r['r']:.2f}\nslope {r['slope']:.2f}",
            transform=ax.transAxes, fontsize=6)
    ax.set_xlabel("Sanger, % of calls")
    ax.set_ylabel("Nanopore, % of calls")
    ax.set_title("Substitution spectrum", fontsize=7)
    ax.scatter([], [], s=12, color="#2c5d8a", label="from G:C")
    ax.scatter([], [], s=12, color="#d9952b", label="from A:T")
    ax.legend(loc="lower right", fontsize=6, handletextpad=0.1)
    letter(ax, "B", dx=-0.32)
    # C: positional profile
    ax = fig.add_subplot(gs[0, 2:])
    edges = np.arange(101, 630 + 25, 25)
    sang = json.load(open(SANGER))
    s_pos = [v[0] for m in ("Wt 2", "Wt 3", "Wt 4") for v in sang[m]["variants"]]
    o_pos = [v["pos"] for v in json.load(open(ONT))["variants"] if v["sample"].startswith("WT")]
    hs, _ = np.histogram(s_pos, edges)
    ho, _ = np.histogram(o_pos, edges)
    mids = (edges[:-1] + edges[1:]) / 2
    ax.plot(mids, 100 * hs / hs.sum(), color=SANG_C, lw=1.1, marker="o", ms=2.2, label=f"Sanger (n = {hs.sum()})")
    ax.plot(mids, 100 * ho / ho.sum(), color=ONT_C, lw=1.1, marker="o", ms=2.2, label=f"Nanopore (n = {ho.sum()})")
    rr = stats.pearsonr(hs, ho)[0]
    ax.text(0.98, 0.78, f"r = {rr:.2f}", transform=ax.transAxes, ha="right", fontsize=6)
    ax.set_xlabel("J$_H$4 position (bp)")
    ax.set_ylabel("Share of calls (%)")
    ax.set_xlim(100, 630)
    ax.set_title("Where calls fall, 25-bp bins", fontsize=7)
    ax.legend(loc="upper center", fontsize=6, ncol=2)
    letter(ax, "C", dx=-0.12)
    # D: hotspots
    ax = fig.add_subplot(gs[1, 0])
    hot = st["beyond"]["hotspot"]
    exp = st["beyond"]["hotspot_expected"]
    for j, (_label, idx) in enumerate((("G:C at WRC/GYW", 0), ("A:T at WA/TW", 1))):
        for k, (plat, col) in enumerate((("sanger", SANG_C), ("ont", ONT_C))):
            vals = [m[idx] for m in hot[plat]]
            x = j * 2.6 + k
            ax.bar(x, np.mean(vals), 0.8, color=col, alpha=0.85)
            ax.scatter([x] * 3, vals, s=8, color=INK, zorder=3)
        ax.plot([j * 2.6 - 0.5, j * 2.6 + 1.5], [exp[idx]] * 2, color=INK, ls="--", lw=0.7)
    ax.set_xticks([0.5, 3.1], ["G:C at\nWRC/GYW", "A:T at\nWA/TW"])
    ax.set_ylabel("Share of calls (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Hotspot use", fontsize=7)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=SANG_C, label="Sanger"), Patch(color=ONT_C, label="Nanopore"),
                       plt.Line2D([], [], color=INK, ls="--", lw=0.7, label="sites")],
              fontsize=6, loc="upper left", ncol=3, bbox_to_anchor=(-0.05, 1.02), columnspacing=0.6, handlelength=1.2)
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    letter(ax, "D", dx=-0.3)
    # E: the same comparison mouse by mouse (crm/fig5_per_mouse.py, which checks it pools back to panel B)
    ax = fig.add_subplot(gs[1, 1:])
    pm = json.load(open(DATA / "fig5_per_mouse.json"))
    marks = {"WT2": "o", "WT3": "s", "WT4": "^"}
    for i, cls in enumerate(pm["classes"]):
        x0 = i * 3
        ax.axvspan(x0 - 0.6, x0 + 1.6, color="#f3f1ec", lw=0, zorder=0)
        for m in pm["mice"]:
            a, b = pm["share"]["sanger"][m][i], pm["share"]["ont"][m][i]
            ax.plot([x0, x0 + 1], [a, b], color=SOFT, lw=0.7, zorder=1)
            ax.scatter([x0], [a], s=13, marker=marks[m], color=SANG_C, edgecolor="white", lw=0.5, zorder=3)
            ax.scatter([x0 + 1], [b], s=13, marker=marks[m], color=ONT_C, edgecolor="white", lw=0.5, zorder=3)
        ax.text(x0 + 0.5, 47.5, cls.replace(" ", "\n", 1), ha="center", va="top", fontsize=6)
    ax.set_xticks([x for i in range(4) for x in (i * 3, i * 3 + 1)], ["S", "N"] * 4, fontsize=6)
    ax.set_xlim(-0.8, 10.8)
    ax.set_ylim(0, 48)
    for t in ax.get_xticklabels():
        t.set_color(SANG_C if t.get_text() == "S" else ONT_C)
    ax.set_ylabel("Share of calls (%)", labelpad=1)
    ax.set_title("The same comparison, mouse by mouse", fontsize=7)
    for m, mk in marks.items():
        ax.scatter([], [], s=12, marker=mk, color=MID, label=f"Wt {m[-1]}")
    ax.legend(loc="lower left", fontsize=6, ncol=3, handletextpad=0.1, columnspacing=0.8, frameon=False)
    ax.text(1.0, -0.2, "S, Sanger; N, nanopore. Platform × type interaction over 12 types, "
            f"p = {pm['summary']['rm_anova_interaction_p']:.2f}", transform=ax.transAxes, ha="right",
            fontsize=6, color=MID)
    letter(ax, "E", dx=-0.1)
    save(fig, "Figure5")


if __name__ == "__main__":
    main()
