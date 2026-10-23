"""Figure S19: class-switch junctions in a public PacBio dataset (PRJNA831666).

A    what each library yields: reads, junction reads and unique junctions per sample
B    microhomology classes of the unique junctions, per genotype
C    mean microhomology per sample, with the wild-type contrast
D    acceptor-side breakpoint ambiguity measured directly, against the module's own flag

Data from crm/switch_pacbio.py (crm/out/switch_pacbio/switch_pacbio.json).
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import CASE, FULL, INK, MID, MM, SOFT, letter, save, setup

GEN = ["WT", "Shld1", "Xlf", "Shld1_Xlf"]
LAB = {"WT": "Wild type", "Shld1": "$\\it{Shld1}^{{-}/{-}}$", "Xlf": "$\\it{Xlf}^{{-}/{-}}$",
       "Shld1_Xlf": "$\\it{Shld1}^{{-}/{-}}$\n$\\it{Xlf}^{{-}/{-}}$"}
COL = {"WT": "#274b7e", "Shld1": "#1f6f8b", "Xlf": "#a3621b", "Shld1_Xlf": "#9a3439"}
CLS = ["blunt", "1-3 nt", ">=4 nt", "insertion"]
CLSC = {"blunt": "#cfd8df", "1-3 nt": "#7fa7cf", ">=4 nt": "#2c5d8a", "insertion": "#d9952b"}


def main():
    setup()
    d = json.load(open(__import__("style").DATA / "switch_pacbio/switch_pacbio.json"))
    S, G = d["samples"], d["groups"]
    by = {g: [v for k, v in S.items() if G[k] == g] for g in GEN}
    comp = {(r["group"], r["feature"]): r for r in d["comparisons"]}
    fig = plt.figure(figsize=(FULL, 104 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.105, right=0.98, top=0.9, bottom=0.12,
                           wspace=0.33, hspace=0.75)
    # A: what each library yields
    ax = fig.add_subplot(gs[0, 0])
    x = 0
    for g in GEN:
        for v in by[g]:
            ax.plot([x, x, x], [v["reads"], v["junction_reads"], v["unique_junctions"]], color=SOFT, lw=0.6,
                    zorder=1)
            ax.scatter([x], [v["reads"]], s=11, facecolor="white", edgecolor=INK, lw=0.7, zorder=3)
            ax.scatter([x], [v["junction_reads"]], s=11, color=SOFT, zorder=3)
            ax.scatter([x], [v["unique_junctions"]], s=13, color=COL[g], zorder=4)
            x += 1
        x += 0.8
    ax.set_yscale("log")
    ax.set_ylabel("Reads or junctions")
    ax.set_xticks([1 + i * 3.8 for i in range(4)], [LAB[g] for g in GEN], fontsize=6)
    ax.set_ylim(top=6e4)                 # headroom: the legend sits above the points, not on them
    ax.set_title("Each library's yield", fontsize=7)
    ax.scatter([], [], s=11, facecolor="white", edgecolor=INK, lw=0.7, label="reads")
    ax.scatter([], [], s=11, color=SOFT, label="junction reads")
    ax.scatter([], [], s=13, color=INK, label="unique junctions")
    ax.legend(fontsize=6, loc="upper left", ncol=3, frameon=False, handletextpad=0.1, columnspacing=0.7,
              bbox_to_anchor=(-0.02, 1.02))
    ax.text(0.5, -0.38, "reads per unique junction range from 1.1 to 23 across samples",
            transform=ax.transAxes, ha="center", fontsize=6, color=MID)
    letter(ax, "A", dx=-0.22)
    # B: microhomology classes, one bar per biological sample (ported from the Codex v2 revision: a group-mean
    # stack hid the between-mouse spread that the tests in C rest on)
    ax = fig.add_subplot(gs[0, 1])
    x, ticks = 0, []
    for g in GEN:
        ticks.append(x + (len(by[g]) - 1) / 2)
        for v in by[g]:
            bottom = 0
            for c in CLS:
                val = v["microhomology_classes_pct"].get(c, 0.0)
                ax.bar(x, val, 0.8, bottom=bottom, color=CLSC[c], edgecolor="white", lw=0.5)
                bottom += val
            x += 1
        x += 0.8
    ax.set_xticks(ticks, [LAB[g] for g in GEN], fontsize=6)
    ax.set_ylabel("Each mouse's unique junctions (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Junction structure, mouse by mouse", fontsize=7, pad=14)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=CLSC[c], label=c.replace(">=", "\u2265")) for c in CLS], fontsize=6, ncol=4,
              loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1, columnspacing=0.8, borderaxespad=0.2)
    letter(ax, "B", dx=-0.22)
    # C: mean microhomology per sample
    ax = fig.add_subplot(gs[1, 0])
    for i, g in enumerate(GEN):
        vals = [v["mean_microhomology"] for v in by[g]]
        ax.scatter(np.full(len(vals), i) + np.linspace(-0.11, 0.11, len(vals)), vals, s=15,
                   facecolor=COL[g], edgecolor="white", lw=0.5, zorder=3)
        ax.plot([i - 0.25, i + 0.25], [np.mean(vals)] * 2, color=COL[g], lw=1.4)
        if g != "WT":
            p = comp[(g, "mean_microhomology")]["p_welch"]
            txt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
            ax.text(i, 2.02, txt, ha="center", fontsize=6, color=MID)
    ax.set_xticks(range(4), [LAB[g] for g in GEN], fontsize=6)
    ax.set_ylabel("Mean microhomology (nt)")
    ax.set_ylim(0, 2.2)
    ax.set_title("Microhomology use, by sample", fontsize=7)
    ax.text(0.02, 0.04, "Welch t against wild type, n = 3 mice per genotype;\n"
                        "the permutation floor at n = 3 is p = 0.1",
            transform=ax.transAxes, fontsize=6, color=MID)
    letter(ax, "C", dx=-0.22)
    # D: breakpoint ambiguity, measured directly against the module's own flag
    ax = fig.add_subplot(gs[1, 1])
    for i, g in enumerate(GEN):
        vals = [v["acceptor_ambiguity_direct"]["ambiguous_pct"] for v in by[g]]
        ax.scatter(np.full(len(vals), i) + np.linspace(-0.11, 0.11, len(vals)), vals, s=15,
                   facecolor=COL[g], edgecolor="white", lw=0.5, zorder=3)
        ax.plot([i - 0.25, i + 0.25], [np.mean(vals)] * 2, color=COL[g], lw=1.4)
    ax.axhline(0, color=CASE, lw=1.0, ls="--")
    ax.text(3.45, 1.5, "the module's own flag:\n0% in every sample", ha="right", fontsize=6, color=CASE)
    ax.set_xticks(range(4), [LAB[g] for g in GEN], fontsize=6)
    ax.set_ylabel("Unique junctions (%)")
    ax.set_ylim(-3, 55)
    ax.set_title("Acceptor breakpoint ambiguity", fontsize=7)
    letter(ax, "D", dx=-0.22)
    save(fig, "FigureS20_switch_pacbio")


if __name__ == "__main__":
    main()
