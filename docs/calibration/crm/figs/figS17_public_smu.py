"""Figure S17: public Sμ mutation libraries under PI3Kδ blockade (GSE77788).

A    what each library reads: the amplicon, the two read layouts, and the two analysis windows
B    unique mutations per library in the window every library reads
C    the same contrast over the Sμ repeats, using only the libraries sequenced the same way

Data from crm/public_smu_pi3k.py (crm/out/public_smu/public_smu.json).
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

TREAT = ["DMSO", "Idelalisib", "Duvelisib"]
WT_C, KO_C = "#274b7e", "#9a3439"


def main():
    setup()
    d = json.load(open(DATA / "public_smu/public_smu.json"))
    fig = plt.figure(figsize=(FULL, 112 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.08, right=0.98, top=0.93, bottom=0.11,
                           wspace=0.32, hspace=0.42, height_ratios=[0.62, 1])
    # A: the amplicon and what each read type reaches
    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(-330, 690)
    ax.set_ylim(0.0, 5.5)
    ax.axis("off")
    ax.plot([1, 651], [4.8, 4.8], color=INK, lw=2.2)
    ax.text(-8, 4.8, "S\u03bc amplicon, 651 bp\n(the bundled reference)", ha="right", va="center",
            fontsize=5.5, color=INK)
    ax.add_patch(plt.Rectangle((1, 4.3), 327, 0.3, color="#d9952b", alpha=0.4, lw=0))
    ax.text(164, 4.12, "S\u03bc repeats: where AID acts", ha="center", va="top", fontsize=5.3, color="#8a5a12")
    rows = [("300-bp single-end read\n(6 of 9 wild-type libraries)", 351, 650, 3.5, WT_C),
            ("R1 of a pair\n(all 9 AID-null, 3 wild-type)", 401, 650, 2.75, WT_C),
            ("R2 of the same pair", 1, 250, 2.0, KO_C)]
    for lab, lo, hi, y, c in rows:
        ax.plot([lo, hi], [y, y], color=c, lw=5, solid_capstyle="butt", alpha=0.85)
        ax.text(-8, y, lab, ha="right", va="center", fontsize=5.3, color=MID)
    ax.annotate("", xy=(399, 1.78), xytext=(252, 1.78),
                arrowprops=dict(arrowstyle="<->", color=SOFT, lw=0.8))
    ax.text(325, 1.68, "mates do not overlap: only 9\u201321% of pairs merge", ha="center", va="top",
            fontsize=5, color=SOFT)
    for lo, hi, y, lab, c in ((410, 625, 1.25, "window B: the stretch every library reads", INK),
                              (20, 240, 0.55, "window C: the repeats, paired libraries only", "#8a5a12")):
        ax.add_patch(plt.Rectangle((lo, y), hi - lo, 0.26, color=c, alpha=0.2, lw=0))
        ax.plot([lo, hi], [y + 0.13, y + 0.13], color=c, lw=1.3)
        ax.text(-8, y + 0.13, lab, ha="right", va="center", fontsize=5.3, color=c)
    letter(ax, "A", dx=-0.055, dy=0.99)   # upper-left of the schematic, not below it
    # B: the window every library reads
    ax = fig.add_subplot(gs[1, 0])
    for i, t in enumerate(TREAT):
        rows_t = d["treatments"][t]["samples"]
        for j, (grp, c) in enumerate((("WT", WT_C), ("AID_KO", KO_C))):
            v = [r["unique_mutations"] for r in rows_t if r["group"] == grp]
            x = i * 2.5 + j
            ax.scatter(x + np.linspace(-0.13, 0.13, len(v)), v, s=14, facecolor=c, edgecolor="white",
                       lw=0.4, zorder=3)
            ax.plot([x - 0.26, x + 0.26], [np.mean(v)] * 2, color=c, lw=1.4)
    ax.set_xticks([i * 2.5 + 0.5 for i in range(3)], TREAT, fontsize=5.5)
    ax.set_ylabel("Unique mutations per library")
    ax.set_title("Window B: the shared 216 bp", fontsize=7)
    ax.scatter([], [], s=12, color=WT_C, label="wild type")
    ax.scatter([], [], s=12, color=KO_C, label="$\\it{Aicda}^{-/-}$")
    ax.legend(fontsize=5.2, loc="upper left", frameon=False, handletextpad=0.1)
    ax.text(0.98, 0.78, "AID-null at background throughout;\nreplicates vary more\nthan treatments do",
            transform=ax.transAxes, ha="right", va="top", fontsize=4.9, color=MID)
    letter(ax, "B", dx=-0.17)
    # C: the repeats, paired libraries only
    ax = fig.add_subplot(gs[1, 1])
    rp = d["repeats_paired"]["samples"]
    for j, (grp, c, lab) in enumerate((("WT", WT_C, "wild type\n(n = 3)"),
                                       ("AID_KO", KO_C, "$\\it{Aicda}^{-/-}$\n(n = 9)"))):
        v = [r["unique_mutations"] for r in rp if r["group"] == grp]
        ax.scatter(np.full(len(v), j) + np.linspace(-0.16, 0.16, len(v)), v, s=14, facecolor=c,
                   edgecolor="white", lw=0.4, zorder=3)
        ax.plot([j - 0.28, j + 0.28], [np.mean(v)] * 2, color=c, lw=1.4)
        ax.text(j, -14, lab, ha="center", va="top", fontsize=5.3, color=c)
    s5 = [r["s5f_index"] for r in rp if r["group"] == "WT" and r["s5f_index"]]
    ct = [r["pct_CT_GA"] for r in rp if r["group"] == "WT" and r["pct_CT_GA"]]
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([])
    ax.set_ylabel("Unique mutations per library")
    ax.set_title("Window C: the Sμ repeats", fontsize=7)
    ax.text(0.5, 0.60, f"calls carry the AID signature:\ntargeting index {min(s5):.2f}–{max(s5):.2f}, "
                       f"C→T/G→A {min(ct):.0f}–{max(ct):.0f}%",
            transform=ax.transAxes, ha="center", va="top", fontsize=5, color=MID)
    # the groups are not balanced: one wild-type library per treatment against three AID-null, so
    # treatment and replicate are confounded on the wild-type side (minor 16)
    ax.text(0.5, -0.24, "every library here is paired, so no layout difference; but one wild-type\n"
                        "library per treatment against three AID-null, pooled as a descriptive contrast",
            transform=ax.transAxes, ha="center", fontsize=5, color=SOFT)
    letter(ax, "C", dx=-0.17)
    save(fig, "FigureS17_public_smu")


if __name__ == "__main__":
    main()
