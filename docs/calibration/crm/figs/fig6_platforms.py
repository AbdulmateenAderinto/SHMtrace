"""Figure 6: the same engine across platforms, loci and questions.

A  Illumina J_H4 (PRJNA973112): mutation frequency per mouse, Wt / Ldha-deleted / AID-deficient
B  CH12F3 Ung-/-Msh2-/-Aicda-/- (PRJNA918596): mutation frequency per replicate by condition
C  spectrum and targeting of AID and R174E/R178E calls
D  Ramos (PRJNA820976): deletions and insertions per clone
E  UMI consensus on R10 reads (Karst et al.): error per base by reads per molecule
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt

from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]
PUB = DATA / "public"


def strip(ax, groups, values, colors, floor):
    for i, (g, (f, e)) in enumerate(zip(groups, colors)):
        v = np.array(values[g], dtype=float)
        shown = np.where(v > 0, v, floor)
        x = i + np.linspace(-0.15, 0.15, len(v))
        for xx, yy, raw in zip(x, shown, v):
            ax.scatter(xx, yy, s=12, facecolor=f, edgecolor=e, lw=0.6, zorder=3, marker="o" if raw > 0 else "v")
        pos = v[v > 0]
        if len(pos):
            g_ = 10 ** np.mean(np.log10(np.where(v > 0, v, floor)))
            ax.plot([i - 0.28, i + 0.28], [g_, g_], color=e, lw=1.4)
    ax.axhline(floor * 1.5, color=SOFT, lw=0.4, ls=":")


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 112 * MM))
    gs = gridspec.GridSpec(2, 6, figure=fig, left=0.085, right=0.99, top=0.93, bottom=0.17, wspace=1.25, hspace=0.95)
    # A
    ax = fig.add_subplot(gs[0, 0:2])
    p = json.load(open(PUB / "ldha/project.json"))
    groups = ["WT", "LDHA_KO", "AID_KO"]
    vals = {g: [s["mutation_frequency_x1e3"] for s in p["samples"] if s["group"] == g] for g in groups}
    strip(ax, groups, vals, [(CASEF, CASE), ("#b9a3d6", "#6a4c93"), (CTRLF, CTRL)], 0.01)
    ax.set_yscale("log")
    ax.set_xticks(range(3), ["Wt", "$Ldha$ deleted", "AID-deficient"])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    ax.set_title("J$_H$4, Illumina (Sharma et al.)", fontsize=7)
    letter(ax, "A", dx=-0.32)
    # B
    ax = fig.add_subplot(gs[0, 2:4])
    p = json.load(open(PUB / "ch12/project.json"))
    groups = ["AID", "R50E", "R174ER178E", "EV"]
    vals = {g: [s["mutation_frequency_x1e3"] for s in p["samples"] if s["group"] == g] for g in groups}
    strip(ax, groups, vals, [(CASEF, CASE), ("#d9dde2", MID), ("#f3c78e", "#a3621b"), (CTRLF, CTRL)], 0.005)
    ax.set_yscale("log")
    ax.set_xticks(range(4), ["AID", "R50E", "R174E/\nR178E", "empty"])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    ax.set_title("CH12F3, Illumina (Wang et al.)", fontsize=7)
    letter(ax, "B", dx=-0.32)
    # C
    ax = fig.add_subplot(gs[0, 4:6])
    s5 = json.load(open(DATA / "targeting_ch12.json"))
    cats = [("C→T or G→A", ("CT", "GA"), "#2c5d8a"), ("other G:C", None, "#9fc0d8"), ("A:T", "AT", "#d9952b")]
    for i, g in enumerate(("AID", "R174ER178E")):
        vs = [v for v in p["variants"] if v["group"] == g]
        n = len(vs)
        ctga = sum(v["ref"] + v["alt"] in ("CT", "GA") for v in vs)
        gc_other = sum(v["ref"] in "GC" for v in vs) - ctga
        at = n - ctga - gc_other
        bottom = 0
        for (lab, _k, col), c in zip(cats, (ctga, gc_other, at)):
            ax.bar(i, 100 * c / n, 0.6, bottom=bottom, color=col, label=lab if i == 0 else None)
            bottom += 100 * c / n
        ax.text(i, 103, f"n = {n}\nS5F {s5[g]['index']:.2f}", ha="center", fontsize=5.5, color=MID)
    ax.set_xticks([0, 1], ["AID", "R174E/\nR178E"])
    ax.set_ylim(0, 125)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Share of calls (%)")
    ax.set_title("Spectrum and targeting", fontsize=7)
    ax.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, handlelength=1, columnspacing=0.7)
    letter(ax, "C", dx=-0.42)
    # D
    ax = fig.add_subplot(gs[1, 0:3])
    p = json.load(open(PUB / "hmces/project.json"))
    groups = [("WT", "Wt", CASEF, CASE), ("HMCES_KO", "HMCES-KO", "#b9a3d6", "#6a4c93"), ("AID_KO", "AID-KO", CTRLF, CTRL)]
    for j, kind in enumerate(("deletions", "insertions")):
        for i, (g, lab, f, e) in enumerate(groups):
            v = [s[kind] for s in p["samples"] if s["group"] == g]
            x = j * 3.6 + i
            ax.bar(x, np.mean(v), 0.75, color=f, edgecolor=e, lw=0.6)
            ax.scatter(x + np.linspace(-0.2, 0.2, len(v)), v, s=7, color=e, zorder=3)
    ax.set_xticks([0, 1, 2, 3.6, 4.6, 5.6], ["Wt", "HMCES\nKO", "AID\nKO"] * 2)
    ax.text(1, -0.36, "deletions", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, fontweight="bold")
    ax.text(4.6, -0.36, "insertions", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, fontweight="bold")
    ax.set_ylabel("Events per clone (20,000 reads)")
    ax.set_title("Indels in Ramos, Illumina (Wu et al.)", fontsize=7)
    letter(ax, "D", dx=-0.16)
    # E
    ax = fig.add_subplot(gs[1, 3:6])
    u = json.load(open(CAL / "umi_karst/results.json"))["summary"]
    bins = ["3-4", "5-9", "10-19", "20-+"]
    x = np.arange(len(bins))
    allerr = [u["by_reads_per_molecule"][b]["error_pct"] for b in bins]
    sub = [u["by_reads_per_molecule"][b]["substitution_error_pct"] for b in bins]
    ax.plot(x, allerr, marker="o", ms=3, color="#9aa5b1", lw=1.1, label="consensus, all errors")
    ax.plot(x, sub, marker="o", ms=3, color="#1f6f8b", lw=1.1, label="consensus, substitutions")
    ax.axhline(u["raw_read_error_pct"], color=CTRL, ls="--", lw=0.8)
    ax.text(3.0, u["raw_read_error_pct"] * 1.25, f"raw reads {u['raw_read_error_pct']:.1f}%", fontsize=5.5,
            color=CTRL, ha="right")
    ax.set_yscale("log")
    ax.set_ylim(0.005, 30)
    ax.set_xticks(x, ["3–4", "5–9", "10–19", "≥20"])
    ax.set_xlabel("Reads per molecule")
    ax.set_ylabel("Error per base (%)")
    ax.set_title(f"UMI consensus, R10 reads ({u['molecules_3plus_reads']} molecules)", fontsize=7)
    ax.legend(fontsize=5.3, loc="lower left")
    letter(ax, "E", dx=-0.16)
    save(fig, "Figure6")


if __name__ == "__main__":
    main()
