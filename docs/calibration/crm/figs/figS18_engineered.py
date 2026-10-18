"""Figure S18: engineered AID diversification, four public datasets.

A    AID against an APOBEC on one GFP substrate: what separates them and what does not (SRP159221)
B    the mutation position follows the guide, five sgRNAs on one amplicon (SRP075898)
C    localised sites against the number of guides carried (SRP075898)
D    an AID-diversified TCR against its own ancestral construct (figshare 12803501)

Only completed analyses are drawn. Panels B and C are per-position mismatch rates, not mutation calls, and
say so; A and D are calls from the project workflow. Data from crm/liu_mutators.py, crm/tam_footprint.py,
crm/tam_saturation.py and crm/tcr_bassan.py.
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

AID_C, APO_C, BASE_C = "#274b7e", "#a3621b", "#9a3439"
GUIDE_C = {"sgRNA_4": "#8ab0d0", "sgRNA_6": "#d9952b", "sgRNA_7": "#b4bec8",
           "sgRNA_9": "#9a3439", "sgRNA_15": "#1f6f8b"}


def main():
    setup()
    liu = json.load(open(DATA / "liu/liu_mutators.json"))["samples"]
    foot = json.load(open(DATA / "tam/footprint.json"))
    sat = json.load(open(DATA / "tam/saturation.json"))["totals"]
    tcr = json.load(open(DATA / "tcr/tcr_bassan.json"))
    fig = plt.figure(figsize=(FULL, 108 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.08, right=0.98, top=0.92, bottom=0.11,
                           wspace=0.33, hspace=0.62)
    # A: what separates the two deaminases
    ax = fig.add_subplot(gs[0, 0])
    groups = [("wild-type 293T", "AID_wt293T", "A3A_wt293T"), ("UNG/MSH2-null", "AID_ungmsh2", "A3A_ungmsh2")]
    for gi, (_lab, a, b) in enumerate(groups):
        for mi, (key, mark) in enumerate((("pct_GC_in_WRC_GYW", "o"), ("pct_CT_GA", "s"))):
            for ci, (pfx, col) in enumerate(((a, AID_C), (b, APO_C))):
                vals = [v[key] for k, v in liu.items() if k.startswith(pfx) and v[key] is not None]
                x = gi * 2.6 + mi * 1.1 + ci * 0.34
                ax.scatter([x] * len(vals), vals, s=11, marker=mark, facecolor=col, edgecolor="white",
                           lw=0.3, zorder=3)
                ax.plot([x - 0.14, x + 0.14], [np.mean(vals)] * 2, color=col, lw=1.2)
    ax.set_xticks([0.17, 1.27, 2.77, 3.87], ["WRC/\nGYW", "C→T/\nG→A", "WRC/\nGYW", "C→T/\nG→A"], fontsize=5)
    for x, lab in ((0.72, "wild-type"), (3.32, "UNG/MSH2-null")):
        ax.text(x, -0.30, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=5.5, color=MID)
    ax.set_ylabel("Share of calls (%)")
    ax.set_ylim(0, 100)
    ax.set_title("One substrate, two deaminases", fontsize=7)
    ax.scatter([], [], s=11, color=AID_C, label="AID")
    ax.scatter([], [], s=11, color=APO_C, label="APOBEC3A")
    ax.legend(fontsize=5, loc="upper left", frameon=False, handletextpad=0.1, ncol=2, columnspacing=0.7)
    ax.text(0.98, 0.04, "the motif separates them in both;\nthe transition share only in wild type",
            transform=ax.transAxes, ha="right", fontsize=4.8, color=MID)
    letter(ax, "A", dx=-0.2)
    # B: the position follows the guide, one track per guide so the peaks do not overlap
    ax = fig.add_subplot(gs[0, 1])
    prof = foot["amplicons"]["amp1_GTCTTGTTGG"]["profiles_per_1e3"]
    order = ["sgRNA_4", "sgRNA_6", "sgRNA_7", "sgRNA_9", "sgRNA_15"]
    # each track is scaled to its own maximum: the panel is about where the excess sits, not how big it is,
    # and a shared scale makes the smallest guide's peak invisible
    for row, g in enumerate(order):
        p_ = prof[g]
        xs = [k for k, r in enumerate(p_) if r is not None]
        top = max(p_[k] for k in xs)
        ys = [row + 0.86 * p_[k] / top for k in xs]
        ax.fill_between(xs, row, ys, color=GUIDE_C[g], lw=0, alpha=0.9)
        ax.axhline(row, color=SOFT, lw=0.4)
        ax.text(-10, row + 0.40, g.replace("sgRNA_", "#"), ha="right", va="center", fontsize=5.5,
                color=GUIDE_C[g])
        peak = max(xs, key=lambda k: p_[k])
        ax.text(peak + 5, row + 0.62, f"{top:.0f}", fontsize=4.6, color=GUIDE_C[g], va="center")
    ax.set_xlim(-42, 245)
    ax.set_ylim(0, 5.35)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)      # no y scale here, and it collided with the track labels
    ax.set_xlabel("Position in the amplicon (bp)")
    ax.set_title("The position follows the guide", fontsize=7)
    ax.text(0.5, -0.30, "five sgRNAs on one amplicon; each track scaled to its own peak, whose\n"
                        "height in mismatches per 10$^3$ is given. Rates, not calls.",
            transform=ax.transAxes, ha="center", va="top", fontsize=4.8, color=MID)
    letter(ax, "B", dx=-0.1)
    # C: sites against guides
    ax = fig.add_subplot(gs[1, 0])
    ns = [int(k.split("_")[0]) for k in sat]
    vals = [sat[k] for k in sat]
    order = np.argsort(ns)
    ax.plot([ns[i] for i in order], [vals[i] for i in order], color=INK, lw=1.1, marker="o", ms=4,
            zorder=3)
    for i in order:
        ax.annotate(str(vals[i]), (ns[i], vals[i]), xytext=(0, 6), textcoords="offset points",
                    ha="center", fontsize=5, color=MID)
    ax.set_xticks(sorted(ns))
    ax.set_xlabel("sgRNAs carried by the library")
    ax.set_ylabel("Localised sites")
    ax.set_ylim(-6, 80)
    ax.set_title("More guides, more sites", fontsize=7)
    ax.text(0.03, 0.88, "a site sits ≥ 5× above the lowest of\nthe other libraries at that base",
            transform=ax.transAxes, fontsize=4.8, color=MID)
    letter(ax, "C", dx=-0.2)
    # D: the TCR against its own ancestor
    ax = fig.add_subplot(gs[1, 1])
    s = tcr["samples"]
    base = "h5_wt"
    others = [k for k in s if k != base]
    ax.scatter([0] * 1, [s[base]["mutation_frequency_x1e3"]], s=22, marker="D", facecolor=BASE_C,
               edgecolor="white", lw=0.4, zorder=4)
    xs = np.linspace(0.75, 1.25, len(others))
    ax.scatter(xs, [s[k]["mutation_frequency_x1e3"] for k in others], s=14, facecolor=AID_C,
               edgecolor="white", lw=0.3, zorder=3)
    ax.plot([0.72, 1.28], [np.mean([s[k]["mutation_frequency_x1e3"] for k in others])] * 2,
            color=AID_C, lw=1.2)
    ax.set_xticks([0, 1], ["the available\nbaseline", "diversified and\nselected (n = 10)"], fontsize=5.3)
    ax.set_xlim(-0.45, 1.6)
    ax.set_ylabel("Mutation frequency per 10$^3$")
    ax.set_title("An AID-diversified TCR", fontsize=7)
    ax.text(0.98, 0.04, "one baseline library, so no\nfalse-discovery rate is claimed",
            transform=ax.transAxes, ha="right", fontsize=4.8, color=SOFT)
    letter(ax, "D", dx=-0.22)
    save(fig, "FigureS18_engineered")


if __name__ == "__main__":
    main()
