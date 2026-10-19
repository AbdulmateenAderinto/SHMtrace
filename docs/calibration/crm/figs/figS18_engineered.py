"""Figure S18: engineered AID diversification.

A    AID against an APOBEC on one GFP substrate: what separates them and what does not (SRP159221)
B    an AID-diversified antibody V_H against the authors' own list of its mutations (Devilder et al.)

A is calls from the project workflow against replicated no-editor controls. B is the only engineered
analysis with an external ground truth: the comparator is the authors' published mutation list, not
another run of this software. Data from crm/liu_mutators.py and crm/devilder_clones.py.
"""
import json
from collections import Counter

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
    dev = json.load(open(DATA / "devilder/devilder_clones.json"))["rounds"]
    fig = plt.figure(figsize=(FULL, 62 * MM))
    gs = gridspec.GridSpec(1, 2, figure=fig, left=0.085, right=0.975, top=0.87, bottom=0.27,
                           wspace=0.30)
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
        ax.text(x, -0.26, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=5.5, color=MID)
    ax.set_ylabel("Share of calls (%)")
    ax.set_ylim(0, 100)
    ax.set_title("One substrate, two deaminases", fontsize=7)
    ax.scatter([], [], s=11, color=AID_C, label="AID")
    ax.scatter([], [], s=11, color=APO_C, label="APOBEC3A")
    ax.legend(fontsize=5, loc="upper left", frameon=False, handletextpad=0.1, ncol=2, columnspacing=0.7)
    ax.text(0.98, 0.04, "the motif separates them in both;\nthe transition share only in wild type",
            transform=ax.transAxes, ha="right", fontsize=4.8, color=MID)
    letter(ax, "A", dx=-0.2)
    # B: every mutation the authors list, against what SHMtrace called in the same cluster. The
    # comparator is their published list, so this is the one engineered panel with outside truth.
    ax = fig.add_subplot(gs[0, 1])
    pts, rounds = Counter(), []
    for rnd, clusters in dev.items():
        n_c = n_m = 0
        for c in clusters:
            if not c["published_mutations"]:
                continue
            pts[(len(c["published_mutations"]), len(c["recovered_substitutions"]), rnd)] += 1
            n_c += 1
            n_m += len(c["published_mutations"])
        rounds.append((rnd, n_c, n_m))
    hi = max(max(a, b) for a, b, _r in pts) + 1
    ax.plot([0, hi], [0, hi], color=SOFT, lw=0.8, ls="--", zorder=1)
    for (x, y, rnd), n in sorted(pts.items()):
        col, mark = (AID_C, "o") if rnd == "R1" else (BASE_C, "s")
        ax.scatter([x], [y], s=9 + 7 * n, marker=mark, facecolor=col, edgecolor="white", lw=0.4, zorder=3)
    for rnd, n_c, n_m in sorted(rounds):
        col, mark = (AID_C, "o") if rnd == "R1" else (BASE_C, "s")
        ax.scatter([], [], s=14, marker=mark, color=col,
                   label=f"round {rnd[-1]}: {n_c} entries, {n_m} mutations")
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xlabel("Mutations the authors list")
    ax.set_ylabel("Mutations SHMtrace called")
    ax.set_title("An external list, cluster by cluster", fontsize=7)
    ax.legend(fontsize=4.8, loc="upper left", frameon=False, handletextpad=0.1)
    tot_c = sum(n for _r, n, _m in rounds)
    tot_m = sum(m for _r, _n, m in rounds)
    ax.text(0.97, 0.06, f"{tot_m}/{tot_m} recovered, none added,\nover {tot_c} cluster entries; "
                        "marker area, entries", transform=ax.transAxes, ha="right", fontsize=4.8, color=MID)
    letter(ax, "B", dx=-0.22)

    save(fig, "FigureS18_engineered")


if __name__ == "__main__":
    main()
