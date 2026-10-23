"""Figure S9: what depth buys, with the draw-to-draw spread shown.

A-B  J_H4: unique variants and density against depth, every draw drawn
C-D  N114: the same
E    the share of each AID-active library's own variant set recovered at every depth, not only at 500 reads

The earlier version stacked two figures, each lettered from A, and drew one random draw per depth, so the
reader could not tell a depth effect from the luck of the draw. Here the whole analysis is repeated five
times at each depth (crm/depth_draws.py), controls redrawn and recalibrated each time.

Zeros are kept as zeros: an AID-null library with no calls has a density of 0, which has no logarithm. Rather
than floor it silently, the log panels draw it on a separate "0" row below the axis break and the legend says
so. Data: out/depth_draws_{jh4,n114}.json.
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, MID, MM, SOFT, letter, save, setup

ASSAYS = {"jh4": ("J$_H$4 intron", ("WT2", "WT3", "WT4"), ("KO1", "KO2", "KO3"), "Wt", "$Aicda^{-/-}$"),
          "n114": ("N114 V$_H$", ("C1", "D3", "D5"), ("C3", "C4", "D1"), "pMX-$Aicda$", "empty vector")}
LINTHRESH = 1.0      # below this the axis is linear, so a zero is drawn at zero and not at an invented floor


def load(key):
    runs = json.load(open(DATA / f"depth_draws_{key}.json"))
    per = {}
    for r in runs["runs"]:
        per.setdefault((r["sample"], r["depth"]), []).append(r)
    return runs, per


def series(ax, per, samples, depths, field, colr, colf):
    """Each draw as a faint point, the mean over draws as a line. Zeros are drawn at zero (symlog axis)."""
    xs = np.arange(len(depths))
    for s in samples:
        means = []
        for i, d in enumerate(depths):
            v = [r[field] for r in per[(s, d)]]
            means.append(np.mean(v))
            ax.scatter(np.full(len(v), i) + np.linspace(-0.13, 0.13, len(v)), v, s=3.2, color=colf,
                       linewidths=0, zorder=2)
        ax.plot(xs, means, color=colr, lw=1.1, marker="o", ms=2.6, zorder=3)
    ax.set_xticks(xs, [str(d) for d in depths], fontsize=6)


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 104 * MM))
    gs = gridspec.GridSpec(2, 3, figure=fig, left=0.075, right=0.985, top=0.9, bottom=0.11, wspace=0.42,
                           hspace=0.62, width_ratios=[1, 1, 1])
    reached = {}
    for row, key in enumerate(("jh4", "n114")):
        title, case, ctl, case_lab, ctl_lab = ASSAYS[key]
        runs, per = load(key)
        depths = runs["depths"]
        for col, (field, ylab, short) in enumerate((("calls", "Unique variants", "unique variants"),
                                                    ("density", "Unique variants per 10$^5$ bases", "density"))):
            ax = fig.add_subplot(gs[row, col])
            series(ax, per, case, depths, field, CASE, CASEF)
            series(ax, per, ctl, depths, field, CTRL, CTRLF)
            # symlog: linear below 1, so a library with no calls sits at zero rather than at an invented floor
            ax.set_yscale("symlog", linthresh=LINTHRESH, linscale=0.45)
            ax.set_ylim(0, None)
            ax.axhline(LINTHRESH, color=SOFT, lw=0.5, ls=":")
            ax.set_xlabel("Reads analyzed")
            ax.set_ylabel(ylab)
            ax.set_title(f"{title}: {short}", fontsize=6.6)
            letter(ax, "ABCD"[row * 2 + col], dx=-0.3)
        # the share of that library's own variant set recovered at each depth. The denominator is the most
        # that library ever yielded, over every depth and draw, so the curve is read against the library's own
        # discovery rather than against the 1,000-read value alone.
        for s_ in case:
            mx = max(r["calls"] for d in depths for r in per[(s_, d)])
            reached[s_] = ({d: [100 * r["calls"] / mx for r in per[(s_, d)]] for d in depths}, key)
    # E
    ax = fig.add_subplot(gs[:, 2])
    depths = load("jh4")[0]["depths"]
    xs = range(len(depths))
    for _s, (by_depth, key) in reached.items():
        col = CASE if key == "jh4" else "#6d9ad0"
        mean = [np.mean(by_depth[d]) for d in depths]
        ax.plot(xs, mean, color=col, lw=0.9, marker="o" if key == "jh4" else "s", ms=2.6,
                markeredgewidth=0, zorder=3)
        ax.vlines(xs, [min(by_depth[d]) for d in depths], [max(by_depth[d]) for d in depths],
                  color=col, lw=0.6, alpha=0.7)
    ax.axhline(100, color=SOFT, lw=0.5, ls=":")
    ax.axvline(depths.index(500), color=SOFT, lw=0.5)
    ax.set_xticks(xs, [str(d) for d in depths], fontsize=6)
    ax.set_xlabel("Reads analyzed")
    ax.set_ylim(0, 112)
    ax.set_ylabel("Share of that library's variants recovered (%)")
    ax.set_title("Share of each library's own total", fontsize=6.6)
    ax.text(0.97, 0.06, "AID-active libraries only\nJ$_H$4 (dark), N114 V$_H$ (light)\nlines, the range of five draws",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6, color=MID, linespacing=1.4)
    letter(ax, "E", dx=-0.26)
    save(fig, "FigureS10")
    out = {s_: {"mean_pct": {str(d): round(float(np.mean(v[d])), 1) for d in v},
                "mean_pct_at_500": round(float(np.mean(v[500])), 1)}
           for s_, (v, _k) in reached.items()}
    print({s: out[s]["mean_pct_at_500"] for s in out})


if __name__ == "__main__":
    main()
