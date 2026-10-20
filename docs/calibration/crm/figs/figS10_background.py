"""Figure S10: where the background sits along the amplicon, and what the reads behind a call look like.

The controls agree only weakly position by position (Spearman r = 0.12-0.45, so one control explains a few
per cent of another's variation), which is why the error model pools over substitution type, strand and base
quality and shrinks each position towards that pooled rate rather than trusting a position's own count.

A-B  mismatch rate (base quality >= 20, germline sites excluded) at every position of the window in each
     AID-null control and in the AID-active libraries (mean), J_H4 and Smu; r, Spearman between controls
C-D  read rasters, J_H4: every mismatch of quality >= 20 in each read of one wild-type and one Aicda-/-
     library; reads sorted by clone and mutation load; called sites in color
Data: out/intuition_<dataset>.json (intuition_data.py).
"""
import json
from itertools import combinations

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from style import CASE, CTRL, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CTRL_SHADES = ("#9a3439", "#c4585e", "#e0999c")


def rates(d, n, lo, hi):
    prof = d["profile"][n]
    return np.array([prof[str(x)][0] / prof[str(x)][1] if str(x) in prof and prof[str(x)][1] else np.nan
                     for x in range(lo, hi + 1)])


def profile_panel(ax, key, title):
    d = json.load(open(DATA / f"intuition_{key}.json"))
    lo, hi = d["window"]
    xs = np.arange(lo, hi + 1)
    ctl = {c: 100 * rates(d, c, lo, hi) for c in d["controls"]}
    case = 100 * np.nanmean([rates(d, n, lo, hi) for n in d["cases"]], axis=0)
    ax.plot(xs, case, color=CASE, lw=0.6, alpha=0.55, label="AID-active (mean)")
    for (c, r), col in zip(ctl.items(), CTRL_SHADES, strict=False):
        ax.plot(xs, r, color=col, lw=0.7, label=f"{c} (AID-null)")
    rs = [spearmanr(ctl[a], ctl[b], nan_policy="omit")[0] for a, b in combinations(ctl, 2)]
    # a Spearman coefficient is a rank statistic; squaring it is not a share of variance, so the panel
    # reports the correlations themselves (minor 18)
    print(f"   {title}: between-control Spearman r = {min(rs):.2f}-{max(rs):.2f} ({len(rs)} pairs)")
    ax.text(0.99, 0.95, f"between controls: Spearman r = {min(rs):.2f}–{max(rs):.2f}\n"
            "(weak rank agreement position by position)",
            transform=ax.transAxes, ha="right", va="top", fontsize=6, color=INK, linespacing=1.3)
    # the controls on their own scale: the AID-active curve is an order of magnitude higher and hides them
    ins = ax.inset_axes([0.055, 0.40, 0.42, 0.32], facecolor="white", zorder=5)
    ins.patch.set_alpha(1.0)
    for r, col in zip(ctl.values(), CTRL_SHADES, strict=False):
        ins.plot(xs, r, color=col, lw=0.5)
    ins.set_xlim(lo, hi)
    ins.set_ylim(0, np.nanpercentile(np.concatenate(list(ctl.values())), 99.9) * 1.15)
    ins.tick_params(labelsize=6, length=1.5, pad=1)
    # a white backing: at 6 pt the label and the legend would otherwise sit on the busy traces
    ins.set_title("controls only, own scale", fontsize=6, color=MID, pad=1.5,
                  bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    for sp in ins.spines.values():
        sp.set_linewidth(0.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, np.nanpercentile(np.concatenate([*ctl.values(), case]), 99.7) * 1.25)
    ax.set_xlabel(f"{title} position")
    ax.set_ylabel("Mismatches, Q ≥ 20 (%)")
    ax.set_title(f"{title}: background along the amplicon ({d['depth']:,} reads per library)", fontsize=6.8)
    ax.legend(fontsize=6, loc="upper left", ncol=4, handlelength=1.2, borderaxespad=0.2, columnspacing=1,
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.9)


def raster_panel(ax, d, n, label):
    lo, hi = d["window"]
    reads = d["raster"][n]
    sizes = {}
    for r in reads:
        if r["clone"] is not None:
            sizes[r["clone"]] = sizes.get(r["clone"], 0) + 1
    order = sorted(range(len(reads)), key=lambda i: (-sizes.get(reads[i]["clone"], 0), reads[i]["clone"] or -1,
                                                     -sum(m[2] for m in reads[i]["muts"]), -len(reads[i]["muts"])))
    gx, gy, cx, cy = [], [], [], []
    for y, i in enumerate(order):
        for pos, _alt, called in reads[i]["muts"]:
            (cx if called else gx).append(pos)
            (cy if called else gy).append(y)
    ax.scatter(gx, gy, s=0.6, color=SOFT, linewidths=0, rasterized=True)
    ax.scatter(cx, cy, s=1.2, color=CASE if label == "wild type" else CTRL, linewidths=0, rasterized=True)
    # clone boundaries for the largest clones
    y = 0
    for _c, k in sorted(sizes.items(), key=lambda x: -x[1])[:5]:
        if k < 10:
            break
        y += k
        ax.axhline(y - 0.5, color=MID, lw=0.5, ls=":")
    ax.set_xlim(lo, hi)
    ax.set_ylim(len(reads), 0)
    ax.set_xlabel("J$_H$4 position")
    ax.set_ylabel("Reads (sorted by clone)")
    n_called = sum(m[2] for r in reads for m in r["muts"])
    n_other = sum(1 - m[2] for r in reads for m in r["muts"])
    ax.set_title(f"{n} ({label}): {n_called:,} at called sites, {n_other:,} other mismatches", fontsize=6.8)


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 150 * MM))
    gs = gridspec.GridSpec(3, 2, figure=fig, left=0.075, right=0.985, top=0.96, bottom=0.06, hspace=0.62,
                           wspace=0.22, height_ratios=[1, 1, 1.7])
    ax = fig.add_subplot(gs[0, :])
    profile_panel(ax, "jh4", "J$_H$4")
    letter(ax, "A", dx=-0.06)
    ax = fig.add_subplot(gs[1, :])
    profile_panel(ax, "su_lps", "Sμ")
    letter(ax, "B", dx=-0.06)
    d = json.load(open(DATA / "intuition_jh4.json"))
    for j, (n, label) in enumerate((("WT3", "wild type"), ("KO1", "Aicda$^{-/-}$"))):
        ax = fig.add_subplot(gs[2, j])
        raster_panel(ax, d, n, label)
        letter(ax, "CD"[j], dx=-0.14)
    save(fig, "FigureS10")


if __name__ == "__main__":
    main()
