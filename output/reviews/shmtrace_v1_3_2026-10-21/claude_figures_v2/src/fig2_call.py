"""Figure 2: how a call is made, and where real SHM sits.

A    a true AID hotspot in one AID-active library: share of reads carrying it, by strand, against each control
B    a recurrent background site in the same library, which a read-level test calls
C    the base qualities of the supporting reads at each site, aligned under A and B
D    the control false-discovery threshold for the same library: candidates kept against the false calls the
     held-out controls predict
E    the allele fractions of SHMtrace's calls, all AID-active libraries of four loci, against the range in which
     (the somatic-caller comparison that used to sit here is in Figure 3C)
Data: out/intuition_<dataset>.json (intuition_data.py).
"""
import json
from bisect import bisect_right

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import SCRATCH  # noqa: E402
from style import (
    CALLERS,
    CASE,
    CASEF,
    CTRL,
    CTRLF,
    DATA,
    FULL,
    INK,
    MID,
    MM,
    SOFT,
    letter,
    save,
    setup,
)

SITE_SET = "su_lps"
LOCI = {"jh4": "J$_H$4", "n114": "N114 V$_H$", "su_lps": "Sμ", "pim1": "$Pim1$"}


def load(key):
    return json.load(open(DATA / f"intuition_{key}.json"))


def sci(p):
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}×10$^{{{int(e)}}}$"


def site_panel(ax, site, case, show_y=True):
    """Share of reads carrying the change in the case and each control, split by strand. A and B are drawn by
    this one function on the same axis range, so the two sites read as a matched pair."""
    libs = list(site["libraries"])
    ys = np.arange(len(libs))[::-1]
    for y, n in zip(ys, libs, strict=True):
        v = site["libraries"][n]
        d = v["depth_plus"] + v["depth_minus"]
        plus, minus = 100 * v["alt_plus"] / d, 100 * v["alt_minus"] / d
        dark, light = (CASE, CASEF) if n == case else (CTRL, CTRLF)
        ax.barh(y, plus, 0.62, color=dark)
        ax.barh(y, minus, 0.62, left=plus, color=light)
        # Bars count all qualities; Q20 counts are a descriptive subset, not a test cutoff.
        hq = sum(q >= 20 for q in v["alt_quals"])
        ax.text(plus + minus + 0.04, y, f"{v['alt_plus'] + v['alt_minus']}/{d} ({hq} at Q≥20)", va="center",
                fontsize=6, color=MID)
    ax.set_yticks(ys, [f"{n} (case)" if n == case else n for n in libs] if show_y else [], fontsize=6)
    ax.set_ylim(-0.6, len(libs) - 0.15)
    ax.set_xlim(0, 1.95)
    ax.set_xticks([0, 0.5, 1, 1.5])
    ax.set_xlabel("Reads carrying the change (%)")
    called = site["called"]
    col = CASE if called else CTRL
    ax.text(0.98, 0.97, f"{'called' if called else 'not called'}\np = {sci(site['p'])}", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5, color=col, fontweight="bold", linespacing=1.25)


def site_title(ax, site, label):
    ctx = site["context"]
    ax.set_title(f"{label}\nposition {site['pos']}, {site['ref']}>{site['alt']} in "
                 f"{ctx[:2]}$\\underline{{{ctx[2]}}}${ctx[3:]}", fontsize=6.6, linespacing=1.3)


def quality_strip(ax, site, case, show_y=True):
    """Base qualities of the reads supporting the change: the case against the controls pooled."""
    rng = np.random.default_rng(0)
    groups = [(case, site["libraries"][case]["alt_quals"], CASE),
              ("controls", [q for n, v in site["libraries"].items() if n != case for q in v["alt_quals"]], CTRL)]
    for i, (_lab, qs, col) in enumerate(groups):
        y = 1 - i
        if qs:
            ax.scatter(qs, y + rng.uniform(-0.2, 0.2, len(qs)), s=5, color=col, alpha=0.8, linewidths=0)
            ax.plot([np.median(qs)] * 2, [y - 0.3, y + 0.3], color=INK, lw=1)
        ax.text(63, y, f"n = {len(qs)}", fontsize=6, color=MID, va="center", ha="right")
    ax.axvline(20, color=SOFT, lw=0.6, ls=":")
    ax.set_yticks([1, 0], [f"{case} (case)", "controls"] if show_y else [], fontsize=6)
    ax.set_ylim(-0.55, 1.55)
    ax.set_xlim(0, 64)
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_xlabel("Base quality of supporting reads")


def fdr_panel(ax, d, case):
    c = d["candidates"][case]
    pv = sorted(c["p"])
    decoys = [(sorted(p), b) for p, b in d["decoys"].values()]
    ts = np.logspace(np.log10(max(min(pv), 1e-30)), 0, 300)
    kept = np.array([bisect_right(pv, t) for t in ts])
    false = np.array([np.mean([bisect_right(p, t) * max(1.0, c["bases"] / b) for p, b in decoys]) for t in ts])
    ax.plot(ts, kept, color=CASE, lw=1.2)
    ax.plot(ts, false, color=CTRL, lw=1.2)
    ax.plot(ts, 0.05 * kept, color=CTRL, lw=0.7, ls=":")
    thr = c["threshold"]
    n_kept = bisect_right(pv, thr)
    ax.axvline(thr, color=INK, lw=0.7, ls="--")
    ax.scatter([thr], [n_kept], s=16, color=CASE, edgecolor="white", linewidth=0.5, zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1)
    top = max(kept) * 4
    ax.set_ylim(0, top)
    ax.set_xlim(ts[0], 1)          # a p-value threshold cannot exceed 1
    # direct labels in place of a legend
    i = np.searchsorted(ts, 1e-14)
    ax.text(ts[i], kept[i] * 1.6, "candidates kept", color=CASE, fontsize=6, ha="center", va="bottom")
    ax.text(1e-13, 0.04, "decoy-estimated false calls", color=CTRL, fontsize=6, ha="center", va="bottom")
    j = np.searchsorted(ts, 1e-24)
    ax.text(ts[j], 0.05 * kept[j] * 1.7, "5% of kept", color=CTRL, fontsize=6, ha="center", va="bottom")
    ax.text(thr * 0.6, top * 0.55, f"threshold\np ≤ {sci(thr)}\n{n_kept} calls", fontsize=6, color=INK,
            ha="right", va="top", fontweight="bold", linespacing=1.2)
    ax.set_xticks([1e-25, 1e-20, 1e-15, 1e-10, 1e-5, 1])
    ax.set_xlabel("p-value threshold")
    ax.set_ylabel("Candidates or false calls")
    ax.set_title(f"Empirical threshold, {LOCI[d['dataset']]} {case}", fontsize=6.6)


def af_panel(ax):
    bins = np.logspace(np.log10(0.04), np.log10(60), 34)
    total, nlib = 0, 0
    bottom = np.zeros(len(bins) - 1)
    colors = {"jh4": "#274b7e", "n114": "#6d9ad0", "su_lps": "#1f6f8b", "pim1": "#9fc0d8"}
    allaf = []
    for key, col in colors.items():
        d = load(key)
        nlib += len(d["cases"])
        af = [100 * c[3] / c[4] for n in d["cases"] for c in d["calls"][n] if c[4]]
        allaf += af
        h, _ = np.histogram(af, bins)
        ax.bar(bins[:-1], h, np.diff(bins), bottom=bottom, align="edge", color=col, label=LOCI[key],
               edgecolor="white", linewidth=0.5)
        bottom += h
        total += len(af)
    allaf = np.array(allaf)
    shares = {x: 100 * np.mean(allaf < x) for x in (1, 5, 20)}
    top = bottom.max() * 1.5
    for x in (1, 5):
        ax.axvline(x, ymax=0.74, color=SOFT, lw=0.5, ls=(0, (2, 2)))
    ax.text(44, top * 0.97, f"{shares[1]:.0f}% below 1% of reads\n{shares[5]:.0f}% below 5%", fontsize=6,
            ha="right", va="top", color=INK)
    ax.axvline(50, ymax=0.74, color=SOFT, lw=0.7, ls=":")
    ax.text(44, top * 0.42, "≥50%:\nexcluded as\ngermline", fontsize=6, color=MID, ha="right", va="center")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 1, 5, 10, 50], ["0.1", "1", "5", "10", "50"])
    ax.minorticks_off()
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, top)
    ax.set_xlabel("Reads carrying the call (%)")
    ax.set_ylabel("SHMtrace calls")
    ax.set_title(f"Allele fractions of {total:,} calls, {nlib} libraries", fontsize=6.6)
    ax.legend(fontsize=6, loc="upper left", handlelength=1, borderaxespad=0.2, ncol=2, columnspacing=0.8)
    return shares, total


def main():
    setup()
    d = load(SITE_SET)
    a = d["anatomy"]
    case = a["case"]
    fig = plt.figure(figsize=(FULL, 118 * MM))
    # left: the worked comparison (A, B and, under each, its qualities C); right: threshold and allele fractions
    outer = gridspec.GridSpec(1, 2, figure=fig, left=0.085, right=0.985, top=0.87, bottom=0.085, wspace=0.3,
                              width_ratios=[1.55, 1])
    left = outer[0].subgridspec(2, 2, wspace=0.08, hspace=0.62, height_ratios=[1.35, 0.62])
    right = outer[1].subgridspec(2, 1, hspace=0.62)
    for k, (key, label) in enumerate((("true", "AID hotspot"), ("error", "Recurrent background"))):
        site = a["sites"][key]
        ax = fig.add_subplot(left[0, k])
        site_panel(ax, site, case, show_y=(k == 0))
        site_title(ax, site, label)
        if k == 0:
            letter(ax, "A", dx=-0.36, dy=1.13)
        else:
            letter(ax, "B", dx=-0.05, dy=1.13)
        axq = fig.add_subplot(left[1, k])
        quality_strip(axq, site, case, show_y=(k == 0))
        if k == 0:
            letter(axq, "C", dx=-0.36, dy=1.12)
            axq.set_title(f"{case} supporting reads have similar quality at both sites; the controls differ",
                          fontsize=6.6, loc="left", pad=6)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (CASE, CASEF, CTRL, CTRLF)]
    fig.legend(handles, [f"{case} (AID-active), + strand", "− strand", "AID-null controls, + strand",
                         "− strand"], loc="upper left", ncol=4, fontsize=6, bbox_to_anchor=(0.02, 1.0),
               handlelength=1.2, columnspacing=1.4)
    ax = fig.add_subplot(right[0])
    fdr_panel(ax, d, case)
    letter(ax, "D", dx=-0.3)
    ax = fig.add_subplot(right[1])
    shares, total = af_panel(ax)
    letter(ax, "E", dx=-0.3)
    save(fig, "Figure2")
    json.dump({"case": case, "dataset": SITE_SET, "af_share_below": shares, "af_calls": total},
              open(SCRATCH / "fig2_numbers.json", "w"), indent=1)


if __name__ == "__main__":
    main()
