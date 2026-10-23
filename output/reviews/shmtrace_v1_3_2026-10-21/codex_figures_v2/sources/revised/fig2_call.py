"""Figure 2: how a call is made, and where real SHM sits.

A    a true AID hotspot in one AID-active library: share of reads carrying it, by strand, against each control
B    a recurrent background site in the same library, which a read-level test calls
C    the base qualities of the alt reads at both sites: quality alone cannot tell them apart
D    the control false-discovery threshold for the same library: candidates kept against the false calls the
     held-out controls predict
E    the allele fractions of SHMtrace's calls, all AID-active libraries of four loci, against the range in which
     the somatic callers were tested on planted mutations (Figure 3C) and their recall there
Data: out/intuition_<dataset>.json (intuition_data.py).
"""
import json
from bisect import bisect_right

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
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


def site_panel(ax, site, case, title):
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
    ax.set_yticks(ys, [f"{n}\n(case)" if n == case else n for n in libs], fontsize=6)
    ax.set_xlim(0, 2.3)
    ax.set_xticks([0, 0.5, 1, 1.5])
    ax.set_xlabel("Reads carrying the change, any quality (%)")
    ctx = site["context"]
    ax.set_title(f"{title}: {site['pos']} {site['ref']}>{site['alt']} "
                 f"({ctx[:2]}$\\underline{{{ctx[2]}}}${ctx[3:]})", fontsize=6.6)
    verdict = "called" if site["called"] else "not called"
    ax.text(0.99, 0.42, f"SHMtrace\np = {sci(site['p'])}\n{verdict}", transform=ax.transAxes, ha="right",
            va="center", fontsize=6, color=INK, linespacing=1.3)


def fdr_panel(ax, d, case):
    c = d["candidates"][case]
    pv = sorted(c["p"])
    decoys = [(sorted(p), b) for p, b in d["decoys"].values()]
    ts = np.logspace(np.log10(max(min(pv), 1e-30)), 0, 300)
    kept = [bisect_right(pv, t) for t in ts]
    false = [np.mean([bisect_right(p, t) * max(1.0, c["bases"] / b) for p, b in decoys]) for t in ts]
    ax.plot(ts, kept, color=CASE, lw=1.1, label=f"{case} candidates with p ≤ t")
    ax.plot(ts, false, color=CTRL, lw=1.1, label="decoy-estimated false calls")
    ax.plot(ts, [0.05 * k for k in kept], color=CTRL, lw=0.7, ls=":", label="5% of the candidates kept")
    thr = c["threshold"]
    n_kept = bisect_right(pv, thr)
    ax.axvline(thr, color=INK, lw=0.6, ls="--", label=f"threshold p = {sci(thr)}: {n_kept} calls")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, max(kept) * 2.5)
    ax.set_xlim(ts[0], 1)          # a p-value threshold cannot exceed 1; do not let the log view run to 10
    ax.set_xlabel("p-value threshold t")
    ax.set_ylabel("Count")
    ax.set_title(f"Empirical threshold, {LOCI[d['dataset']]} {case}", fontsize=6.6)
    ax.legend(fontsize=6.5, loc="upper left", handlelength=1.3, borderaxespad=0.2)


def af_panel(ax):
    bins = np.logspace(np.log10(0.04), np.log10(60), 34)
    shares, total = {}, 0
    bottom = np.zeros(len(bins) - 1)
    colors = {"jh4": "#274b7e", "n114": "#6d9ad0", "su_lps": "#1f6f8b", "pim1": "#9fc0d8"}
    allaf = []
    for key, col in colors.items():
        d = load(key)
        af = [100 * c[3] / c[4] for n in d["cases"] for c in d["calls"][n] if c[4]]
        allaf += af
        h, _ = np.histogram(af, bins)
        ax.bar(bins[:-1], h, np.diff(bins), bottom=bottom, align="edge", color=col, label=LOCI[key],
               edgecolor="white", linewidth=0.5)
        bottom += h
        total += len(af)
    allaf = np.array(allaf)
    shares = {x: 100 * np.mean(allaf < x) for x in (1, 5, 20)}
    top = bottom.max() * 1.45
    ax.axvline(50, color=SOFT, lw=0.7, ls=":")
    ax.text(46, top * 0.38, "heterozygous\ngermline", fontsize=6, color=MID, ha="right", va="center")
    ax.text(0.05, top * 0.97, f"{shares[1]:.0f}% of {total:,} calls in < 1% of reads\n"
            f"{shares[5]:.0f}% in < 5%", fontsize=6, va="top", color=INK)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 1, 10, 50], ["0.1", "1", "10", "50"])
    ax.minorticks_off()
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, top)
    ax.set_xlabel("Reads carrying the call (%)")
    ax.set_ylabel("SHMtrace calls")
    ax.set_title("Allele fractions of SHMtrace calls", fontsize=6.6)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0.0, 0.84), handlelength=1, borderaxespad=0.2,
              ncol=2, columnspacing=0.8)
    return shares, total


def quality_panel(ax, a, case):
    groups = []
    for label in ("true", "error"):
        s = a["sites"][label]
        groups.append((f"{s['pos']}\n{case}", s["libraries"][case]["alt_quals"], CASE))
        ctl = [q for n, v in s["libraries"].items() if n != case for q in v["alt_quals"]]
        groups.append((f"{s['pos']}\ncontrols", ctl, CTRL))
    rng = np.random.default_rng(0)
    for i, (_lab, qs, col) in enumerate(groups):
        if qs:
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(qs)), qs, s=4, color=col, alpha=0.75, linewidths=0)
            ax.plot([i - 0.28, i + 0.28], [np.median(qs)] * 2, color=INK, lw=1)
    ax.axhline(20, color=SOFT, lw=0.6, ls=":")
    ax.text(1, 20.8, "Q20 reference", fontsize=6, color=MID, ha="center", va="bottom")
    ax.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=6)
    ax.set_ylabel("Base quality of alt reads")
    ax.set_ylim(0, 52)
    ax.set_title("Quality cannot tell them apart", fontsize=6.6)


def main():
    setup()
    d = load(SITE_SET)
    a = d["anatomy"]
    case = a["case"]
    fig = plt.figure(figsize=(FULL, 112 * MM))
    outer = gridspec.GridSpec(2, 1, figure=fig, left=0.075, right=0.985, top=0.87, bottom=0.09, hspace=0.62)
    top = outer[0].subgridspec(1, 3, wspace=0.62, width_ratios=[1, 1, 0.8])
    bottom = outer[1].subgridspec(1, 2, wspace=0.28, width_ratios=[1, 1.15])
    ax = fig.add_subplot(top[0])
    site_panel(ax, a["sites"]["true"], case, "AID hotspot")
    letter(ax, "A", dx=-0.3)
    ax = fig.add_subplot(top[1])
    site_panel(ax, a["sites"]["error"], case, "Recurrent background")
    letter(ax, "B", dx=-0.3)
    ax = fig.add_subplot(top[2])
    quality_panel(ax, a, case)
    letter(ax, "C", dx=-0.32)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (CASE, CASEF, CTRL, CTRLF)]
    fig.legend(handles, [f"{case} (AID-active), + strand", "− strand", "AID-null controls, + strand",
                         "− strand"], loc="upper center", ncol=4, fontsize=6, bbox_to_anchor=(0.5, 1.0),
               handlelength=1.2, columnspacing=1.4)
    ax = fig.add_subplot(bottom[0])
    fdr_panel(ax, d, case)
    letter(ax, "D", dx=-0.16)
    ax = fig.add_subplot(bottom[1])
    shares, total = af_panel(ax)
    letter(ax, "E", dx=-0.14)
    save(fig, "Figure2")
    json.dump({"case": case, "dataset": SITE_SET, "af_share_below": shares, "af_calls": total},
              open('/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21/codex_version/scratch/fig2_numbers.json', "w"), indent=1)


if __name__ == "__main__":
    main()
