"""Supplementary figures whose corrections come from the Codex v2 revision, ported to this revision's style.

Ported from codex_figures_v2/sources/revised/custom_figures.py (matrix, matrices, s17, background), with the same
data sources and counting; the spectra are drawn in the style of Figures 5 and 7 rather than as From -> To matrices:

FigureS9_n114           N114 density per subclone on a linear axis, so a zero-call subclone sits at zero rather than
                        at a one-event value on a log axis, with the separation; the pMX-Aicda spectrum by subclone
FigureS17_spectra_full_length  the Wt spectrum in the shared window and the full interval, and their shares plotted
                        against each other, with the windows' nesting stated
FigureS11               A-B: controls on separate tracks sharing one scale per locus (the AID-active mean on its
                        own scale); C-D: the read maps of the original figure (figS11_background.raster_panel)
A spectrum shows counts under its shares, and the AID-null libraries are named, not plotted: a share of one or two
control calls is not interpretable.
"""
import json
from collections import Counter

import numpy as np
from matplotlib import pyplot as plt
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, mouse, save, setup

CAL = DATA.parents[1]
BASES = "ATCG"
# substitution types in the order of Figure 5: from G:C first, then from A:T; coloured by the base pair mutated
ORDER = ["C>T", "C>A", "C>G", "G>A", "G>T", "G>C", "A>G", "A>C", "A>T", "T>C", "T>A", "T>G"]
GCC, ATC = "#2c5d8a", "#d9952b"
GCL, ATL = "#a9bfd6", "#f0cf9c"          # the same two, lighter, for the nested window in Figure S17


def counts(vs):
    return Counter(v["ref"] + ">" + v["alt"] for v in vs)


def head(ax, title, sub):
    """Bold title with a grey subtitle beneath it, as in Figure 7."""
    ax.set_title(title, fontsize=7, pad=10)
    ax.text(0.0, 1.012, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=MID)


def type_axis(ax, rows):
    """Substitution types on x; under each, one row of counts per entry of rows ((label, Counter, colour)); then
    the from-G:C and from-A:T brackets."""
    ax.set_xticks(range(12), [t.replace(">", "→") for t in ORDER])
    ax.set_xlim(-0.6, 11.6)
    for k, (lab, c, col) in enumerate(rows):
        y = -0.115 - 0.075 * k
        for x, t in enumerate(ORDER):
            ax.text(x, y, str(c[t]), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6, color=col)
        if lab:
            ax.text(-0.75, y, lab, transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=6, color=col)
    yb = -0.2 - 0.075 * (len(rows) - 1)
    for x0, x1, lab, col in ((0, 5, "from G:C", GCC), (6, 11, "from A:T", ATC)):
        ax.annotate("", xy=(x0 - 0.3, yb), xytext=(x1 + 0.3, yb), xycoords=("data", "axes fraction"),
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.0))
        ax.text((x0 + x1) / 2, yb - 0.03, lab, transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6, color=col)


def spectrum(ax, libs, title, unit, who, one, null_note):
    """The AID-active libraries' spectrum: bars, pooled shares by type, with the count under each; points, each
    library (circle, square, triangle). The AID-null libraries are named in red, since a handful of calls is not a
    spectrum."""
    pooled = sum((c for _, c in libs), Counter())
    tot = sum(pooled.values())
    for x, t in enumerate(ORDER):
        ax.bar(x, 100 * pooled[t] / tot, width=0.66, color=GCC if t[0] in "GC" else ATC, alpha=0.85, lw=0)
    top = 0
    for (name, c), mk, dx in zip(libs, "os^", (-0.17, 0, 0.17), strict=True):
        n = sum(c.values())
        ys = [100 * c[t] / n for t in ORDER]
        top = max(top, *ys)
        ax.scatter([x + dx for x in range(12)], ys, s=9, marker=mk, facecolor="white", edgecolor=INK, lw=0.5,
                   zorder=3, label=f"{name} ({n} {unit})")
    ax.set_ylim(0, 5 * np.ceil(1.6 * top / 5))           # head room for the legend and the red note
    ax.set_ylabel(f"Share of {unit} (%)")
    type_axis(ax, [("", pooled, MID)])
    gc = sum(pooled[t] for t in ORDER if t[0] in "GC")
    head(ax, title, f"{tot} {unit} pooled over three {who}, {100 * gc / tot:.0f}% at G:C; {unit} under each type")
    ax.legend(title=f"each {one}", title_fontsize=6, loc="upper right",
              fontsize=6, handletextpad=0.2)
    ax.text(0.62, 0.96, null_note, transform=ax.transAxes, ha="center", va="top", fontsize=6, color=CTRL)


def null_note(group, libs, unit):
    """The AID-null libraries in words: "<group>: 0, 2 and 0 calls (KO1–KO3),\ntoo few for a spectrum"."""
    ns = [sum(c.values()) for _, c in libs]
    if sum(ns) == 1:
        (t,) = [t for _, c in libs for t in c]
        return f"{group}: 1 {unit[:-1]} ({t.replace('>', '→')}),\ntoo few for a spectrum"
    nums = ", ".join(str(n) for n in ns[:-1]) + f" and {ns[-1]}"
    return f"{group}: {nums} {unit} ({libs[0][0]}–{libs[-1][0]}),\ntoo few for a spectrum"


def project(key):
    return json.load(open(CAL / f"figures/{key}_shmtrace/project.json"))


def s9():
    """N114: density per subclone (A) and the pMX-Aicda spectrum (B), in the style of Figures 5 and 7. (A second page,
    the J_H4 spectrum by mouse, was dropped, as was the Wt Sanger spectrum of Figure S13: Figures 5B-C and S17A
    already show them.)"""
    p = project("n114")
    fig = plt.figure(figsize=(FULL, 74 * MM))
    ax = fig.add_axes([0.08, 0.24, 0.17, 0.58])
    logs = []
    for x, (g, f, e) in enumerate((("pMX-Aicda", CASEF, CASE), ("pMX-Empty", CTRLF, CTRL))):
        ss = [s for s in p["samples"] if s["group"] == g]
        ax.scatter(x + np.linspace(-0.12, 0.12, len(ss)), [s["unique_density_x1e5"] for s in ss], s=14, facecolor=f,
                   edgecolor=e, lw=0.6, zorder=3)
        logs.append(np.log10([max(s["unique_density_x1e5"], s["one_event_density_x1e5"]) for s in ss]))
        gm = 10 ** logs[-1].mean()
        ax.plot([x - 0.25, x + 0.25], [gm, gm], color=e, lw=1.4)
    ax.set_xticks([0, 1], ["pMX-$\\it{Aicda}$", "empty\nvector"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-2, 72)
    # the separation the paper reports (Figure 4, Table S1): Welch on log density, zero-call subclones at one event
    from scipy import stats
    fold = 10 ** (logs[0].mean() - logs[1].mean())
    pv = stats.ttest_ind(*logs, equal_var=False).pvalue
    assert (round(fold), f"{pv:.2g}") == (85, "0.0071"), (fold, pv)
    ax.text(0.5, 70, f"{fold:.0f}×\np = {pv:.2g}", ha="center", va="top", fontsize=6, color=MID)
    ax.set_ylabel("Unique variants per 10$^5$ bases")
    head(ax, "Density per subclone", "linear axis; bars, geometric means")
    letter(ax, "A", dx=-0.45, dy=1.08)
    ax = fig.add_axes([0.36, 0.24, 0.62, 0.58])
    lib = lambda n: (n, counts([v for v in p["variants"] if v["sample"] == n]))   # noqa: E731
    nulls = [lib(n) for n in ("C3", "C4", "D1")]
    spectrum(ax, [lib(n) for n in ("C1", "D3", "D5")], "pMX-$\\it{Aicda}$ substitution spectrum", "calls",
             "subclones", "subclone", null_note("empty vector", nulls, "calls"))
    letter(ax, "B", dx=-0.09, dy=1.08)
    save(fig, "FigureS9_n114")
    plt.close(fig)


def s17():
    """The Wt spectrum in the shared window and in the full interval: paired bars by type (lighter, the window) with
    both counts under each type (A), and the two sets of shares against each other (B)."""
    c1 = counts([v for v in project("jh4")["variants"] if v["sample"].startswith("WT")])
    q = json.load(open(CAL / "figures/jh4_shmtrace_full/project.json"))
    c2 = counts([v for v in q["variants"] if v["sample"].startswith("WT")])
    n1, n2 = sum(c1.values()), sum(c2.values())
    fig = plt.figure(figsize=(FULL, 80 * MM))
    ax = fig.add_axes([0.1, 0.3, 0.55, 0.52])
    for x, t in enumerate(ORDER):
        g = t[0] in "GC"
        ax.bar(x - 0.18, 100 * c1[t] / n1, width=0.34, color=GCL if g else ATL, lw=0)
        ax.bar(x + 0.18, 100 * c2[t] / n2, width=0.34, color=GCC if g else ATC, alpha=0.85, lw=0)
    ax.set_ylabel("Share of Wt calls (%)")
    ax.set_ylim(0, 25)
    type_axis(ax, [("window", c1, SOFT), ("interval", c2, MID)])
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GCL, label=f"shared 529-bp window ({n1} calls)"),
                       Patch(color=GCC, alpha=0.85, label=f"full 1,141-bp interval ({n2} calls)")],
              loc="upper right", fontsize=6, handlelength=1, handletextpad=0.4)
    head(ax, "J$_H$4 Wt spectrum in nested windows", "lighter bar, the window; darker, the interval that contains it; "
         "calls under each type")
    letter(ax, "A", dx=-0.13, dy=1.08)
    ax = fig.add_axes([0.76, 0.3, 0.21, 0.52])
    xs = [100 * c1[t] / n1 for t in ORDER]
    ys = [100 * c2[t] / n2 for t in ORDER]
    ax.plot([0, 22], [0, 22], ls="--", lw=0.6, color=SOFT)
    ax.scatter(xs, ys, s=14, c=[GCC if t[0] in "GC" else ATC for t in ORDER], edgecolor="white", lw=0.5, zorder=3)
    ax.set(xlim=(0, 22), ylim=(0, 22))
    ax.set_xlabel("529-bp window (% of calls)")
    ax.set_ylabel("1,141-bp interval (% of calls)")
    r = float(np.corrcoef(xs, ys)[0, 1])
    ax.text(0.04, 0.96, f"r = {r:.3f}\nr² = {r * r:.3f}", transform=ax.transAxes, fontsize=6, va="top")
    ax.text(0.96, 0.04, "nested windows:\ndescriptive, not\nindependent", transform=ax.transAxes, fontsize=6,
            color=MID, va="bottom", ha="right")
    head(ax, "Shares against each other", "dashed line, identity")
    letter(ax, "B", dx=-0.42, dy=1.08)
    save(fig, "FigureS17_spectra_full_length")
    plt.close(fig)


def s11():
    from figS11_background import raster_panel, rates
    fig = plt.figure(figsize=(FULL, 170 * MM))
    for j, (key, title) in enumerate((("jh4", "J$_H$4"), ("su_lps", "Sμ"))):
        d = json.load(open(DATA / f"intuition_{key}.json"))
        lo, hi = d["window"]
        xs = np.arange(lo, hi + 1)
        arrays = [100 * np.nanmean([rates(d, n, lo, hi) for n in d["cases"]], axis=0)]
        arrays += [100 * rates(d, n, lo, hi) for n in d["controls"]]
        names = ["AID-active\nmean"] + list(d["controls"])
        ctrlmax = max(np.nanmax(v) for v in arrays[1:]) * 1.08
        for k, (vals, n) in enumerate(zip(arrays, names, strict=True)):
            ax = fig.add_axes([0.14 + j * 0.49, 0.85 - k * 0.085, 0.33, 0.068])
            ax.plot(xs, vals, color=CASE if k == 0 else CTRL, lw=0.6)
            top = np.nanmax(vals) * 1.08 if k == 0 else ctrlmax
            ax.set_xlim(lo, hi)
            ax.set_ylim(0, top)
            ax.set_yticks([0, round(top * 0.8, 1)])
            ax.set_ylabel(n, rotation=0, ha="right", va="center", labelpad=6)
            if k < len(arrays) - 1:
                ax.set_xticks([])
            else:
                ax.set_xlabel(f"{title} position (bp)", labelpad=2)
            if k == 0:
                ax.set_title(f"{title}: mismatches at Q ≥ 20 (% of bases)", fontsize=6.6)
                letter(ax, "AB"[j], dx=-0.36, dy=1.15)
            if k == 1:
                ax.text(1.0, 1.04, "controls: one shared scale", transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=6, color=MID)
    d = json.load(open(DATA / "intuition_jh4.json"))
    for j, (n, label) in enumerate((("WT3", "wild type"), ("KO1", "$\\it{Aicda}^{{-}/{-}}$"))):
        ax = fig.add_axes([0.1 + j * 0.49, 0.07, 0.37, 0.36])
        raster_panel(ax, d, n, label)
        for t in list(ax.texts):
            if t.get_text() in ("C", "D"):
                t.remove()
        letter(ax, "CD"[j], dx=-0.2, dy=1.04)
    save(fig, "FigureS11")
    plt.close(fig)


def main():
    setup()
    s9()
    s17()
    s11()


if __name__ == "__main__":
    main()
