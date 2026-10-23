"""Supplementary figures whose corrections come from the Codex v2 revision, ported to this revision's style.

Ported from codex_figures_v2/sources/revised/custom_figures.py (matrix, matrices, s17, background), with the same
data sources and counting, laid out more compactly:

FigureS9_n114           N114 density per subclone on a linear axis, so a zero-call subclone sits at zero rather than
                        at a one-event value on a log axis; pooled From -> To counts for each vector
FigureS9_spectra_jh4    one count matrix per J_H4 library, one shared colour scale; zero-call spectra marked undefined
FigureS9_spectra_n114   the same for N114 V_H
FigureS13_sanger_matrices  per-mouse Sanger Q>=20 counts, one shared scale; tiny control totals marked descriptive
FigureS17_spectra_full_length  pooled Wt counts in the shared window and the full interval, and their shares plotted
                        against each other, with the windows' nesting stated
FigureS11               A-B: controls on separate tracks sharing one scale per locus (the AID-active mean on its
                        own scale); C-D: the read maps of the original figure (figS11_background.raster_panel)
Matrices give counts, not percentages: a percentage spectrum of one or two control calls is not interpretable.
"""
import json
from collections import Counter

import numpy as np
from matplotlib import pyplot as plt
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, letter, mouse, save, setup

CAL = DATA.parents[1]
BASES = "ATCG"
MUT = "#6a4c93"


def counts(vs):
    return Counter(v["ref"] + ">" + v["alt"] for v in vs)


def matrix(ax, cnt, title, vmax):
    total = sum(cnt.values())
    arr = np.full((4, 4), np.nan)
    for i, x in enumerate(BASES):
        for j, y in enumerate(BASES):
            if i != j:
                arr[i, j] = cnt.get(x + ">" + y, 0)
    im = ax.imshow(np.ma.masked_invalid(arr), cmap="Blues", vmin=0, vmax=max(1, vmax))
    for i in range(4):
        for j in range(4):
            if i != j:
                v = arr[i, j]
                ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6,
                        color="white" if v and v > vmax * 0.55 else INK)
    ax.set_xticks(range(4), list(BASES))
    ax.set_yticks(range(4), list(BASES))
    ax.tick_params(length=0, pad=1.5)
    ax.set_xlabel("To", labelpad=1)
    ax.set_ylabel("From", labelpad=1)
    ax.set_title(f"{title}: {total:,} call{'' if total == 1 else 's'}", fontsize=6.6, pad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)
    note = "no calls: spectrum undefined" if total == 0 else "descriptive only" if total < 5 else ""
    if note:
        ax.text(0.5, -0.3, note, transform=ax.transAxes, fontsize=6, ha="center", color=MID)
    return im


def colourbar(fig, im, rect, label):
    """The colour-scale legend of the original matrices: one bar per page, since every matrix on it shares the
    scale."""
    cax = fig.add_axes(rect)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    cb.set_label(label, fontsize=6, labelpad=2)


def grid(names, cnts, name, groups, h=110, titles=None, label="unique variants (count)"):
    fig = plt.figure(figsize=(FULL, h * MM))
    vmax = max(max(c.values(), default=0) for c in cnts)
    for i, (n, c) in enumerate(zip(names, cnts, strict=True)):
        ax = fig.add_axes([0.07 + (i % 3) * 0.33, 0.63 - (i // 3) * 0.42, 0.22, 0.26])
        im = matrix(ax, c, (titles or names)[i], vmax)
        if i % 3 == 0:
            letter(ax, "AB"[i // 3], dx=-0.32, dy=1.16)
            ax.text(-0.12, 1.16, groups[i // 3], transform=ax.transAxes, fontsize=7.5, fontweight="bold",
                    va="bottom", ha="left")
    colourbar(fig, im, [0.3, 0.085, 0.4, 0.018], label + "; one scale for every matrix on this page")
    save(fig, name)
    plt.close(fig)


def project(key):
    return json.load(open(CAL / f"figures/{key}_shmtrace/project.json"))


def s9():
    for key, names in (("jh4", ["WT2", "WT3", "WT4", "KO1", "KO2", "KO3"]),
                       ("n114", ["C1", "D3", "D5", "C3", "C4", "D1"])):
        p = project(key)
        groups = ("Wt", "$\\it{Aicda}^{{-}/{-}}$") if key == "jh4" else ("pMX-$\\it{Aicda}$", "empty vector")
        grid(names, [counts([v for v in p["variants"] if v["sample"] == n]) for n in names],
             f"FigureS9_spectra_{key}", groups, titles=[mouse(n) for n in names])
    # page 1 in the style of Figures 5 and 7: points with their group's geometric-mean bar, the separation in grey
    # above, and the spectrum as shares of calls by substitution type, coloured by the base pair mutated
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
    ax.set_title("Density per subclone", fontsize=7, pad=10)
    ax.text(0.0, 1.012, "linear axis; bars, geometric means", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=6, color=MID)
    letter(ax, "A", dx=-0.45, dy=1.08)
    # B: the pMX-Aicda spectrum, pooled shares as bars (count above), each subclone as a point
    ax = fig.add_axes([0.36, 0.24, 0.62, 0.58])
    order = ["C>T", "C>A", "C>G", "G>A", "G>T", "G>C", "A>G", "A>C", "A>T", "T>C", "T>A", "T>G"]
    GCC, ATC = "#2c5d8a", "#d9952b"
    case = [v for v in p["variants"] if v["group"] == "pMX-Aicda"]
    pooled, tot = counts(case), len(case)
    for x, t in enumerate(order):
        h = 100 * pooled[t] / tot
        ax.bar(x, h, width=0.66, color=GCC if t[0] in "GC" else ATC, alpha=0.85, lw=0)
        ax.text(x, -0.115, str(pooled[t]), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6,
                color=MID)
    for (sub, mk), dx in zip((("C1", "o"), ("D3", "s"), ("D5", "^")), (-0.17, 0, 0.17), strict=True):
        c = counts([v for v in case if v["sample"] == sub])
        n = sum(c.values())
        ax.scatter([x + dx for x in range(12)], [100 * c[t] / n for t in order], s=9, marker=mk, facecolor="white",
                   edgecolor=INK, lw=0.5, zorder=3, label=f"{sub} ({n} calls)")
    ax.set_xticks(range(12), [t.replace(">", "→") for t in order])
    ax.set_xlim(-0.6, 11.6)
    ax.set_ylim(0, 45)
    ax.set_ylabel("Share of calls (%)")
    gc = sum(pooled[t] for t in order if t[0] in "GC")
    ax.set_title("pMX-$\\it{Aicda}$ substitution spectrum", fontsize=7, pad=10)
    ax.text(0.0, 1.012, f"{tot} calls pooled over three subclones, {100 * gc / tot:.0f}% at G:C; calls under each type", transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=MID)
    for x0, x1, lab, col in ((0, 5, "from G:C", GCC), (6, 11, "from A:T", ATC)):
        ax.annotate("", xy=(x0 - 0.3, -0.2), xytext=(x1 + 0.3, -0.2), xycoords=("data", "axes fraction"),
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.0))
        ax.text((x0 + x1) / 2, -0.23, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6,
                color=col)
    from matplotlib.patches import Patch
    hs = ax.get_legend_handles_labels()[0]
    ax.legend(handles=hs, title="each subclone", title_fontsize=6, loc="upper right", fontsize=6, handletextpad=0.2)
    empty = [v for v in p["variants"] if v["group"] == "pMX-Empty"]
    ax.text(0.62, 0.96, f"empty vector: {len(empty)} call "
            f"({', '.join(v['ref'] + '→' + v['alt'] for v in empty)}),\ntoo few for a spectrum",
            transform=ax.transAxes, ha="center", va="top", fontsize=6, color=CTRL)
    letter(ax, "B", dx=-0.09, dy=1.08)
    save(fig, "FigureS9_n114")
    plt.close(fig)


def s13():
    d = json.load(open(CAL / "manuscript_corrected/sanger_calls/corrected_sanger.json"))
    names = [n for n in d if n != "source"]
    # each variant is [position, ref, alt, clones carrying it]: one count per distinct substitution
    cnts = [Counter(v[1] + ">" + v[2] for v in d[n]["variants"]) for n in names]
    grid(names, cnts, "FigureS13_sanger_matrices", ("Wt", "$\\it{Aicda}^{{-}/{-}}$"),
         titles=[mouse(n) for n in names],
         label="distinct Q≥20 substitutions (count)")


def s17():
    c1 = counts([v for v in project("jh4")["variants"] if v["sample"].startswith("WT")])
    q = json.load(open(CAL / "figures/jh4_shmtrace_full/project.json"))
    c2 = counts([v for v in q["variants"] if v["sample"].startswith("WT")])
    fig = plt.figure(figsize=(FULL, 74 * MM))
    vmax = max(max(c1.values()), max(c2.values()))
    for i, (c, t) in enumerate(((c1, "Shared 529-bp window"), (c2, "Full 1,141-bp interval"))):
        ax = fig.add_axes([0.07 + i * 0.3, 0.33, 0.2, 0.52])
        im = matrix(ax, c, t, vmax)
        letter(ax, "AB"[i], dx=-0.32, dy=1.06)
    colourbar(fig, im, [0.1, 0.17, 0.45, 0.025], "Wt unique variants (count); one scale for A and B")
    ax = fig.add_axes([0.75, 0.2, 0.22, 0.65])
    ts = [x + ">" + y for x in BASES for y in BASES if x != y]
    xs = [100 * c1[t] / sum(c1.values()) for t in ts]
    ys = [100 * c2[t] / sum(c2.values()) for t in ts]
    ax.plot([0, 22], [0, 22], ls="--", lw=0.6, color="#8a96a3")
    ax.scatter(xs, ys, s=14, color=MUT, zorder=3)
    ax.set(xlim=(0, 22), ylim=(0, 22))
    ax.set_xlabel("529-bp window (% of calls)")
    ax.set_ylabel("1,141-bp interval (% of calls)")
    ax.set_title("Shares, nested windows", fontsize=6.6)
    r = float(np.corrcoef(xs, ys)[0, 1])
    ax.text(0.04, 0.96, f"Pearson r = {r:.3f}\nr² = {r * r:.3f}", transform=ax.transAxes, fontsize=6, va="top")
    ax.text(0.96, 0.04, "nested windows:\ndescriptive, not\nindependent", transform=ax.transAxes, fontsize=6,
            color=MID, va="bottom", ha="right")
    letter(ax, "C", dx=-0.3, dy=1.06)
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
    s13()
    s17()
    s11()


if __name__ == "__main__":
    main()
