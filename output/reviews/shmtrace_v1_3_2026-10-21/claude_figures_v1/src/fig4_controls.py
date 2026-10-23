"""Figure 4: on real libraries with matched AID-null controls, the controls set the background.

A-E  unique-variant density per library (case blue, AID-null control red) for each caller, five datasets
F    separation: geometric-mean case density over control density (controls without calls at one event)
G    calls made in the AID-null control libraries, summed per dataset
H    leave-one-control-out calibration: share of held-out control tests with p < 0.01 against the margin
(Claude revision: the read accounting that was panel I is drawn by read_fates() as Figure S1E.)
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy import stats
from style import SCRATCH  # noqa: E402
from style import (
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

FIGS = DATA.parents[1] / "figures"
SHOW = ["recommended", "geneious", "lofreq", "lofreq_nobaq", "lofreq_permissive"]
SHORT = {"recommended": "SHMtrace", "geneious": "Geneious", "lofreq": "LoFreq",
         "lofreq_nobaq": "LoFreq −BAQ", "lofreq_permissive": "LoFreq perm."}
SHORTSET = {"jh4": "J$_H$4", "n114": "N114 V$_H$", "su_lps": "Sμ LPS", "su_ch12": "Sμ CH12", "pim1": "Pim1",
            "ot1": "Off-target 1", "ot2": "Off-target 2"}
# the per-library panels (A-E) show the five loci of the first pass; the summary panels F and G take every
# dataset, including the two CH12 AID off-targets added later
PANEL_SETS = [("jh4", "J$_H$4 intron (ONT)"), ("n114", "N114 V$_H$ (ONT)"), ("su_lps", "Sμ, LPS blasts (ONT)"),
              ("su_ch12", "Sμ, CH12 (ONT)"), ("pim1", "Pim1, LPS blasts (ONT)")]
SETS = PANEL_SETS + [("ot1", "CH12 off-target 1 (ONT)"), ("ot2", "CH12 off-target 2 (ONT)")]


def load(key):
    f = DATA / f"benchmark_real_{key}.json"
    return json.load(open(f))["rows"] if f.exists() else None


def geo(xs):
    return 10 ** np.mean(np.log10(xs))


# OT1 and OT2 carry a third biological group beside wild type, the MLLT1/MLLT3 double knockout, which is
# expected to sit at the control level. Pooling it into the case group would average the perturbation together
# with wild type and hide the very effect the panel is meant to show, so these two datasets get one row per
# AID-competent group. Every other dataset has a single case group.
SPLIT = {"ot1": [("Off-target 1 Wt", "WT"), ("Off-target 1 dKO", "dKO")],
         "ot2": [("Off-target 2 Wt", "WT"), ("Off-target 2 dKO", "dKO")]}


def separation(rows, mode, prefix=None):
    """Geometric-mean case density over control density, for one AID-competent group (prefix) or for all."""
    rs = [r for r in rows if r["mode"] == mode]
    dens = lambda r: r["density_x1e5"] if r["calls"] else 1e5 / r["bases"]   # noqa: E731
    case = [dens(r) for r in rs
            if not r["control"] and (prefix is None or r["sample"].startswith(prefix))]
    ctrl = [dens(r) for r in rs if r["control"]]
    return geo(case) / geo(ctrl), sum(r["calls"] for r in rs if r["control"]), len(ctrl)


def stars(p):
    """Significance marker. ns is spelled out, because an unmarked column reads as untested."""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def case_control(rows, mode):
    """Case and control densities for one caller, a zero-call library entered at one event on its exposure."""
    rs = [r for r in rows if r["mode"] == mode]
    dens = lambda r: r["density_x1e5"] if r["calls"] else 1e5 / r["bases"]   # noqa: E731
    return ([dens(r) for r in rs if not r["control"]], [dens(r) for r in rs if r["control"]])


def panel_tests(rows):
    """Welch on log10 density, AID-competent against AID-null, per caller; Holm-corrected within the panel."""
    raw = {}
    for c in SHOW:
        case, ctrl = case_control(rows, c)
        if len(case) < 2 or len(ctrl) < 2:
            continue
        a, b = np.log10(case), np.log10(ctrl)
        if np.ptp(a) == 0 and np.ptp(b) == 0:       # both groups constant: Welch is undefined
            continue
        raw[c] = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    out, m = {}, len(raw)
    for i, (c, p) in enumerate(sorted(raw.items(), key=lambda kv: kv[1])):
        out[c] = min(1.0, max(p * (m - i), *(out.values() or [0])))     # Holm, kept monotone
    return out


GROUPS = {"jh4": "3 Wt vs 3 null mice", "n114": "3 vs 3 subclones", "su_lps": "3 vs 3 libraries",
          "su_ch12": "4 vs 3 libraries", "pim1": "3 vs 3 libraries"}
PANEL_TITLES = {"jh4": "J$_H$4 intron", "n114": "N114 V$_H$", "su_lps": "Sμ, LPS blasts", "su_ch12": "Sμ, CH12",
                "pim1": "$Pim1$, LPS blasts"}
METHOD_ROWS = {"recommended": "SHMtrace", "geneious": "Geneious rule", "lofreq": "LoFreq, defaults",
               "lofreq_nobaq": "LoFreq, BAQ off", "lofreq_permissive": "LoFreq, permissive"}


def density_panel(ax, rows, title, sub, show_y):
    """One locus: a row per caller, AID-competent libraries above the row's centre and AID-null controls below.
    A library with no call has no density on a log axis, so it is drawn in a separate shaded '0' strip left of an
    axis break, never at a small positive value."""
    pos = [r["density_x1e5"] for r in rows if r["calls"]]
    lo_d, hi_d = min(pos), max(pos)
    xmin = lo_d / 4                      # left edge of the measured range
    zero_x = xmin / 7                    # centre of the zero strip
    strip = (zero_x / 2.2, zero_x * 2.2)
    rng = np.random.default_rng(1)
    ys = np.arange(len(SHOW))[::-1]
    pv = panel_tests(rows)
    for y, c in zip(ys, SHOW, strict=True):
        for r in rows:
            if r["mode"] != c:
                continue
            ctl = r["control"]
            yy = y + (-0.17 if ctl else 0.17) + rng.uniform(-0.07, 0.07)
            x = r["density_x1e5"] if r["calls"] else zero_x * 10 ** rng.uniform(-0.12, 0.12)
            ax.scatter(x, yy, s=9, marker="o", linewidth=0.5, facecolor=(CTRLF if ctl else CASEF),
                       edgecolor=(CTRL if ctl else CASE), zorder=3)
        if c in pv:
            ax.text(1.0, y, stars(pv[c]), transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=6, color=MID)
    for y in ys[:-1]:
        ax.axhline(y - 0.5, color=SOFT, lw=0.5, alpha=0.35)
    ax.set_xscale("log")
    xmax = hi_d * 9
    ax.set_xlim(strip[0], xmax)
    ax.axvspan(strip[0], strip[1], color=SOFT, alpha=0.15, lw=0)
    # break marks on the axis between the zero strip and the measured range
    for dx in (0.93, 1.07):
        ax.plot([strip[1] * dx * 1.08] * 2, [-0.62, -0.38], color=INK, lw=0.6, clip_on=False, zorder=5)
    decades = [10.0 ** k for k in range(int(np.ceil(np.log10(lo_d / 1.5))), int(np.ceil(np.log10(xmax))))]
    if len(decades) > 4:
        decades = decades[::-1][::2][::-1]
    ax.set_xticks([zero_x] + decades, ["0"] + [f"$10^{{{int(np.log10(d))}}}$" for d in decades])
    ax.minorticks_off()
    ax.set_ylim(-0.5, len(SHOW) - 0.5)
    ax.set_yticks(ys, [METHOD_ROWS[c] for c in SHOW] if show_y else [])
    ax.tick_params(axis="y", length=0)
    ax.set_title(f"{title}\n{sub}", fontsize=6.6, linespacing=1.25)


def heat_panel(ax, mat, row_labels, cmap, log, title, cbar_label, cax, fmt):
    show = np.log10(np.clip(mat, 1, None)) if log else mat
    vmax = np.nanmax(show)
    im = ax.imshow(show, cmap=cmap, aspect="auto", vmin=0, vmax=vmax * 1.1)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                continue
            dark = show[i, j] > vmax * 0.6
            ax.text(j, i, fmt(v), ha="center", va="center", fontsize=6, color="white" if dark else INK)
    ax.set_xticks(range(len(SHOW)), [SHORT[c] for c in SHOW], rotation=45, ha="right")
    ax.set_yticks(range(len(row_labels)), row_labels)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=6.6)
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    if log:
        ticks = [t for t in (1, 3, 10, 30, 100) if np.log10(t) <= vmax * 1.1]
        cb.set_ticks([np.log10(t) for t in ticks], labels=[f"{t}" for t in ticks])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    cb.set_label(cbar_label, fontsize=6, labelpad=2)


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 150 * MM))
    outer = gridspec.GridSpec(2, 1, figure=fig, left=0.115, right=0.985, top=0.885, bottom=0.075,
                              hspace=0.42, height_ratios=[1, 1.08])
    row1 = outer[0].subgridspec(1, 5, wspace=0.1)
    row2 = outer[1].subgridspec(2, 3, wspace=0.75, width_ratios=[1, 1, 1.25], height_ratios=[1, 0.045],
                                hspace=0.95)
    avail = [(k, lbl, pre) for k, _t in SETS if load(k)
             for lbl, pre in SPLIT.get(k, [(SHORTSET[k], None)])]
    panels = [(k, t) for k, t in PANEL_SETS if load(k)]
    tests = {}
    for i, (key, _title) in enumerate(panels):
        ax = fig.add_subplot(row1[0, i])
        rows = load(key)
        density_panel(ax, rows, PANEL_TITLES[key], GROUPS[key], show_y=(i == 0))
        tests[key] = {c: round(v, 5) for c, v in panel_tests(rows).items()}
        letter(ax, "ABCDE"[i], dx=-0.08, dy=1.2)
    fig.text(0.55, 0.52, "Unique variants per 10$^5$ analyzed bases (log scale; each panel has its own range)",
             ha="center", va="top", fontsize=6.5)
    h = [plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc=CASEF, mec=CASE, label="AID-competent library"),
         plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc=CTRLF, mec=CTRL, label="AID-null control library"),
         plt.Rectangle((0, 0), 1, 1, color=SOFT, alpha=0.3, label="0 calls (left of axis break)")]
    fig.legend(handles=h, loc="upper center", ncol=3, bbox_to_anchor=(0.55, 0.995), fontsize=6)
    # F: separation fold, one row per case group
    ax = fig.add_subplot(row2[0, 0])
    mat = np.array([[separation(load(k), c, pre)[0] for c in SHOW] for k, _l, pre in avail])
    labs = [lbl.replace("Off-target 1", "OT1").replace("Off-target 2", "OT2") for _k, lbl, _p in avail]
    heat_panel(ax, mat, labs, "Blues", True, "Case / control separation", "fold (geometric means)",
               fig.add_subplot(row2[1, 0]), lambda v: f"{v:.0f}" if v >= 10 else f"{v:.1f}")
    letter(ax, "F", dx=-0.62, dy=1.04)
    # G: calls in the AID-null controls, one row per dataset
    ax = fig.add_subplot(row2[0, 1])
    keys = [k for k, _t in SETS if load(k)]
    mat = np.array([[separation(load(k), c)[1] for c in SHOW] for k in keys], dtype=float)
    heat_panel(ax, mat, [SHORTSET[k].replace("Off-target ", "OT") for k in keys], "Reds", False,
               "Calls in AID-null controls", "calls, summed over 3 controls", fig.add_subplot(row2[1, 1]),
               lambda v: f"{v:.0f}")
    letter(ax, "G", dx=-0.62, dy=1.04)
    # H: leave-one-control-out trace
    ax = fig.add_subplot(row2[:, 2])
    runs = [("jh4_shmtrace", "J$_H$4, 529-bp window", "#2c5d8a", "-"),
            ("jh4_shmtrace_full", "J$_H$4, 1,141 bp", "#2c5d8a", "--"),
            ("n114_shmtrace", "N114 V$_H$", "#6a4c93", "-")]
    for run, lab, col, ls in runs:
        cal = json.load(open(FIGS / run / "project.json"))["calibration"]
        tr = cal["qscale_trace"]
        xs = sorted(tr, key=float)
        n = tr[xs[0]][1]
        ax.plot([float(x) for x in xs], [100 * tr[x][0] / tr[x][1] for x in xs], marker="o", ms=2.5, color=col,
                ls=ls, lw=1, label=f"{lab} ({n:,} tests)")
        sel = cal["qscale"]
        ax.scatter([sel], [100 * tr[f"{sel:.1f}"][0] / tr[f"{sel:.1f}"][1]], s=34, facecolor="none",
                   edgecolor=col, linewidth=0.9, zorder=4)
    ax.axhline(1, color=CTRL, lw=0.7, ls=":")
    ax.text(1.0, 0.25, "1%: expected if calibrated", fontsize=6, color=CTRL, va="bottom")
    ax.scatter([], [], s=34, facecolor="none", edgecolor=MID, linewidth=0.9, label="margin selected")
    ax.set_xlabel("Margin on control error rates (×)")
    ax.set_ylabel("Observed held-out tests with p < 0.01 (%)")
    ax.set_xscale("log")
    ticks = sorted({float(x) for run, *_ in runs
                    for x in json.load(open(FIGS / run / "project.json"))["calibration"]["qscale_trace"]})
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(0, 10.5)
    ax.set_title("Held-out control calibration", fontsize=6.6)
    ax.legend(fontsize=6, loc="upper right", handlelength=1.8)
    letter(ax, "H", dx=-0.25, dy=1.02)
    json.dump(tests, open(SCRATCH / "fig4_tests.json", "w"), indent=1)
    print("Figure 4A-E Welch p (Holm within panel):",
          {k: {c: f"{v:.3g}" for c, v in t.items()} for k, t in tests.items()})
    save(fig, "Figure4")


def read_fates():
    """The read accounting that was Figure 4I, now Figure S1E (page 2 of Figure S1)."""
    setup()
    fig = plt.figure(figsize=(FULL, 70 * MM))
    ax = fig.add_axes([0.075, 0.27, 0.69, 0.62])
    parts = [("too_long", "over-length", "#5f6b75"), ("low_quality", "expected error > 0.01", "#9aa5b1"),
             ("not_target", "unmapped or paralog", "#d9952b"), ("chimeric", "chimeric", "#c0392b"),
             ("short_alignment", "alignment too short", "#e6b98a"), ("paralogue_alleles", "paralog alleles", "#6a4c93"),
             ("usable", "usable", "#1f6f8b")]
    labels = []
    for run, name in (("jh4_shmtrace", "J$_H$4"), ("n114_shmtrace", "N114 V$_H$")):
        acct = json.load(open(FIGS / run / "table_s1.json"))
        for lib, a in acct.items():
            labels.append(f"{name} {lib}")
            left = 0
            for key, lab, col in parts:
                v = 100 * a[key] / a["delivered"]
                ax.bar(len(labels) - 1, v, bottom=left, color=col, width=0.75, label=lab if len(labels) == 1 else None)
                left += v
            ax.text(len(labels) - 1, 101.5, f"{a['delivered']:,}", ha="center", va="bottom", fontsize=6, color=MID,
                    rotation=90)
    ax.set_xticks(range(len(labels)), labels, fontsize=6, rotation=40, ha="right")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylabel("Delivered reads (%)")
    ax.set_title("Fate of every delivered read (numbers above bars: reads delivered)", fontsize=6.6, pad=26)
    ax.legend([h for h in ax.get_legend_handles_labels()[0]][::-1], [lab for _k, lab, _c in parts][::-1],
              fontsize=6, loc="center left", bbox_to_anchor=(1.01, 0.5), handlelength=1)
    letter(ax, "E", dx=-0.075, dy=1.12)
    save(fig, "FigureS1_read_fates")


if __name__ == "__main__":
    main()
    read_fates()
