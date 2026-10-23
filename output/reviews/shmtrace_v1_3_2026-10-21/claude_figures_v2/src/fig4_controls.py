"""Figure 4: on real libraries with matched AID-null controls, the controls set the background.

A-E  unique-variant density per library (case blue, AID-null control red) for each caller, five datasets
F    separation: geometric-mean case density over control density (controls without calls at one event)
G    calls made in the AID-null control libraries, summed per dataset
H    leave-one-control-out calibration: share of held-out control tests with p < 0.01 against the margin
(Claude revision: the read accounting that was panel I is drawn by read_fates() as its own Figure S2.)
"""
import json
from pathlib import Path

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
    mouse,
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


GROUPS = {"jh4": "3 Wt vs 3 $\\it{Aicda}^{{-}/{-}}$ mice", "n114": "3 vs 3 subclones", "su_lps": "3 vs 3 libraries",
          "su_ch12": "4 vs 3 libraries", "pim1": "3 vs 3 libraries"}
PANEL_TITLES = {"jh4": "J$_H$4 intron", "n114": "N114 V$_H$", "su_lps": "Sμ, LPS blasts", "su_ch12": "Sμ, CH12",
                "pim1": "$Pim1$, LPS blasts"}
METHOD_ROWS = {"recommended": "SHMtrace", "geneious": "Geneious rule", "lofreq": "LoFreq, defaults",
               "lofreq_nobaq": "LoFreq, BAQ off", "lofreq_permissive": "LoFreq, permissive"}


def density_panel(ax, rows, title, sub, show_y):
    """One locus in the original orientation: a column pair per caller, AID-competent libraries left and AID-null
    controls right. A library with no call has no density on a log axis, so it is drawn in a separate shaded '0'
    strip below an axis break, never at a small positive value."""
    pos = [r["density_x1e5"] for r in rows if r["calls"]]
    lo_d, hi_d = min(pos), max(pos)
    zero_y = lo_d / 4.5                  # centre of the zero strip
    strip = (zero_y / 1.8, zero_y * 1.8)
    rng = np.random.default_rng(1)
    pv = panel_tests(rows)
    ymax = hi_d * 12
    for j, c in enumerate(SHOW):
        for r in rows:
            if r["mode"] != c:
                continue
            ctl = r["control"]
            x = j + (0.17 if ctl else -0.17) + rng.uniform(-0.07, 0.07)
            y = r["density_x1e5"] if r["calls"] else zero_y * 10 ** rng.uniform(-0.12, 0.12)
            ax.scatter(x, y, s=9, marker="o", linewidth=0.5, facecolor=(CTRLF if ctl else CASEF),
                       edgecolor=(CTRL if ctl else CASE), zorder=3)
        if c in pv:
            ax.text(j, ymax * 0.85, stars(pv[c]), ha="center", va="top", fontsize=6, color=MID)
    ax.set_yscale("log")
    ax.set_ylim(strip[0], ymax)
    ax.axhspan(strip[0], strip[1], color=SOFT, alpha=0.15, lw=0)
    # break marks on the y axis between the zero strip and the measured range
    for dy in (0.93, 1.07):
        ax.plot([-0.62, -0.38], [strip[1] * dy * 1.08] * 2, color=INK, lw=0.6, clip_on=False, zorder=5)
    decades = [10.0 ** k for k in range(int(np.ceil(np.log10(lo_d / 1.5))), int(np.ceil(np.log10(ymax))))]
    if len(decades) > 4:
        decades = decades[::-1][::2][::-1]
    ax.set_yticks([zero_y] + decades, ["0"] + [f"$10^{{{int(np.log10(d))}}}$" for d in decades])
    ax.minorticks_off()
    ax.set_xlim(-0.6, len(SHOW) - 0.4)
    ax.set_xticks(range(len(SHOW)), [SHORT[c] for c in SHOW], rotation=45, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=2)
    if show_y:
        ax.set_ylabel("Unique variants per 10$^5$ bases")
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
    row1 = outer[0].subgridspec(1, 5, wspace=0.42)
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
        letter(ax, "ABCDE"[i], dx=-0.32, dy=1.2)
    h = [plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc=CASEF, mec=CASE, label="AID-competent library"),
         plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc=CTRLF, mec=CTRL, label="AID-null control library"),
         plt.Rectangle((0, 0), 1, 1, color=SOFT, alpha=0.3, label="0 calls (below axis break)")]
    fig.legend(handles=h, loc="upper center", ncol=3, bbox_to_anchor=(0.55, 0.995), fontsize=6)
    # F: separation fold, one row per case group
    ax = fig.add_subplot(row2[0, 0])
    mat = np.array([[separation(load(k), c, pre)[0] for c in SHOW] for k, _l, pre in avail])
    labs = [lbl.replace("Off-target 1", "CH12 OT1").replace("Off-target 2", "CH12 OT2") for _k, lbl, _p in avail]
    heat_panel(ax, mat, labs, "Blues", True, "Case / control separation", "fold (geometric means)",
               fig.add_subplot(row2[1, 0]), lambda v: f"{v:.0f}" if v >= 10 else f"{v:.1f}")
    letter(ax, "F", dx=-0.62, dy=1.04)
    # G: calls in the AID-null controls, one row per dataset
    ax = fig.add_subplot(row2[0, 1])
    keys = [k for k, _t in SETS if load(k)]
    mat = np.array([[separation(load(k), c)[1] for c in SHOW] for k in keys], dtype=float)
    heat_panel(ax, mat, [SHORTSET[k].replace("Off-target ", "CH12 OT") for k in keys], "Reds", False,
               "Calls in AID-null controls", "calls, summed over 3 controls", fig.add_subplot(row2[1, 1]),
               lambda v: f"{v:.0f}")
    letter(ax, "G", dx=-0.62, dy=1.04)
    # H: leave-one-control-out trace
    ax = fig.add_subplot(row2[:, 2])
    runs = [("jh4_shmtrace", "J$_H$4, 529-bp window", "#2c5d8a", "-"),
            ("jh4_shmtrace_full", "J$_H$4, 1,141 bp", "#2c5d8a", "--"),
            ("n114_shmtrace", "N114 V$_H$", "#6a4c93", "-")]
    dodge = {"jh4_shmtrace": 0.93, "jh4_shmtrace_full": 1.0, "n114_shmtrace": 1.08}   # display offset only
    for run, lab, col, ls in runs:
        cal = json.load(open(FIGS / run / "project.json"))["calibration"]
        tr = cal["qscale_trace"]
        xs = sorted(tr, key=float)
        n = tr[xs[0]][1]
        ax.plot([float(x) * dodge[run] for x in xs], [100 * tr[x][0] / tr[x][1] for x in xs], marker="o", ms=2.5, color=col,
                ls=ls, lw=1, label=f"{lab} ({n:,})")
        sel = cal["qscale"]
        ax.scatter([sel * dodge[run]], [100 * tr[f"{sel:.1f}"][0] / tr[f"{sel:.1f}"][1]], s=34, facecolor="none",
                   edgecolor=col, linewidth=0.9, zorder=4)
    ax.axhline(1, color=CTRL, lw=0.7, ls=":")
    ax.text(0.9, 0.3, "1%: if calibrated", fontsize=6, color=CTRL, va="bottom")
    ax.scatter([], [], s=34, facecolor="none", edgecolor=MID, linewidth=0.9, label="margin selected")
    ax.set_xlabel("Margin on control error rates (×)")
    ax.set_ylabel("Observed held-out tests with p < 0.01 (%)")
    ax.set_xscale("log")
    ticks = sorted({float(x) for run, *_ in runs
                    for x in json.load(open(FIGS / run / "project.json"))["calibration"]["qscale_trace"]})
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(0, 15)
    ax.set_yticks([0, 2, 4, 6, 8, 10])
    ax.set_xlim(0.85, 12)
    ax.set_title("Held-out control calibration", fontsize=6.6)
    ax.legend(fontsize=6, loc="upper right", handlelength=1.8, borderaxespad=0.1, title="held-out tests",
              title_fontsize=6)
    letter(ax, "H", dx=-0.25, dy=1.02)
    json.dump(tests, open(SCRATCH / "fig4_tests.json", "w"), indent=1)
    print("Figure 4A-E Welch p (Holm within panel):",
          {k: {c: f"{v:.3g}" for c, v in t.items()} for k, t in tests.items()})
    save(fig, "Figure4")


def read_fates():
    """The read accounting that was Figure 4I, now Figure S2 (it was S1E-F until it got its own number)."""
    setup()
    fig = plt.figure(figsize=(FULL, 128 * MM))
    ax = fig.add_axes([0.085, 0.56, 0.68, 0.34])
    parts = [("too_long", "over-length", "#5f6b75"), ("low_quality", "expected error > 0.01", "#9aa5b1"),
             ("not_target", "unmapped or paralog", "#d9952b"), ("chimeric", "chimeric", "#c0392b"),
             ("short_alignment", "alignment too short", "#e6b98a"), ("paralogue_alleles", "paralog alleles", "#6a4c93"),
             ("usable", "usable", "#1f6f8b")]
    labels = []
    for run, name in (("jh4_shmtrace", "J$_H$4"), ("n114_shmtrace", "N114 V$_H$")):
        acct = json.load(open(FIGS / run / "table_s1.json"))
        for lib, a in acct.items():
            labels.append(f"{name} {mouse(lib)}")
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
    letter(ax, "A", dx=-0.075, dy=1.12)
    # F: competitive alignment, moved here from the engine figure (S23). Identity of each J_H4 read to the target
    # and to its paralogue, as frozen in scicanvas/figS23_engine.js DATA.p2 (engine_data.py)
    d = json.load(open(Path(__file__).resolve().parents[1] / "scicanvas/s23_paralogue_identity.json"))
    ax = fig.add_axes([0.085, 0.08, 0.27, 0.3])
    # engine_data.py keeps every second target read (DATA.p2) and every paralogue read, of the first 400 reads of
    # the bundled Wt3 example library
    for pts, f, e, lab in ((d["target"], CASEF, CASE, "assigned to the target (every second drawn)"),
                           (d["paralogue"], "#bda3d4", "#7a4f9a", f"assigned to the paralogue (all {len(d['paralogue'])})")):
        ax.scatter([a for a, _b in pts], [b for _a, b in pts], s=8, facecolor=f, edgecolor=e, linewidth=0.5,
                   label=lab, zorder=3)
    lo = min(min(a for a, _b in d["target"] + d["paralogue"]), 0.8)
    ax.plot([lo, 1.0], [lo, 1.0], color=SOFT, lw=0.6, ls="--")
    ax.set_xlim(lo, 1.005)
    ax.set_ylim(lo, 1.005)
    ax.set_xlabel("Identity to the J$_H$4 target")
    ax.set_ylabel("Identity to its paralogue")
    ax.set_title("Each read goes to its best match", fontsize=6.6)
    ax.legend(fontsize=6, loc="upper left", handletextpad=0.2, borderaxespad=0.2, markerscale=1.3)
    letter(ax, "B", dx=-0.27, dy=1.06)
    fig.text(0.44, 0.3, "Every read is aligned to the target, its paralogue and any contaminant at once,\n"
             "and goes to whichever it matches best: reads below the diagonal are closer to\n"
             "the target. A paralogue can also be analysed as a target of its own.",
             fontsize=6.5, color=MID, va="center", linespacing=1.5)
    save(fig, "FigureS1_read_fates")


if __name__ == "__main__":
    main()
    read_fates()
