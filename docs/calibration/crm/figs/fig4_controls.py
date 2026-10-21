"""Figure 4: on real libraries with matched AID-null controls, the controls set the background.

A-E  unique-variant density per library (case blue, AID-null control red) for each caller, five datasets
F    separation: geometric-mean case density over control density (controls without calls at one event)
G    calls made in the AID-null control libraries, summed per dataset
H    leave-one-control-out calibration: share of held-out control tests with p < 0.01 against the margin
I    what SHMtrace did with every delivered read (J_H4 and N114 totals)
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy import stats
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


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 168 * MM))
    outer = gridspec.GridSpec(3, 1, figure=fig, left=0.07, right=0.99, top=0.93, bottom=0.06,
                              hspace=0.62, height_ratios=[1, 0.95, 0.72])
    row1 = outer[0].subgridspec(1, 5, wspace=0.42)
    row2 = outer[1].subgridspec(1, 3, wspace=0.55, width_ratios=[1, 1, 0.85])
    # rows of F and G: one per dataset, except OT1/OT2 which contribute one row per AID-competent group
    avail = [(k, lbl, pre) for k, _t in SETS if load(k)
             for lbl, pre in SPLIT.get(k, [(SHORTSET[k], None)])]
    panels = [(k, t) for k, t in PANEL_SETS if load(k)]
    tests = {}
    for i, (key, title) in enumerate(panels):
        ax = fig.add_subplot(row1[0, i])
        rows = load(key)
        floor = min(1e5 / r["bases"] for r in rows) * 0.5
        rng = np.random.default_rng(1)
        for j, c in enumerate(SHOW):
            for r in rows:
                if r["mode"] != c:
                    continue
                y = r["density_x1e5"] if r["calls"] else floor
                ctl = r["control"]
                x = j + (0.17 if ctl else -0.17) + rng.uniform(-0.07, 0.07)
                ax.scatter(x, y, s=9, marker="o" if r["calls"] else "v", linewidth=0.5,
                           facecolor=(CTRLF if ctl else CASEF), edgecolor=(CTRL if ctl else CASE), zorder=3)
        ax.set_yscale("log")
        ax.axhline(floor * 1.6, color=SOFT, lw=0.5, ls=":")
        # the test each caller's two groups support, marked over its own pair of columns
        pv = panel_tests(rows)
        tests[key] = {c: round(v, 5) for c, v in pv.items()}
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi ** 1.22 / lo ** 0.22)
        top = ax.get_ylim()[1]
        for j, c in enumerate(SHOW):
            if c in pv:
                ax.text(j, top ** 0.995 * lo ** 0.005, stars(pv[c]), ha="center", va="top",
                        fontsize=6, color=MID)
        ax.set_xticks(range(len(SHOW)), [SHORT[c] for c in SHOW], rotation=55, ha="right")
        ax.set_xlim(-0.6, len(SHOW) - 0.4)
        ax.set_title(title, fontsize=6.5)
        if i == 0:
            ax.set_ylabel("Unique variants per 10$^5$ bases")
        letter(ax, "ABCDE"[i], dx=-0.3 if i == 0 else -0.22)
    # legend for case / control
    h = [plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc=CASEF, mec=CASE, label="AID-competent"),
         plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc=CTRLF, mec=CTRL, label="AID-null control"),
         plt.Line2D([], [], marker="v", ls="", ms=3.5, mfc="white", mec=MID, label="no calls (below dotted line)")]
    fig.legend(handles=h, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.995), fontsize=6)
    # F: separation fold, G: control calls (dot/heat tables)
    for p, (idx, title, which) in enumerate(((0, "Separation, case / control", 0), (1, "Calls in AID-null controls", 1))):
        ax = fig.add_subplot(row2[0, idx])
        mat = np.full((len(avail), len(SHOW)), np.nan)
        # G counts calls in a dataset's control libraries, which do not depend on the case group, so a dataset
        # split across two rows in F reports them once, on its first row
        seen = set()
        for i, (key, _lbl, pre) in enumerate(avail):
            for j, c in enumerate(SHOW):
                fold, ctl_calls, n_ctl = separation(load(key), c, pre)
                if which == 0:
                    mat[i, j] = fold
                elif key not in seen:
                    mat[i, j] = ctl_calls
            seen.add(key)
        show = np.log10(np.clip(mat, 1, None)) if which == 0 else mat
        cmap = "Blues" if which == 0 else "Reds"
        ax.imshow(show, cmap=cmap, aspect="auto", vmin=0, vmax=np.nanmax(show) * 1.15)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isnan(v):
                    continue
                txt = (f"{v:.0f}" if v >= 10 else f"{v:.1f}") if which == 0 else f"{v:.0f}"
                dark = show[i, j] > np.nanmax(show) * 0.6
                ax.text(j, i, txt, ha="center", va="center", fontsize=6, color="white" if dark else INK)
        ax.set_xticks(range(len(SHOW)), [SHORT[c] for c in SHOW], rotation=55, ha="right")
        # F is per case group ("OT1 Wt", "OT1 dKO"); G is per dataset, so a split dataset is named once
        if which == 0:
            labs = [lbl for _k, lbl, _p in avail]
        else:
            first, labs = set(), []
            for k, _lbl, _p in avail:
                labs.append("" if k in first else SHORTSET[k])
                first.add(k)
        ax.set_yticks(range(len(avail)), labs)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(title, fontsize=6.5)
        letter(ax, "FG"[p], dx=-0.3)
    # H: leave-one-control-out trace
    ax = fig.add_subplot(row2[0, 2])
    runs = [("jh4_shmtrace", "J$_H$4, 529 bp", "#2c5d8a", "-"), ("jh4_shmtrace_full", "J$_H$4, 1,141 bp", "#2c5d8a", "--"),
            ("n114_shmtrace", "N114 V$_H$", "#6a4c93", "-")]
    for run, lab, col, ls in runs:
        cal = json.load(open(FIGS / run / "project.json"))["calibration"]
        tr = cal["qscale_trace"]
        xs = sorted(tr, key=float)
        ax.plot([float(x) for x in xs], [100 * tr[x][0] / tr[x][1] for x in xs], marker="o", ms=2.5, color=col, ls=ls,
                lw=1, label=lab)
    ax.axhline(1, color=CTRL, lw=0.6, ls=":")
    ax.text(1.0, 0.55, "1% target", fontsize=6, color=CTRL)
    ax.set_xlabel("Margin on control error rates (×)")
    ax.set_ylabel("Held-out control tests, p < 0.01 (%)")
    ax.set_xscale("log")
    ticks = sorted({float(x) for run, *_ in runs for x in json.load(open(FIGS / run / "project.json"))["calibration"]["qscale_trace"]})
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_title("Leave-one-control-out", fontsize=6.5)
    ax.legend(fontsize=6, loc="upper right")
    letter(ax, "H", dx=-0.3)
    # I: read accounting
    ax = fig.add_subplot(outer[2])
    parts = [("too_long", "over-length", "#5f6b75"), ("low_quality", "expected error > 0.01", "#9aa5b1"),
             ("not_target", "unmapped or paralog", "#d9952b"), ("chimeric", "chimeric", "#c0392b"),
             ("short_alignment", "alignment too short", "#e6b98a"), ("paralogue_alleles", "paralog alleles", "#6a4c93"),
             ("usable", "usable", "#1f6f8b")]
    labels = []
    for k, (run, name) in enumerate((("jh4_shmtrace", "J$_H$4"), ("n114_shmtrace", "N114 V$_H$"))):
        acct = json.load(open(FIGS / run / "table_s1.json"))
        for lib, a in acct.items():
            labels.append(f"{name} {lib}")
            left = 0
            for key, lab, col in parts:
                v = 100 * a[key] / a["delivered"]
                ax.bar(len(labels) - 1, v, bottom=left, color=col, width=0.75, label=lab if len(labels) == 1 else None)
                left += v
    ax.set_xticks(range(len(labels)), labels, fontsize=6, rotation=40, ha="right")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylabel("Delivered reads (%)")
    ax.set_title("Fate of every delivered read (SHMtrace read preparation)", fontsize=6.5)
    ax.legend(fontsize=6, loc="center left", bbox_to_anchor=(1.005, 0.5), handlelength=1)
    ax.set_position([0.07, ax.get_position().y0, 0.72, ax.get_position().height])
    letter(ax, "I", dx=-0.075)
    json.dump(tests, open(DATA / "fig4_tests.json", "w"), indent=1)
    print("Figure 4A-E Welch p (Holm within panel):",
          {k: {c: f"{v:.3g}" for c, v in t.items()} for k, t in tests.items()})
    save(fig, "Figure4")


if __name__ == "__main__":
    main()
