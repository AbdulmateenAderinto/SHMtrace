"""Figure S4: the control-based false-discovery threshold against planted truth, at six targets.

A  realised false-discovery proportion (false / all calls, pooled over a dataset's libraries) against the target
B  false calls per library: below the default they stop falling, while recall (C) does not
C  recall of the planted mutations at each target: what a looser or stricter target buys

Libraries are the spike-in libraries of Figure 3 (spikein.py --tag _fdrcal): each held-out AID-null control in turn,
60 substitutions planted into 1-20 reads, three draws, called with only the other controls as the null and decoys.
The targets are fdr_alpha 0.01, 0.02, 0.05 (the default), 0.10, 0.20 and 0.30, everything else as the default.
Writes out/fdr_calibration.json, the numbers the text quotes.
"""
import json

from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

TARGETS = {"recommended_fdr01": 0.01, "recommended_fdr02": 0.02, "recommended": 0.05,
           "recommended_fdr10": 0.10, "recommended_fdr20": 0.20, "recommended_fdr30": 0.30}
SETS = [("jh4", "J$_H$4", "#274b7e", "o"), ("n114", "N114", "#6a4c93", "s"),
        ("su_lps", "Sμ", "#1f6f8b", "^"), ("pim1", "Pim1", "#a3621b", "D")]
ZERO = 0.003                     # where a pooled proportion of 0 is drawn on the log axis


def summary():
    out = {}
    for key, *_ in SETS:
        f = DATA / f"spikein_{key}_fdrcal.json"
        if not f.exists():
            continue
        rows = json.load(open(f))["rows"]
        per = {}
        for caller, t in TARGETS.items():
            rs = [r for r in rows if r["caller"] == caller]
            if not rs:
                continue
            tp, fp = sum(r["true_calls"] for r in rs), sum(r["false_calls"] for r in rs)
            fdps = [r["false_calls"] / (r["true_calls"] + r["false_calls"]) if r["true_calls"] + r["false_calls"] else 0.0
                    for r in rs]
            per[t] = {"libraries": len(rs), "true": tp, "false": fp, "fdp_pooled": fp / (tp + fp) if tp + fp else 0.0,
                      "fdp_libraries": fdps, "over_target": sum(x > t for x in fdps),
                      "false_per_library": fp / len(rs),
                      "fdp_library_mean": sum(fdps) / len(fdps),
                      "recall": sum(r["true_calls"] for r in rs) / sum(r["planted"] for r in rs)}
        out[key] = per
    return out


def main():
    setup()
    s = summary()
    json.dump(s, open(__import__("pathlib").Path(__file__).resolve().parents[2] / "scratch" / "fdr_calibration.json", "w"), indent=1)
    ts = sorted(TARGETS.values())
    X = {t: i for i, t in enumerate(ts)}            # the targets, evenly spaced
    labels = [f"{100 * t:g}%" for t in ts]
    fig = plt.figure(figsize=(FULL, 62 * MM))
    gs = gridspec.GridSpec(1, 3, figure=fig, left=0.085, right=0.985, top=0.9, bottom=0.2, wspace=0.42)

    def axis(ax):
        ax.set_xticks(range(len(ts)), labels)
        ax.set_xlim(-0.4, len(ts) - 0.6)
        ax.axvline(X[0.05], color=SOFT, linewidth=0.5, zorder=0)
        ax.set_xlabel("Target false-discovery rate")

    def series(ax, value, legend=False):
        for key, label, col, mk in SETS:
            if key not in s:
                continue
            xs = [t for t in ts if t in s[key]]
            ax.plot([X[t] for t in xs], [value(s[key][t]) for t in xs], color=col, linewidth=0.9, marker=mk,
                    markersize=3.4, markeredgewidth=0, zorder=3, label=label if legend else None)

    # A: realised against target; the dashed line is where they agree, shaded above it is over target
    ax = fig.add_subplot(gs[0])
    xs = list(range(len(ts)))
    ax.plot(xs, ts, color=INK, linestyle=(0, (3, 2)), linewidth=0.6, zorder=2)
    # the agreement line carried on to the panel's edges at the slope of its end segments (log scale)
    edges = [ts[0] * (ts[0] / ts[1]) ** 0.4, ts[-1] * (ts[-1] / ts[-2]) ** 0.4]
    ax.fill_between([-0.4] + xs + [len(ts) - 0.6], [edges[0]] + ts + [edges[1]], 0.6, color="#f3e1e2", linewidth=0, zorder=0)
    ax.text(-0.25, 0.42, "over target", color="#9a3439", fontsize=6, va="center")
    ax.text(len(ts) - 0.75, ZERO * 1.6, "under target", color=MID, fontsize=6, ha="right", va="center")
    series(ax, lambda v: v["fdp_pooled"])
    ax.set_yscale("symlog", linthresh=0.01, linscale=0.5); ax.set_ylim(-0.001, 0.6)
    yt = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    ax.set_yticks([0] + yt, ["0"] + [f"{100 * t:g}%" for t in yt]); ax.set_yticks([], minor=True)
    axis(ax)
    ax.set_ylabel("Realised false-discovery\nproportion (pooled)")
    letter(ax, "A", dx=-0.3)

    # B: false calls per library, which stop falling below the default: a floor no threshold reaches
    ax = fig.add_subplot(gs[1])
    series(ax, lambda v: v["false_per_library"], legend=True)
    ax.set_yscale("log"); ax.set_ylim(0.1, 40)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_minor_formatter(plt.NullFormatter())
    axis(ax)
    ax.set_ylabel("False calls per library")
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.2, labelspacing=0.3)
    letter(ax, "B", dx=-0.3)

    # C: what each target costs in recall
    ax = fig.add_subplot(gs[2])
    series(ax, lambda v: 100 * v["recall"])
    ax.set_ylim(0, 100); ax.set_yticks([0, 25, 50, 75, 100])
    axis(ax)
    ax.text(X[0.05] + 0.08, 3, "default", color=MID, fontsize=6, rotation=90, va="bottom")
    ax.set_ylabel("Planted mutations recalled (%)")
    letter(ax, "C", dx=-0.3)
    save(fig, "FigureS4_fdr_calibration")


if __name__ == "__main__":
    main()
