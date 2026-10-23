"""Figure S8 (Claude revision): downstream inference in simulation, and UMI consensus.

Page 1 (FigureS8_downstream): what the read-library route's estimates and group tests deliver once calls are made.
A  lineage simulation on J_H4 control reads: recall by read support (was Figure 3E)
B  clone-weighted frequency recovers the truth when one unmutated clone dominates (was Figure 3F)
C  quantification: estimate / truth against depth (was Figure 3G)
D  single-read estimate: change in the estimate against one-read mutations planted (was Figure S8A)
E  group tests: false-positive rate with no true difference, power with one (was Figure 3H)
F  power by animals per group and size of difference (was Figure S8B)
Page 2 (FigureS8_umi): UMI consensus, simulated (G, was S8C) and on real dual-UMI reads (H, was S8D).
The drawing code is moved from fig3_calibration.py and figS8_estimates.py; the data and statistics are unchanged.
"""
import json
from collections import defaultdict

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import SCRATCH  # noqa: E402
from style import CALLERS, CTRL, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]


def lineage_panel(ax, lin, runs):
    bins = ["1", "2", "3-5", "6-10", ">10"]
    for scen, ls, dx0 in (("polyclonal", "-", -0.05), ("low", "--", 0.05)):
        for mode in ("recommended", "geneious"):
            r = next(x for x in lin if x["scenario"] == scen and x["mode"] == mode)
            main = mode == "recommended"
            ax.plot(range(5), [r["recall_by_support"][b] for b in bins], ls=ls, marker=CALLERS[mode][2], ms=3,
                    color=CALLERS[mode][1], lw=1.3 if main else 0.8)
            per = [x["recall_by_support"] for x in runs if x["scenario"] == scen and x["mode"] == mode]
            lo = [min(p[b]["recall"] for p in per if p[b]["planted"]) for b in bins]
            hi = [max(p[b]["recall"] for p in per if p[b]["planted"]) for b in bins]
            dx = dx0 + (0.02 if mode == "geneious" else -0.02)
            ax.vlines([i + dx for i in range(5)], lo, hi, color=CALLERS[mode][1], lw=0.6, alpha=0.8)
    ax.set_xticks(range(5), ["1", "2", "3–5", "6–10", ">10"])
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads carrying the mutation")
    ax.set_ylabel("Recall")
    ax.set_title("Simulated lineages, J$_H$4 (3 libraries each)", fontsize=7)
    nulls = {m: next(x for x in lin if x["scenario"] == "null" and x["mode"] == m)["false_calls_per_library"]
             for m in ("recommended", "geneious")}
    ax.text(0.98, 0.04, f"all-null libraries, false calls\nper library: SHMtrace {nulls['recommended']:.1f},\n"
            f"Geneious {nulls['geneious']:.1f}", transform=ax.transAxes, ha="right", fontsize=6, color=MID)
    h = [plt.Line2D([], [], color=CALLERS[m][1], marker=CALLERS[m][2], ms=3, lw=1.3 if m == "recommended" else 0.8,
                    label=CALLERS[m][0]) for m in ("recommended", "geneious")]
    h += [plt.Line2D([], [], color=SOFT, ls="-", label="polyclonal"),
          plt.Line2D([], [], color=SOFT, ls="--", label="low load")]
    ax.legend(handles=h, loc="center right", bbox_to_anchor=(1.0, 0.45), fontsize=6, handlelength=1.6)


def clone_panel(ax, lin, runs):
    scen = ["polyclonal", "expanded", "low"]
    xs = np.arange(len(scen))
    for j, (kind, col) in enumerate((("read", "#9fc0d8"), ("clone", "#274b7e"))):
        tru = [next(x for x in lin if x["scenario"] == s and x["mode"] == "recommended")[f"true_{kind}_weighted_x1e3"]
               for s in scen]
        est = [next(x for x in lin if x["scenario"] == s and x["mode"] == "recommended")[f"est_{kind}_weighted_x1e3"]
               for s in scen]
        ax.bar(xs + (j - 0.5) * 0.36, est, 0.34, color=col, label=f"{kind}-weighted estimate")
        ax.scatter(xs + (j - 0.5) * 0.36, tru, marker="_", s=60, color=INK, linewidth=1.2, zorder=3)
        for i, sc in enumerate(scen):
            rs = [x for x in runs if x["scenario"] == sc and x["mode"] == "recommended"]
            px = xs[i] + (j - 0.5) * 0.36 + np.linspace(-0.09, 0.09, len(rs))
            ax.vlines(px, [x[f"true_{kind}_weighted_x1e3"] for x in rs], [x[f"est_{kind}_weighted_x1e3"] for x in rs],
                      color=INK, lw=0.5, zorder=4)
            ax.scatter(px, [x[f"est_{kind}_weighted_x1e3"] for x in rs], s=3, color="white", edgecolor=INK,
                       linewidth=0.5, zorder=5)
    ax.set_xticks(xs, ["polyclonal", "expanded\nunmutated clone", "low load"])
    ax.set_ylim(0, 27)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    ax.set_title("Read- and clone-weighted estimates", fontsize=7)
    ax.scatter([], [], marker="_", s=40, color=INK, label="truth (mean)")
    ax.legend(fontsize=6, loc="upper left", handlelength=1, ncol=1)


def quant_panel(ax):
    q = json.load(open(DATA / "quant_accuracy_jh4/runs.json")) + json.load(open(DATA / "quant_accuracy_n114/runs.json"))
    metrics = {"read-weighted frequency": ("est_read_weighted_x1e3", "true_read_weighted_x1e3", "#2c5d8a"),
               "unique density": ("est_unique_density_x1e5", "true_unique_density_x1e5", "#a3621b"),
               "combined density": ("est_combined_density_x1e5", "true_unique_density_x1e5", "#1f6f8b")}
    depths = [125, 250, 500, 1000]
    for amp, scen_, ls in (("JH4", "polyclonal", "-"), ("N114", "subclonal", "--")):
        for name, (e, t, col) in metrics.items():
            ys, lo, hi = [], [], []
            for d in depths:
                rs = [r for r in q if r["amplicon"] == amp and r["scenario"] == scen_ and r["depth"] == d
                      and r["mode"] == "recommended" and r.get(t)]
                v = [r[e] / r[t] for r in rs]
                ys.append(np.mean(v))
                lo.append(min(v))
                hi.append(max(v))
            ax.plot(depths, ys, ls=ls, marker="o", ms=2.5, color=col, lw=1, label=name if amp == "JH4" else None)
            off = 1.04 if amp == "N114" else 0.96
            ax.vlines([d * off for d in depths], lo, hi, color=col, lw=0.6, alpha=0.8)
    ax.axhline(1, color=SOFT, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(depths, [str(d) for d in depths])
    ax.minorticks_off()
    ax.set_ylim(0.3, 1.2)
    ax.set_xlabel("Reads analyzed")
    ax.set_ylabel("Estimate / truth")
    ax.set_title("Quantification against depth", fontsize=7)
    ax.plot([], [], color=SOFT, ls="-", label="J$_H$4 polyclonal")
    ax.plot([], [], color=SOFT, ls="--", label="N114 V$_H$ subclonal")
    ax.legend(fontsize=6, loc="lower left", handlelength=1.6, ncol=1)


def single_read_panel(ax):
    rows = json.load(open(DATA / "single_read_spike.json"))
    by = defaultdict(list)
    for r in rows:
        by[r["planted"]].append(r["delta"])
    ks = sorted(by)
    for k in ks:
        ax.scatter([k] * len(by[k]), by[k], s=10, color="#1f6f8b", zorder=3)
    ax.plot([0, max(ks) * 1.05], [0, max(ks) * 1.05], color=SOFT, ls="--", lw=0.7)
    rec = np.mean([d / k for k in ks if k for d in by[k]])
    ax.text(0.04, 0.9, f"mean recovery {100 * rec:.0f}%\ndashed: identity", transform=ax.transAxes, fontsize=6,
            va="top")
    ax.set_xlabel("One-read mutations planted")
    ax.set_ylabel("Change in single-read estimate")
    ax.set_title("Aggregate single-read estimate", fontsize=7)


def group_panel(ax):
    """Rejection rates: with no true difference these are false-positive rates (calibration, target 5%); with a
    two- or four-fold difference they are power. The two questions are drawn as two separate blocks."""
    rep = json.load(open(CAL / "replicates_calibration/runs.json"))
    tests = [("p_pooled_reads", "pooled reads (not animals)", "#c0392b"), ("p_welch_log", "Welch, log", "#2c5d8a"),
             ("p_negbinomial", "negative binomial", "#1f6f8b"), ("p_permutation", "permutation", "#9aa5b1")]
    cells = [(3, 1.0), (5, 1.0), None, (3, 2.0), (3, 4.0), (5, 2.0), (5, 4.0)]
    w = 0.19
    n_runs = set()
    for j, (k, lab, col) in enumerate(tests):
        xs, vals = [], []
        for x, cell in enumerate(cells):
            if cell is None:
                continue
            n, f = cell
            sel = [r[k] < 0.05 for r in rep if r["n"] == n and r["fold"] == f]
            n_runs.add(len(sel))
            xs.append(x + (j - 1.5) * w)
            vals.append(np.mean(sel))
        ax.bar(xs, vals, w * 0.95, color=col, label=lab)
    ax.axvspan(-0.5, 1.5, color=SOFT, alpha=0.1, lw=0)
    ax.plot([-0.5, 1.5], [0.05, 0.05], color=INK, lw=0.6, ls=":")
    ax.text(1.45, 0.07, "5%", fontsize=6, color=INK, ha="right", va="bottom")
    ax.text(0.5, 1.04, "no true difference:\nfalse-positive rate", ha="center", va="bottom", fontsize=6, color=INK)
    ax.text(4.5, 1.04, "true difference: power", ha="center", va="bottom", fontsize=6, color=INK)
    ax.set_xticks([0, 1, 3, 4, 5, 6], ["3 vs 3", "5 vs 5", "2×\n3 vs 3", "4×\n3 vs 3", "2×\n5 vs 5", "4×\n5 vs 5"],
                  fontsize=6)
    ax.set_xlim(-0.5, 6.5)
    ax.set_ylabel("Share of experiments, p < 0.05")
    ax.set_xlabel(f"Animals per group, true fold difference ({'/'.join(map(str, sorted(n_runs)))} simulations each)")
    ax.set_ylim(0, 1.55)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title("Group tests: animals, not reads, as the unit", fontsize=7)
    ax.legend(fontsize=6, ncol=2, loc="upper right", bbox_to_anchor=(1.0, 1.02), columnspacing=0.8,
              handlelength=1.2, borderaxespad=0.1)


def power_panel(ax):
    pw = json.load(open(DATA / "power_s3.json"))
    animals = pw["animals"]
    folds = pw["folds"]
    shades = {"2.0": "#9aa5b1", "5.0": "#6d9ad0", "10.0": "#274b7e", "50.0": "#1f6f8b"}
    for f, col in shades.items():
        if f not in folds["power"]:
            continue
        ax.plot(animals, folds["power"][f], marker="o", ms=2.4, color=col, lw=1, label=f"{float(f):g}-fold")
    ax.axhline(0.8, color=INK, ls=":", lw=0.6)
    ax.set_xscale("log")
    ax.set_xticks(animals, [str(a) for a in animals], fontsize=6)
    ax.minorticks_off()
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Animals per group")
    ax.set_ylabel("Power at α = 0.05")
    ax.set_title("Power by size of difference", fontsize=7, pad=13)
    ax.text(0.0, 1.005, f"unique density, CV {folds['cv']:.2f}; this paper's contrasts are 38–141-fold",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=MID)
    ax.legend(fontsize=6, loc="lower right", handlelength=1.4, borderaxespad=0.2, labelspacing=0.3)


def downstream():
    fig = plt.figure(figsize=(FULL, 150 * MM))
    gs = gridspec.GridSpec(3, 2, figure=fig, left=0.085, right=0.975, top=0.955, bottom=0.07, wspace=0.3,
                           hspace=0.62)
    lin = json.load(open(DATA / "lineage_sim_jh4/lineage_sim.json"))
    summary, runs = lin["summary"], lin["runs"]
    panels = [(lineage_panel, (summary, runs)), (clone_panel, (summary, runs)), (quant_panel, ()),
              (single_read_panel, ()), (group_panel, ()), (power_panel, ())]
    for i, (fn, args) in enumerate(panels):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        fn(ax, *args)
        letter(ax, "ABCDEF"[i], dx=-0.14, dy=1.06)
    save(fig, "FigureS8_downstream")


def umi():
    fig = plt.figure(figsize=(FULL, 58 * MM))
    gs = gridspec.GridSpec(1, 2, figure=fig, left=0.085, right=0.975, top=0.88, bottom=0.2, wspace=0.3)
    # G
    ax = fig.add_subplot(gs[0])
    runs = json.load(open(CAL / "umi_sim/runs.json"))
    cells = [("JH4", "polyclonal"), ("JH4", "low"), ("N114", "subclonal"), ("N114", "low")]
    x = np.arange(len(cells))
    for j, (kind, col, lab) in enumerate((("reads", "#9aa5b1", "every read"), ("umi", "#1f6f8b", "UMI consensus"))):
        fp = [np.mean([r[kind]["false_calls"] for r in runs if (r["amplicon"], r["scenario"]) == c]) for c in cells]
        ax.bar(x + (j - 0.5) * 0.38, fp, 0.36, color=col, label=lab)
        if kind == "umi":
            for xi, v in zip(x, fp, strict=True):
                ax.text(xi + 0.19, v + 0.25, "0" if v == 0 else f"{v:.1f}", ha="center", va="bottom",
                        fontsize=6, color=col)
    ax.set_xticks(x, ["J$_H$4\npolyclonal", "J$_H$4\nlow load", "N114 V$_H$\nsubclonal", "N114 V$_H$\nlow load"],
                  fontsize=6)
    ax.set_ylabel("False calls per library")
    ax.set_title("UMI consensus, simulated libraries", fontsize=7)
    ax.legend(fontsize=6, loc="upper right")
    letter(ax, "G", dx=-0.14, dy=1.04)
    # H
    ax = fig.add_subplot(gs[1])
    u = json.load(open(CAL / "umi_karst/results.json"))["summary"]
    bins = ["3-4", "5-9", "10-19", "20-+"]
    x = np.arange(len(bins))
    allerr = [u["by_reads_per_molecule"][b]["error_pct"] for b in bins]
    sub = [u["by_reads_per_molecule"][b]["substitution_error_pct"] for b in bins]
    ax.plot(x, allerr, marker="o", ms=3, color="#9aa5b1", lw=1.1, label="all errors")
    ax.plot(x, sub, marker="o", ms=3, color="#1f6f8b", lw=1.1, label="substitutions")
    ax.axhline(u["raw_read_error_pct"], color=CTRL, ls="--", lw=0.8)
    ax.text(3.0, u["raw_read_error_pct"] * 1.25, f"raw reads {u['raw_read_error_pct']:.1f}%", fontsize=6,
            color=CTRL, ha="right")
    ax.set_yscale("log")
    ax.set_ylim(0.0005, 30)
    ax.set_xticks(x, [f"{lab}\nn = {u['by_reads_per_molecule'][b]['molecules']}"
                      for lab, b in zip(["3–4", "5–9", "10–19", "≥20"], bins, strict=True)])
    ax.set_xlabel("Reads per molecule (n, molecules)")
    ax.set_ylabel("Error per base (%)")
    mol = [m for m in json.load(open(CAL / "umi_karst/results.json"))["molecules"] if m["reads"] >= 3]
    nm, bl = np.array([m["nm"] for m in mol]), np.array([m["blen"] for m in mol])
    rng = np.random.default_rng(7)
    boot = [100 * nm[i].sum() / bl[i].sum() for i in (rng.integers(0, len(mol), len(mol)) for _ in range(2000))]
    cons = 100 * nm.sum() / bl.sum()
    lo, hi = np.percentile(boot, [2.5, 97.5])
    raw = u["raw_read_error_pct"]
    ax.errorbar([3.45], [cons], yerr=[[cons - lo], [hi - cons]], fmt="D", ms=3, color=INK, capsize=1.5, lw=0.8)
    ax.text(3.3, cons * 2.6, f"all {len(mol)}: {cons:.2f}%\n{raw / cons:.0f}-fold\n({raw / hi:.0f}–{raw / lo:.0f})",
            fontsize=6, color=INK, ha="center", va="bottom", linespacing=1.15)
    ax.set_xlim(-0.3, 3.75)
    ax.set_title("UMI consensus on real dual-UMI rRNA reads", fontsize=7)
    ax.legend(fontsize=6, loc="lower left", title="after consensus", title_fontsize=6, handlelength=1.4)
    letter(ax, "H", dx=-0.14, dy=1.04)
    json.dump({"consensus_error_pct": round(cons, 4), "ci": [round(lo, 4), round(hi, 4)],
               "fold": round(raw / cons, 1), "fold_ci": [round(raw / hi, 1), round(raw / lo, 1)]},
              open(SCRATCH / "figS8_umi.json", "w"), indent=1)
    save(fig, "FigureS8_umi")


def main():
    setup()
    downstream()
    umi()


if __name__ == "__main__":
    main()
