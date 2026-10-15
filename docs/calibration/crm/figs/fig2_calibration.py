"""Figure 2: calibration against known truth.

A-B  spike-ins into held-out AID-null libraries (J_H4, N114): recall by planted read count, per caller
C    recall against false calls per library, every caller and dataset with spike-ins
D    lineage simulation on J_H4 control reads: recall by read support, and false calls in null libraries
E    clone-weighted frequency recovers the truth when one unmutated clone dominates
F    quantification: estimate / truth against depth (J_H4 polyclonal, N114 subclonal)
G    group tests: share of simulated experiments with p < 0.05, animals as the unit vs pooled reads
H    ClairS recall against the share of reads carrying a mutation: the J_H4 spike-ins (1-20 of 500 reads) and
     clean mutations planted at 5, 10 and 20% of reads, with SHMtrace on the same spike-ins
"""
import json
from collections import defaultdict

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt

from style import CALLERS, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CALIB = DATA.parents[1]          # docs/calibration
SHOW = ["recommended", "geneious", "lofreq", "lofreq_nobaq", "lofreq_permissive"]
DATASET_LABEL = {"jh4": "J$_H$4", "n114": "N114", "su_lps": "Sμ (LPS)", "pim1": "Pim1"}
DATASET_MARK = {"jh4": "o", "n114": "s", "su_lps": "^", "pim1": "D"}


def spike_summary(key):
    f = DATA / f"spikein_{key}.json"
    if not f.exists():
        return None
    rows = json.load(open(f))["rows"]
    out = {}
    for c in SHOW:
        rs = [r for r in rows if r["caller"] == c]
        rec = defaultdict(lambda: [0, 0])
        for r in rs:
            for k, (got, n) in r["recall_by_reads"].items():
                rec[int(k)][0] += got
                rec[int(k)][1] += n
        tp, fp = sum(r["true_calls"] for r in rs), sum(r["false_calls"] for r in rs)
        out[c] = {"by_k": {k: v[0] / v[1] for k, v in sorted(rec.items()) if v[1]},
                  "recall": sum(v[0] for v in rec.values()) / max(1, sum(v[1] for v in rec.values())),
                  "fp": fp / len(rs), "fdr": fp / (tp + fp) if tp + fp else 0.0}
    return out


SPIKE_READS = 500                # reads per spike-in library (spikein.py)


def clairs_by_fraction():
    """{% of reads: ClairS recall}, J_H4: the spike-ins (out/spikein_jh4_clairs.json) and the clean high-fraction
    plants (out/clairs_highaf_jh4.json, 500 reads). Its two settings gave the same recall at every fraction, so
    the defaults are drawn."""
    rec = defaultdict(lambda: [0, 0])
    for f, n in (("spikein_jh4_clairs.json", SPIKE_READS), ("clairs_highaf_jh4.json", SPIKE_READS)):
        for r in json.load(open(DATA / f))["rows"]:
            if r["caller"] != "clairs":
                continue
            for k, (got, tot) in r["recall_by_reads"].items():
                rec[100 * int(k) / n][0] += got
                rec[100 * int(k) / n][1] += tot
    return {x: g / t for x, (g, t) in sorted(rec.items())}


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 128 * MM))
    outer = gridspec.GridSpec(2, 1, figure=fig, left=0.065, right=0.99, top=0.88, bottom=0.08, hspace=0.55)
    top = outer[0].subgridspec(1, 4, wspace=0.55, width_ratios=[1, 1, 1.15, 1])
    bottom = outer[1].subgridspec(1, 4, wspace=0.55, width_ratios=[1, 1, 1.45, 0.95])
    # A, B: recall by planted read count
    for i, key in enumerate(("jh4", "n114")):
        ax = fig.add_subplot(top[i])
        s = spike_summary(key)
        for c in SHOW:
            lab, col, mk = CALLERS[c]
            ks = list(s[c]["by_k"])
            ax.plot(range(len(ks)), [s[c]["by_k"][k] for k in ks], marker=mk, ms=3, color=col, lw=1, label=lab)
        ax.set_xticks(range(len(ks)), [str(k) for k in ks])
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("Reads carrying the planted mutation")
        ax.set_ylabel("Recall")
        ax.set_title(f"{DATASET_LABEL[key]} spike-ins", fontsize=7)
        letter(ax, "AB"[i], dx=-0.32)
    handles = [plt.Line2D([], [], color=CALLERS[c][1], marker=CALLERS[c][2], ms=3.5, lw=1, label=CALLERS[c][0])
               for c in SHOW]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.995), fontsize=6,
               columnspacing=1.6, handlelength=1.8)
    # C: recall vs false calls
    ax = fig.add_subplot(top[2])
    for key, mk in DATASET_MARK.items():
        s = spike_summary(key)
        if not s:
            continue
        for c in SHOW:
            ax.scatter(max(s[c]["fp"], 0.05), s[c]["recall"], marker=mk, s=16, color=CALLERS[c][1],
                       edgecolor="white", linewidth=0.3, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.04, 20)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("False calls per library (log)")
    ax.set_ylabel("Recall, all planted mutations")
    ax.set_title("Recall against false calls", fontsize=7)
    for j, (key, mk) in enumerate(DATASET_MARK.items()):
        if spike_summary(key):
            ax.scatter([], [], marker=mk, s=14, color=SOFT, label=DATASET_LABEL[key])
    ax.legend(loc="lower right", fontsize=5.5, handletextpad=0.2, borderaxespad=0.2)
    ax.text(0.045, 0.02, "0 = plotted at 0.05", fontsize=5, color=SOFT, transform=ax.transData)
    letter(ax, "C", dx=-0.26)
    # D: lineage simulation recall by support, polyclonal and low
    lin = json.load(open(DATA / "lineage_sim_jh4/lineage_sim.json"))["summary"]
    ax = fig.add_subplot(top[3])
    bins = ["1", "2", "3-5", "6-10", ">10"]
    for scen, ls in (("polyclonal", "-"), ("low", "--")):
        for mode in ("recommended", "geneious"):
            r = next(x for x in lin if x["scenario"] == scen and x["mode"] == mode)
            ax.plot(range(5), [r["recall_by_support"][b] for b in bins], ls=ls, marker=CALLERS[mode][2], ms=3,
                    color=CALLERS[mode][1], lw=1)
    ax.set_xticks(range(5), ["1", "2", "3–5", "6–10", ">10"])
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads carrying the mutation")
    ax.set_ylabel("Recall")
    ax.set_title("Simulated lineages", fontsize=7)
    nulls = {m: next(x for x in lin if x["scenario"] == "null" and x["mode"] == m)["false_calls_per_library"]
             for m in ("recommended", "geneious")}
    ax.text(0.98, 0.04, f"null libraries, false calls:\nSHMtrace {nulls['recommended']:.1f}, "
            f"Geneious {nulls['geneious']:.1f}", transform=ax.transAxes, ha="right", fontsize=5.3, color=MID)
    ax.plot([], [], color=SOFT, ls="-", label="polyclonal")
    ax.plot([], [], color=SOFT, ls="--", label="low load")
    ax.legend(loc="center right", fontsize=5.5, handlelength=1.6)
    letter(ax, "D", dx=-0.3)
    # E: clone-weighted vs read-weighted, estimate against truth (recommended)
    ax = fig.add_subplot(bottom[0])
    scen = ["polyclonal", "expanded", "low"]
    xs = np.arange(len(scen))
    for j, (kind, col) in enumerate((("read", "#9fc0d8"), ("clone", "#274b7e"))):
        tru = [next(x for x in lin if x["scenario"] == s and x["mode"] == "recommended")[f"true_{kind}_weighted_x1e3"] for s in scen]
        est = [next(x for x in lin if x["scenario"] == s and x["mode"] == "recommended")[f"est_{kind}_weighted_x1e3"] for s in scen]
        ax.bar(xs + (j - 0.5) * 0.36, est, 0.34, color=col, label=f"{kind}-weighted")
        ax.scatter(xs + (j - 0.5) * 0.36, tru, marker="_", s=60, color=INK, linewidth=1.2, zorder=3)
    ax.set_xticks(xs, ["poly-\nclonal", "expanded\nclone", "low\nload"])
    ax.set_ylim(0, 21)
    ax.set_ylabel("Mutations per 10$^3$ bp")
    ax.set_title("Clones: estimate and truth", fontsize=7)
    ax.scatter([], [], marker="_", s=40, color=INK, label="truth")
    ax.legend(fontsize=5.5, loc="upper left", handlelength=1, ncol=1)
    letter(ax, "E", dx=-0.32)
    # F: quantification accuracy against depth
    ax = fig.add_subplot(bottom[1])
    q = json.load(open(DATA / "quant_accuracy_jh4/runs.json")) + json.load(open(DATA / "quant_accuracy_n114/runs.json"))
    metrics = {"read-weighted frequency": ("est_read_weighted_x1e3", "true_read_weighted_x1e3", "#2c5d8a"),
               "unique density": ("est_unique_density_x1e5", "true_unique_density_x1e5", "#a3621b"),
               "combined density": ("est_combined_density_x1e5", "true_unique_density_x1e5", "#1f6f8b")}
    depths = [125, 250, 500, 1000]
    for amp, scen_, ls in (("JH4", "polyclonal", "-"), ("N114", "subclonal", "--")):
        for name, (e, t, col) in metrics.items():
            ys = []
            for d in depths:
                rs = [r for r in q if r["amplicon"] == amp and r["scenario"] == scen_ and r["depth"] == d
                      and r["mode"] == "recommended" and r.get(t)]
                ys.append(np.mean([r[e] / r[t] for r in rs]))
            ax.plot(depths, ys, ls=ls, marker="o", ms=2.5, color=col, lw=1, label=name if amp == "JH4" else None)
    ax.axhline(1, color=SOFT, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(depths, [str(d) for d in depths])
    ax.minorticks_off()
    ax.set_ylim(0.3, 1.2)
    ax.set_xlabel("Reads analyzed")
    ax.set_ylabel("Estimate / truth")
    ax.set_title("Quantification", fontsize=7)
    ax.plot([], [], color=SOFT, ls="-", label="J$_H$4")
    ax.plot([], [], color=SOFT, ls="--", label="N114")
    ax.legend(fontsize=5.2, loc="lower left", handlelength=1.6, ncol=1)
    letter(ax, "F", dx=-0.32)
    # G: group tests
    ax = fig.add_subplot(bottom[2])
    rep = json.load(open(CALIB / "replicates_calibration/runs.json"))
    tests = [("p_pooled_reads", "pooled reads", "#c0392b"), ("p_welch_log", "Welch, log", "#2c5d8a"),
             ("p_negbinomial", "negative binomial", "#1f6f8b"), ("p_permutation", "permutation", "#9aa5b1")]
    cells = [(3, 1.0), (3, 2.0), (3, 4.0), (5, 1.0), (5, 2.0), (5, 4.0)]
    w = 0.19
    for j, (k, lab, col) in enumerate(tests):
        vals = [np.mean([x[k] < 0.05 for x in rep if x["n"] == n and x["fold"] == f]) for n, f in cells]
        ax.bar(np.arange(len(cells)) + (j - 1.5) * w, vals, w * 0.95, color=col, label=lab)
    ax.axhline(0.05, color=INK, lw=0.6, ls=":")
    ax.set_xticks(range(len(cells)), ["none" if f == 1 else f"{f:.0f}-fold" for n, f in cells], fontsize=5.5)
    for n, x in ((3, 1), (5, 4)):
        ax.text(x, -0.2, f"{n} vs {n} animals", ha="center", fontsize=5.8, color=MID,
                transform=ax.get_xaxis_transform())
    ax.set_ylabel("Share of experiments, p < 0.05")
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title("Group tests (300 simulations)", fontsize=7)
    ax.legend(fontsize=5.3, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.02), columnspacing=0.8, handlelength=1.2)
    letter(ax, "G", dx=-0.2)
    # H: ClairS needs a mutation in a large share of reads
    ax = fig.add_subplot(bottom[3])
    spiked = spike_summary("jh4")["recommended"]["by_k"]
    cs = clairs_by_fraction()
    lab, col, mk = CALLERS["recommended"]
    ax.plot([100 * k / SPIKE_READS for k in spiked], list(spiked.values()), marker=mk, ms=3, color=col, lw=1)
    ax.text(0.5, 0.42, lab, fontsize=5.8, color=col)
    lab, col, mk = CALLERS["clairs"]
    ax.plot(list(cs), list(cs.values()), marker=mk, ms=3.5, color=col, lw=1)
    ax.text(11, 0.62, lab, fontsize=5.8, color=col, ha="right")
    ax.axvspan(100 / SPIKE_READS * 0.8, 100 * 20 / SPIKE_READS * 1.15, color=SOFT, alpha=0.12, lw=0)
    ax.text(1, 0.1, "spike-ins (A)", fontsize=5.2, color=MID, ha="center")
    ax.set_xscale("log")
    ax.set_xticks([0.2, 1, 5, 20], ["0.2", "1", "5", "20"])
    ax.minorticks_off()
    ax.set_xlim(0.15, 30)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads with the mutation (%)")
    ax.set_ylabel("Recall")
    ax.set_title("ClairS on J$_H$4", fontsize=7)
    letter(ax, "H", dx=-0.32)
    save(fig, "Figure2")


if __name__ == "__main__":
    main()
