"""Figure 3: calibration against known truth.

A-B  spike-ins into held-out AID-null libraries (J_H4, N114): recall by planted read count, per caller
C    recall against false calls per library, every caller and dataset with spike-ins
D    lineage simulation on J_H4 control reads: recall by read support, and false calls in null libraries
E    clone-weighted frequency recovers the truth when one unmutated clone dominates
F    quantification: estimate / truth against depth (J_H4 polyclonal, N114 subclonal)
G    group tests: share of simulated experiments with p < 0.05, animals as the unit vs pooled reads
H    ClairS and DeepSomatic recall against the share of reads carrying a mutation: the J_H4 spike-ins (1-20 of 500 reads) and
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
DATASET_LABEL = {"jh4": "J$_H$4", "n114": "N114 V$_H$", "su_lps": "Sμ (LPS)", "pim1": "Pim1"}
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


def spike_by_control(key, caller):
    """Recall by planted read count for each held-out control (pooled over sites and replicates):
    {k: [recall per control]}, the spread behind the pooled curve."""
    rows = [r for r in json.load(open(DATA / f"spikein_{key}.json"))["rows"] if r["caller"] == caller]
    per = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in rows:
        for k, (got, n) in r["recall_by_reads"].items():
            per[int(k)][r["control"]][0] += got
            per[int(k)][r["control"]][1] += n
    return {k: [g / n for g, n in v.values() if n] for k, v in sorted(per.items())}


SPIKE_READS = 500                # reads per spike-in library (spikein.py)


def by_fraction(caller):
    """J_H4 recall by % of reads for a somatic caller ("clairs" or "deepsomatic") at its permissive setting, the
    more sensitive of its two, as two series: the spike-ins (out/spikein_jh4_<caller>.json) and the clean
    high-fraction plants (out/<caller>_highaf_jh4.json), 500 reads each. [{%: recall}, {%: recall}]"""
    out = []
    for f in (f"spikein_jh4_{caller}.json", f"{caller}_highaf_jh4.json"):
        rec = defaultdict(lambda: [0, 0])
        if (DATA / f).exists():
            for r in json.load(open(DATA / f))["rows"]:
                if r["caller"] != f"{caller}_permissive":
                    continue
                for k, (got, tot) in r["recall_by_reads"].items():
                    rec[100 * int(k) / SPIKE_READS][0] += got
                    rec[100 * int(k) / SPIKE_READS][1] += tot
        out.append({x: g / t for x, (g, t) in sorted(rec.items())})
    return out


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 128 * MM))
    outer = gridspec.GridSpec(2, 1, figure=fig, left=0.065, right=0.975, top=0.87, bottom=0.08, hspace=0.6)
    top = outer[0].subgridspec(1, 4, wspace=0.55, width_ratios=[1, 1, 1.15, 1])
    bottom = outer[1].subgridspec(1, 4, wspace=0.55, width_ratios=[1, 1, 1.35, 1.1])
    # A, B: recall by planted read count
    for i, key in enumerate(("jh4", "n114")):
        ax = fig.add_subplot(top[i])
        s = spike_summary(key)
        for j, c in enumerate(SHOW):
            lab, col, mk = CALLERS[c]
            ks = list(s[c]["by_k"])
            ax.plot(range(len(ks)), [s[c]["by_k"][k] for k in ks], marker=mk, ms=3, color=col, lw=1, label=lab)
            # range over the three held-out controls, offset so the callers' bars do not overlap
            spread = spike_by_control(key, c)
            dx = (j - 2) * 0.07
            ax.vlines([i + dx for i in range(len(ks))], [min(spread[k]) for k in ks], [max(spread[k]) for k in ks],
                      color=col, lw=0.7, alpha=0.8)
        ax.set_xticks(range(len(ks)), [str(k) for k in ks])
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("Reads carrying the planted mutation")
        ax.set_ylabel("Recall")
        ax.set_title(f"{DATASET_LABEL[key]} spike-ins", fontsize=7, pad=12)
        ax.text(0.5, 1.005, f"FDR: SHMtrace {s['recommended']['fdr']:.2f}, Geneious {s['geneious']['fdr']:.2f}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=6, color=MID)
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
                       edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.04, 20)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("False calls per library (log)")
    ax.set_ylabel("Recall, all planted mutations")
    ax.set_title("Recall against false calls", fontsize=7)
    for j, (key, mk) in enumerate(DATASET_MARK.items()):
        if spike_summary(key):
            ax.scatter([], [], marker=mk, s=14, color=SOFT, label=DATASET_LABEL[key])
    ax.legend(loc="lower right", fontsize=6, handletextpad=0.2, borderaxespad=0.2)
    ax.text(0.045, 0.13, "0 plotted\nat 0.05", fontsize=6, color=SOFT, transform=ax.transData)
    letter(ax, "C", dx=-0.26)
    # D: lineage simulation recall by support, polyclonal and low
    lin = json.load(open(DATA / "lineage_sim_jh4/lineage_sim.json"))["summary"]
    ax = fig.add_subplot(top[3])
    bins = ["1", "2", "3-5", "6-10", ">10"]
    runs = json.load(open(DATA / "lineage_sim_jh4/lineage_sim.json"))["runs"]
    for scen, ls, dx0 in (("polyclonal", "-", -0.05), ("low", "--", 0.05)):
        for mode in ("recommended", "geneious"):
            r = next(x for x in lin if x["scenario"] == scen and x["mode"] == mode)
            ax.plot(range(5), [r["recall_by_support"][b] for b in bins], ls=ls, marker=CALLERS[mode][2], ms=3,
                    color=CALLERS[mode][1], lw=1)
            # range over the three simulated libraries
            per = [x["recall_by_support"] for x in runs if x["scenario"] == scen and x["mode"] == mode]
            lo = [min(p[b]["recall"] for p in per if p[b]["planted"]) for b in bins]
            hi = [max(p[b]["recall"] for p in per if p[b]["planted"]) for b in bins]
            dx = dx0 + (0.02 if mode == "geneious" else -0.02)
            ax.vlines([i + dx for i in range(5)], lo, hi, color=CALLERS[mode][1], lw=0.7, alpha=0.8)
    ax.set_xticks(range(5), ["1", "2", "3–5", "6–10", ">10"])
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads carrying the mutation")
    ax.set_ylabel("Recall")
    ax.set_title("Simulated lineages", fontsize=7)
    nulls = {m: next(x for x in lin if x["scenario"] == "null" and x["mode"] == m)["false_calls_per_library"]
             for m in ("recommended", "geneious")}
    ax.text(0.98, 0.04, f"null libraries, false calls:\nSHMtrace {nulls['recommended']:.1f}\n"
            f"Geneious {nulls['geneious']:.1f}", transform=ax.transAxes, ha="right", fontsize=6, color=MID)
    ax.plot([], [], color=SOFT, ls="-", label="polyclonal")
    ax.plot([], [], color=SOFT, ls="--", label="low load")
    ax.legend(loc="center right", fontsize=6, handlelength=1.6)
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
        # each simulated library: a stem from its own truth to its estimate, since the truth varies by library
        for i, sc in enumerate(scen):
            rs = [x for x in runs if x["scenario"] == sc and x["mode"] == "recommended"]
            px = xs[i] + (j - 0.5) * 0.36 + np.linspace(-0.09, 0.09, len(rs))
            ax.vlines(px, [x[f"true_{kind}_weighted_x1e3"] for x in rs], [x[f"est_{kind}_weighted_x1e3"] for x in rs],
                      color=INK, lw=0.5, zorder=4)
            ax.scatter(px, [x[f"est_{kind}_weighted_x1e3"] for x in rs], s=3, color="white", edgecolor=INK,
                       linewidth=0.5, zorder=5)
    ax.set_xticks(xs, ["poly-\nclonal", "expanded\nclone", "low\nload"])
    ax.set_ylim(0, 27)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    ax.set_title("Clones: estimate and truth", fontsize=7)
    ax.scatter([], [], marker="_", s=40, color=INK, label="truth")
    ax.legend(fontsize=6, loc="upper left", handlelength=1, ncol=1)
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
            ys, lo, hi = [], [], []
            for d in depths:
                rs = [r for r in q if r["amplicon"] == amp and r["scenario"] == scen_ and r["depth"] == d
                      and r["mode"] == "recommended" and r.get(t)]
                v = [r[e] / r[t] for r in rs]
                ys.append(np.mean(v))
                lo.append(min(v))
                hi.append(max(v))
            ax.plot(depths, ys, ls=ls, marker="o", ms=2.5, color=col, lw=1, label=name if amp == "JH4" else None)
            # range over the three held-out controls
            off = 1.04 if amp == "N114" else 0.96
            ax.vlines([d * off for d in depths], lo, hi, color=col, lw=0.7, alpha=0.8)
    ax.axhline(1, color=SOFT, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(depths, [str(d) for d in depths])
    ax.minorticks_off()
    ax.set_ylim(0.3, 1.2)
    ax.set_xlabel("Reads analyzed")
    ax.set_ylabel("Estimate / truth")
    ax.set_title("Quantification", fontsize=7)
    ax.plot([], [], color=SOFT, ls="-", label="J$_H$4")
    ax.plot([], [], color=SOFT, ls="--", label="N114 V$_H$")
    ax.legend(fontsize=6, loc="lower left", handlelength=1.6, ncol=1)
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
    ax.set_xticks(range(len(cells)), ["none" if f == 1 else f"{f:.0f}×" for n, f in cells], fontsize=6)
    for n, x in ((3, 1), (5, 4)):
        ax.text(x, -0.2, f"{n} vs {n} animals", ha="center", fontsize=6, color=MID,
                transform=ax.get_xaxis_transform())
    ax.set_ylabel("Share of experiments, p < 0.05")
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title("Group tests (300 simulations)", fontsize=7)
    ax.legend(fontsize=6, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.02), columnspacing=0.8, handlelength=1.2)
    letter(ax, "G", dx=-0.2)
    # H: the somatic callers need a mutation in a large share of reads
    ax = fig.add_subplot(bottom[3])
    spiked = spike_summary("jh4")["recommended"]["by_k"]
    lab, col, mk = CALLERS["recommended"]
    ax.plot([100 * k / SPIKE_READS for k in spiked], list(spiked.values()), marker=mk, ms=3, color=col, lw=1)
    ax.text(2.2, 1.04, lab, fontsize=6, color=col, va="bottom", ha="center")
    for caller, (tx, ty, ha) in (("clairs", (4.4, 0.06, "right")), ("deepsomatic", (1.9, 0.765, "center"))):
        spike, clean = by_fraction(caller)
        lab, col, mk = CALLERS[caller]
        ax.plot(list(spike), list(spike.values()), marker=mk, ms=3.5, color=col, lw=1)
        ax.plot(list(clean), list(clean.values()), marker=mk, ms=3.5, color=col, lw=1, ls="--", mfc="white")
        ax.text(tx, ty, lab, fontsize=6, color=col, ha=ha)
    ax.axvspan(100 / SPIKE_READS * 0.8, 100 * 20 / SPIKE_READS * 1.15, color=SOFT, alpha=0.12, lw=0)
    ax.set_xscale("log")
    ax.set_xticks([0.2, 1, 5, 20], ["0.2", "1", "5", "20"])
    ax.minorticks_off()
    ax.set_xlim(0.15, 35)
    ax.set_ylim(-0.03, 1.12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Reads with the mutation (%)")
    ax.set_ylabel("Recall")
    ax.set_title("Somatic callers on J$_H$4", fontsize=7)
    letter(ax, "H", dx=-0.32)
    save(fig, "Figure3")


if __name__ == "__main__":
    main()
