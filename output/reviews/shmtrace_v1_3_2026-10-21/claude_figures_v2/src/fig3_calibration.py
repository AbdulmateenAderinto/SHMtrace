"""Figure 3: calling performance against planted truth (Claude revision: four panels; the downstream simulations
that were E-H are in Figure S8, drawn by figS8_downstream.py).

A-B  spike-ins into held-out AID-null libraries (J_H4, N114): recall by planted read count, per caller
C    ClairS and DeepSomatic recall against the share of reads carrying a mutation (J_H4 spike-ins and clean plants
     at 5, 10 and 20%), with SHMtrace on the same spike-ins, over the cumulative allele-fraction distribution of
     SHMtrace's calls in real libraries (Figure 2E)
D    recall against false calls per library, every caller and dataset with spike-ins
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


def real_call_fractions():
    """Share of reads carrying each SHMtrace call in the AID-active libraries of Figure 2E (%)."""
    af = []
    for key in ("jh4", "n114", "su_lps", "pim1"):
        d = json.load(open(DATA / f"intuition_{key}.json"))
        af += [100 * c[3] / c[4] for n in d["cases"] for c in d["calls"][n] if c[4]]
    return np.sort(af)


def spike_panel(ax, key, show_y=True):
    s = spike_summary(key)
    for j, c in enumerate(SHOW):
        lab, col, mk = CALLERS[c]
        main = c == "recommended"
        ks = list(s[c]["by_k"])
        ax.plot(range(len(ks)), [s[c]["by_k"][k] for k in ks], marker=mk, ms=3.4 if main else 2.6, color=col,
                lw=1.5 if main else 0.8, label=lab, zorder=4 if main else 3)
        # range over the three held-out controls, offset so the callers' bars do not overlap
        spread = spike_by_control(key, c)
        dx = (j - 2) * 0.07
        ax.vlines([i + dx for i in range(len(ks))], [min(spread[k]) for k in ks], [max(spread[k]) for k in ks],
                  color=col, lw=0.6, alpha=0.8)
    ax.set_xticks(range(len(ks)), [str(k) for k in ks])
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads carrying the planted mutation")
    if show_y:
        ax.set_ylabel("Recall")
    ax.set_title(f"{DATASET_LABEL[key]} spike-ins", fontsize=7, pad=12)
    ax.text(0.0, 1.005, f"observed FDP: SHMtrace {s['recommended']['fdr']:.2f}, Geneious {s['geneious']['fdr']:.2f}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=MID)


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 122 * MM))
    outer = gridspec.GridSpec(1, 2, figure=fig, left=0.08, right=0.975, top=0.885, bottom=0.085, wspace=0.24,
                              width_ratios=[1.5, 1])
    left = outer[0].subgridspec(2, 2, wspace=0.3, hspace=0.62)
    # A, B: recall by planted read count
    for i, key in enumerate(("jh4", "n114")):
        ax = fig.add_subplot(left[0, i])
        spike_panel(ax, key, show_y=(i == 0))
        letter(ax, "AB"[i], dx=-0.3 if i == 0 else -0.16, dy=1.1)
    handles = [plt.Line2D([], [], color=CALLERS[c][1], marker=CALLERS[c][2], ms=3.5,
                          lw=1.5 if c == "recommended" else 0.8, label=CALLERS[c][0]) for c in SHOW]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0), fontsize=6,
               columnspacing=1.6, handlelength=1.8)

    # C: the somatic callers need a mutation in a large share of reads; real SHMtrace calls mostly sit below that
    ax = fig.add_subplot(left[1, :])
    af = real_call_fractions()
    xs = np.logspace(np.log10(0.04), np.log10(60), 200)
    cum = np.searchsorted(af, xs, side="right") / len(af)
    ax.fill_between(xs, 0, cum, color=SOFT, alpha=0.18, lw=0, step="post", zorder=1)
    ax.plot(xs, cum, color=SOFT, lw=0.8, drawstyle="steps-post", zorder=2)
    ax.text(0.045, 1.04, f"shaded: cumulative\nshare of {len(af):,} real\nSHMtrace calls (Fig. 2E)",
            fontsize=6, color=MID, ha="left", va="top")
    spiked = spike_summary("jh4")["recommended"]["by_k"]
    lab, col, mk = CALLERS["recommended"]
    ax.plot([100 * k / SPIKE_READS for k in spiked], list(spiked.values()), marker=mk, ms=3.4, color=col, lw=1.5,
            zorder=4)
    ax.text(5.0, 1.0, "SHMtrace", fontsize=6.5, color=col, va="center", ha="left", fontweight="bold")
    for caller, (tx, ty, ha) in (("clairs", (12, 0.07, "right")), ("deepsomatic", (3.4, 0.72, "right"))):
        spike, clean = by_fraction(caller)
        lab, col, mk = CALLERS[caller]
        ax.plot(list(spike), list(spike.values()), marker=mk, ms=3, color=col, lw=0.9, zorder=3)
        ax.plot(list(clean), list(clean.values()), marker=mk, ms=3, color=col, lw=0.9, ls="--", mfc="white",
                zorder=3)
        ax.text(tx, ty, lab, fontsize=6.5, color=col, ha=ha, fontweight="bold")
    ax.text(24, 0.5, "dashed: clean\nplants at\n5, 10, 20%", fontsize=6, color=MID, ha="left", va="center")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 1, 5, 20, 50], ["0.1", "0.2", "1", "5", "20", "50"])
    ax.minorticks_off()
    ax.set_xlim(0.04, 60)
    ax.set_ylim(-0.03, 1.08)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Reads carrying the mutation (%)")
    ax.set_ylabel("Recall, or cumulative share")
    ax.set_title("Somatic callers on J$_H$4, against where real calls fall", fontsize=7)
    letter(ax, "C", dx=-0.135, dy=1.06)

    # D: recall against false calls, the trade-off the panels on the left resolve by read support
    ax = fig.add_subplot(outer[1])
    for key, mk in DATASET_MARK.items():
        s = spike_summary(key)
        if not s:
            continue
        for c in SHOW:
            main = c == "recommended"
            ax.scatter(s[c]["fp"], s[c]["recall"], marker=mk, s=40 if main else 26, color=CALLERS[c][1],
                       edgecolor=INK if main else "white", linewidth=0.6, zorder=4 if main else 3)
    ax.set_xscale("symlog", linthresh=0.1, linscale=0.3)
    ax.set_xlim(-0.02, 20)
    ax.axvspan(-0.02, 0.1, color=SOFT, alpha=0.1, lw=0)
    ax.set_xticks([0, 0.1, 1, 10], ["0", "0.1", "1", "10"])
    ax.minorticks_off()
    ax.set_ylim(-0.03, 1.03)
    ax.text(0.04, 1.0, "linear\n0–0.1", fontsize=6, color=MID, ha="center", va="top")
    ax.text(1.0, 1.0, "logarithmic above 0.1", fontsize=6, color=MID, ha="center", va="top")
    ax.set_xlabel("False calls per library")
    ax.set_ylabel("Recall, all planted mutations")
    ax.set_title("Recall against false calls, four loci", fontsize=7)
    for j, (key, mk) in enumerate(DATASET_MARK.items()):
        if spike_summary(key):
            ax.scatter([], [], marker=mk, s=22, color=SOFT, label=DATASET_LABEL[key])
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.4), fontsize=6, handletextpad=0.3, borderaxespad=0.2,
              title="locus (shape)", title_fontsize=6)
    ax.annotate("", xy=(0.16, 0.66), xytext=(0.4, 0.54), xycoords="data", textcoords="data",
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color=MID, mutation_scale=6))
    ax.text(0.4, 0.52, "fewer false calls,\nhigher recall", fontsize=6, color=MID, ha="center", va="top")
    letter(ax, "D", dx=-0.2, dy=1.03)
    save(fig, "Figure3")


if __name__ == "__main__":
    main()
