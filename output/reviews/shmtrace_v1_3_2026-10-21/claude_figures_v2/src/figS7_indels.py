"""Figure S7: indels planted into held-out AID-null libraries (indel_spike.py).

A  recall by the number of reads carrying the planted indel, J_H4 and N114 V_H
B  recall by indel class among events carried by three or more reads, the breakdown the Results quote

The indel counterpart of Figure 3A and 3B, which is why it stands on its own rather than beside the
single-read estimate and the power curves. Recall below p.indel_min_reads is zero by definition rather than
by measurement, so the line stops at the floor instead of running through it and the excluded counts are
drawn open.
"""
import json
from collections import defaultdict

from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import CALLERS, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

DATASET_LABEL = {"jh4": "J$_H$4", "n114": "N114 V$_H$"}
FLOOR = 3                       # p.indel_min_reads: one- and two-read indels are not called
# the classes indel_spike.py plants, in the order they are discussed
KIND_LABEL = [("del1", "1-bp deletion, no homopolymer"), ("del1_hp", "1-bp deletion, in homopolymer"),
              ("ins1", "1-bp insertion"), ("del2", "2-bp deletion"), ("ins2", "2-bp insertion"),
              ("del3", "3-bp deletion"), ("dellong", "deletion \u2265 10 bp")]


def recall(key):
    """Planted-indel recall by read count and by class, with the per-control spread and false calls."""
    d = json.load(open(DATA / f"indel_spike_{key}.json"))
    pooled, per_ctl = defaultdict(lambda: [0, 0]), defaultdict(lambda: defaultdict(lambda: [0, 0]))
    kind, kind_n = defaultdict(lambda: [0, 0]), defaultdict(int)
    for r in d["rows"]:
        for p_ in r["planted"]:
            pooled[p_["reads"]][0] += p_["found"]
            pooled[p_["reads"]][1] += 1
            per_ctl[p_["reads"]][r["control"]][0] += p_["found"]
            per_ctl[p_["reads"]][r["control"]][1] += 1
            if p_["reads"] >= FLOOR:
                kind[p_["kind"]][0] += p_["found"]
                kind[p_["kind"]][1] += 1
                kind_n[p_["kind"]] += 1
    ks = sorted(pooled)
    false = [len(r["false"]) for r in d["rows"]]
    ge3 = [sum(pooled[k][0] for k in ks if k >= FLOOR), sum(pooled[k][1] for k in ks if k >= FLOOR)]
    return {"ks": ks, "recall": [pooled[k][0] / pooled[k][1] for k in ks],
            "spread": [[v[0] / v[1] for v in per_ctl[k].values()] for k in ks],
            "by_kind": {k: v[0] / v[1] for k, v in kind.items()}, "kind_n": dict(kind_n),
            "false_per_library": sum(false) / len(false), "recall_ge3": ge3[0] / ge3[1],
            "draws": d["draws"]}


def by_reads(ax, data):
    teal = CALLERS["recommended"][1]
    notes = []
    for j, (key, ls, mk) in enumerate((("jh4", "-", "o"), ("n114", "--", "s"))):
        r = data[key]
        xs = list(range(len(r["ks"])))
        below = [i for i, k in enumerate(r["ks"]) if k < FLOOR]
        above = [i for i, k in enumerate(r["ks"]) if k >= FLOOR]
        ax.plot([xs[i] for i in above], [r["recall"][i] for i in above], color=teal, ls=ls, marker=mk,
                ms=3.4, lw=1, mfc=teal if j == 0 else "white",
                label=f"{DATASET_LABEL[key]}  {100 * r['recall_ge3']:.0f}% at ≥ {FLOOR}")
        # each held-out control in its own right, rather than a min-max bar over three values
        for i in above:
            for v in r["spread"][i]:
                ax.plot([xs[i] + (j - 0.5) * 0.16], [v], marker=mk, ms=2.0, color=teal, alpha=0.5,
                        markeredgewidth=0, linestyle="none")
        notes.append(f"{DATASET_LABEL[key]} {r['false_per_library']:.2f}".rstrip("0").rstrip("."))
    ks = data["jh4"]["ks"]
    # below the floor nothing is measured: the calling rule (indel_min_reads) does not call an indel there, so the
    # region is marked as outside the rule and no recall is drawn in it
    nb = sum(k < FLOOR for k in ks)
    ax.axvspan(-0.5, nb - 0.5, color=SOFT, alpha=0.15, lw=0)
    ax.text((nb - 1) / 2, 0.5, f"outside the\ncalling rule:\nindels need\n≥ {FLOOR} reads,\nso none are\ncalled here",
            ha="center", va="center", fontsize=6, color=MID, linespacing=1.15)
    ax.set_xticks(range(len(ks)), [str(k) for k in ks])
    ax.set_xlim(-0.5, len(ks) - 0.5)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads carrying the planted indel")
    ax.set_ylabel("Recall")
    ax.set_title("Recall by supporting reads", fontsize=7, loc="left", pad=12)
    ax.legend(loc="lower right", fontsize=6, handlelength=2.0, borderaxespad=0.3, edgecolor="none",
              framealpha=0.9)
    ax.text(0.5, 1.005, "false calls per library: " + ", ".join(notes), transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6, color=MID)


def by_kind(ax, data):
    """Horizontal, so each class keeps its name in full rather than an abbreviation."""
    teal = CALLERS["recommended"][1]
    kinds = [k for k, _ in KIND_LABEL if k in data["jh4"]["by_kind"]]
    h = 0.38
    for j, key in enumerate(("jh4", "n114")):
        r = data[key]
        ys = [i - (j - 0.5) * h for i in range(len(kinds))]
        vals = [r["by_kind"].get(k, 0) for k in kinds]
        ax.barh(ys, vals, h * 0.9, color=teal if j == 0 else "white", edgecolor=teal, linewidth=0.7,
                label=DATASET_LABEL[key], zorder=2)
        for y, k, v in zip(ys, kinds, vals):
            n = r["kind_n"].get(k, 0)
            if n:
                ax.text(v + 0.015, y, str(n), ha="left", va="center", fontsize=6, color=MID, zorder=3)
    ax.set_yticks(range(len(kinds)), [dict(KIND_LABEL)[k] for k in kinds], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.14)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylim(len(kinds) - 0.5, -0.65)
    ax.set_xlabel("Recall at \u2265 3 reads (number beside bar: indels planted)")
    ax.set_title("Recall by indel class", fontsize=7, loc="left", pad=12)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6, handlelength=1.1,
              frameon=False, columnspacing=1.2, borderaxespad=0.1)


def main():
    setup()
    data = {k: recall(k) for k in ("jh4", "n114")}
    fig = plt.figure(figsize=(FULL, 64 * MM))
    gs = gridspec.GridSpec(1, 2, figure=fig, left=0.072, right=0.988, top=0.82, bottom=0.17,
                           wspace=0.95, width_ratios=[1, 1.05])
    ax = fig.add_subplot(gs[0]); by_reads(ax, data); letter(ax, "A", dx=-0.19)
    ax = fig.add_subplot(gs[1]); by_kind(ax, data); letter(ax, "B", dx=-0.42)
    for a in fig.axes:
        a.tick_params(labelcolor=INK)
    save(fig, "FigureS7_indels")


if __name__ == "__main__":
    main()
