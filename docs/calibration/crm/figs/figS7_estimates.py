"""Figure S7: single-read estimate, power, and UMI consensus on simulated libraries.

A  single-read estimate: change in the estimate against one-read mutations planted (J_H4, 3 held-out controls)
B  power at alpha 0.05 by animals per group and by the size of the difference sought (J_H4 Wt pilot; power_s3.py)
C  UMI consensus against every read: false calls per library (umi_sim)
D  UMI consensus on real dual-UMI rRNA reads (umi_karst), the measured counterpart of C
E  indels planted into held-out AID-null libraries (indel_spike.py): recall by planted read count, J_H4 and N114
"""
import json
from collections import defaultdict

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import CALLERS, CTRL, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]


DATASET_LABEL = {"jh4": "J$_H$4", "n114": "N114 V$_H$"}


def indel_recall(key):
    """Planted-indel recall by read count, pooled and per held-out control, with the false calls it cost."""
    d = json.load(open(DATA / f"indel_spike_{key}.json"))
    pooled, per_ctl = defaultdict(lambda: [0, 0]), defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in d["rows"]:
        for p_ in r["planted"]:
            pooled[p_["reads"]][0] += p_["found"]
            pooled[p_["reads"]][1] += 1
            per_ctl[p_["reads"]][r["control"]][0] += p_["found"]
            per_ctl[p_["reads"]][r["control"]][1] += 1
    ks = sorted(pooled)
    false = [len(r["false"]) for r in d["rows"]]
    ge3 = [sum(pooled[k][0] for k in ks if k >= 3), sum(pooled[k][1] for k in ks if k >= 3)]
    return {"ks": ks, "recall": [pooled[k][0] / pooled[k][1] for k in ks],
            "spread": [[v[0] / v[1] for v in per_ctl[k].values()] for k in ks],
            "false_per_library": sum(false) / len(false), "recall_ge3": ge3[0] / ge3[1]}


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 104 * MM))
    gs = gridspec.GridSpec(2, 6, figure=fig, left=0.075, right=0.978, top=0.93, bottom=0.10,
                           wspace=1.25, hspace=0.46)
    # A
    ax = fig.add_subplot(gs[0, 0:3])
    rows = json.load(open(DATA / "single_read_spike.json"))
    by = defaultdict(list)
    for r in rows:
        by[r["planted"]].append(r["delta"])
    ks = sorted(by)
    for k in ks:
        ax.scatter([k] * len(by[k]), by[k], s=10, color="#1f6f8b", zorder=3)
    ax.plot([0, max(ks) * 1.05], [0, max(ks) * 1.05], color=SOFT, ls="--", lw=0.7)
    rec = np.mean([d / k for k in ks if k for d in by[k]])
    ax.text(0.04, 0.9, f"mean recovery {100 * rec:.0f}%", transform=ax.transAxes, fontsize=6)
    ax.set_xlabel("One-read mutations planted")
    ax.set_ylabel("Change in single-read estimate")
    ax.set_title("Single-read estimate", fontsize=7)
    letter(ax, "A", dx=-0.13)
    # B: power by effect size. A 2-fold difference between two AID-proficient groups is a far harder question
    # than the AID-versus-null contrasts this paper makes, and reporting only that invites the reading that the
    # method needs tens of animals.
    ax = fig.add_subplot(gs[0, 3:6])
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
    ax.text(0.5, 1.005, f"unique density, CV {folds['cv']:.2f}; this paper's own contrasts are 38–141-fold",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=6, color=MID)
    ax.legend(fontsize=6, loc="lower right", handlelength=1.4, borderaxespad=0.2, labelspacing=0.3)
    letter(ax, "B", dx=-0.13)
    # C
    ax = fig.add_subplot(gs[1, 0:2])
    runs = json.load(open(CAL / "umi_sim/runs.json"))
    cells = [("JH4", "polyclonal"), ("JH4", "low"), ("N114", "subclonal"), ("N114", "low")]
    x = np.arange(len(cells))
    for j, (kind, col, lab) in enumerate((("reads", "#9aa5b1", "every read"), ("umi", "#1f6f8b", "UMI consensus"))):
        fp = [np.mean([r[kind]["false_calls"] for r in runs if (r["amplicon"], r["scenario"]) == c]) for c in cells]
        ax.bar(x + (j - 0.5) * 0.38, fp, 0.36, color=col, label=lab)
        if kind == "umi":
            # three of the four are exactly zero, and a zero-height bar reads as missing data
            for xi, v in zip(x, fp, strict=True):
                ax.text(xi + 0.19, v + 0.25, "0" if v == 0 else f"{v:.1f}", ha="center", va="bottom",
                        fontsize=6, color=col)
    ax.set_xticks(x, ["J$_H$4\npolyclonal", "J$_H$4\nlow load", "N114 V$_H$\nsubclonal", "N114 V$_H$\nlow load"],
                  fontsize=6)
    ax.set_ylabel("False calls per library")
    ax.set_title("UMI consensus (simulated)", fontsize=7)
    ax.legend(fontsize=6, loc="upper right")
    letter(ax, "C", dx=-0.13)

    # D: the same question on real data, beside the simulation in C
    ax = fig.add_subplot(gs[1, 2:4])
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
    # the aggregate reduction, with a bootstrap interval over molecules (the raw-read rate is held fixed)
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
    ax.set_title("UMI consensus on real reads", fontsize=7)
    ax.legend(fontsize=6, loc="lower left", title="consensus", title_fontsize=6, handlelength=1.4)
    json.dump({"consensus_error_pct": round(cons, 4), "ci": [round(lo, 4), round(hi, 4)],
               "fold": round(raw / cons, 1), "fold_ci": [round(raw / hi, 1), round(raw / lo, 1)]},
              open(DATA / "figS7_umi.json", "w"), indent=1)
    letter(ax, "D", dx=-0.13)
    # E: planted indels, the indel counterpart of Figure 3A and 3B
    ax = fig.add_subplot(gs[1, 4:6])
    teal = CALLERS["recommended"][1]
    notes = []
    for j, (key, ls, mk) in enumerate((("jh4", "-", "o"), ("n114", "--", "s"))):
        r = indel_recall(key)
        xs = range(len(r["ks"]))
        ax.plot(xs, r["recall"], color=teal, ls=ls, marker=mk, ms=3, lw=1, mfc=teal if j == 0 else "white",
                label=f"{DATASET_LABEL[key]}  {100 * r['recall_ge3']:.0f}% at \u2265 3")
        ax.vlines([x + (j - 0.5) * 0.12 for x in xs], [min(v) for v in r["spread"]], [max(v) for v in r["spread"]],
                  color=teal, lw=0.7, alpha=0.8)
        notes.append(f"{DATASET_LABEL[key]} {r['false_per_library']:.2f}".rstrip("0").rstrip("."))
    ax.axvspan(-0.5, 1.5, color=SOFT, alpha=0.12, lw=0)
    ax.text(0.5, 0.52, "below the\n3-read floor", ha="center", va="center", fontsize=6, color=MID)
    ax.set_xticks(range(len(r["ks"])), [str(k) for k in r["ks"]])
    ax.set_xlim(-0.5, len(r["ks"]) - 0.5)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Reads with the indel")
    ax.set_ylabel("Recall")
    ax.set_title("Planted indels", fontsize=7, pad=12)
    ax.legend(loc="lower right", fontsize=6, handlelength=2.0, borderaxespad=0.1)
    ax.text(0.5, 1.005, "false calls/library: " + ", ".join(notes), transform=ax.transAxes, ha="center",
            va="bottom", fontsize=6, color=MID)
    letter(ax, "E", dx=-0.13)
    save(fig, "FigureS7")


if __name__ == "__main__":
    main()
