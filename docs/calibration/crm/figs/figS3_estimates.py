"""Figure S3: single-read estimate, power, and UMI consensus on simulated libraries.

A  single-read estimate: change in the estimate against one-read mutations planted (J_H4, 3 held-out controls)
B  power at alpha 0.05 to detect a 2-fold change, by animals per group (J_H4 Wt pilot; power_s3.py), with a band
   for the uncertainty of a CV estimated from three mice
C  UMI consensus against every read: false calls and recall of single-molecule mutations (umi_sim)
"""
import json
from collections import defaultdict

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt

from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 58 * MM))
    gs = gridspec.GridSpec(1, 3, figure=fig, left=0.085, right=0.99, top=0.86, bottom=0.2, wspace=0.42)
    # A
    ax = fig.add_subplot(gs[0, 0])
    rows = json.load(open(DATA / "single_read_spike.json"))
    by = defaultdict(list)
    for r in rows:
        by[r["planted"]].append(r["delta"])
    ks = sorted(by)
    for k in ks:
        ax.scatter([k] * len(by[k]), by[k], s=10, color="#1f6f8b", zorder=3)
    ax.plot([0, max(ks) * 1.05], [0, max(ks) * 1.05], color=SOFT, ls="--", lw=0.7)
    rec = np.mean([d / k for k in ks if k for d in by[k]])
    ax.text(0.04, 0.9, f"mean recovery {100 * rec:.0f}%", transform=ax.transAxes, fontsize=5.8)
    ax.set_xlabel("One-read mutations planted")
    ax.set_ylabel("Change in single-read estimate")
    ax.set_title("Single-read estimate", fontsize=7)
    letter(ax, "A", dx=-0.3)
    # B: power from the J_H4 wild-type pilot (power_s3.py), with the pilot CV's own uncertainty
    ax = fig.add_subplot(gs[0, 1])
    pw = json.load(open(DATA / "power_s3.json"))
    animals = pw["animals"]
    cols = {"density": "#1f6f8b", "clone_weighted": "#274b7e", "frequency": "#9aa5b1"}
    d = pw["measures"]["density"]
    ax.fill_between(animals, d["power_cv_double"], d["power_cv_low"], color=cols["density"], alpha=0.15, lw=0)
    for m, col in cols.items():
        r = pw["measures"][m]
        ax.plot(animals, r["power"], marker="o", ms=2.4, color=col, lw=1, label=f"{r['label']} (CV {r['cv']:.2f})")
    # where the unique-density estimate reaches 80%
    n80 = float(np.interp(0.8, d["power"], animals))
    ax.axhline(0.8, color=INK, ls=":", lw=0.6)
    ax.annotate(f"≈{n80:.0f} mice", (n80, 0.8), (n80 * 1.35, 0.6), fontsize=5.2, color=cols["density"],
                arrowprops={"arrowstyle": "-", "lw": 0.5, "color": cols["density"]})
    ax.set_xscale("log")
    ax.set_xticks(animals, [str(a) for a in animals], fontsize=5.4)
    ax.minorticks_off()
    ax.set_ylim(0, 1.32)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Animals per group")
    ax.set_ylabel("Power, 2-fold change")
    ax.set_title("Power (J$_H$4 wild-type pilot)", fontsize=7)
    ax.legend(fontsize=5, loc="upper left", bbox_to_anchor=(0, 1.02), handlelength=1.4)
    letter(ax, "B", dx=-0.3)
    # C
    ax = fig.add_subplot(gs[0, 2])
    runs = json.load(open(CAL / "umi_sim/runs.json"))
    cells = [("JH4", "polyclonal"), ("JH4", "low"), ("N114", "subclonal"), ("N114", "low")]
    x = np.arange(len(cells))
    for j, (kind, col, lab) in enumerate((("reads", "#9aa5b1", "every read"), ("umi", "#1f6f8b", "UMI consensus"))):
        fp = [np.mean([r[kind]["false_calls"] for r in runs if (r["amplicon"], r["scenario"]) == c]) for c in cells]
        ax.bar(x + (j - 0.5) * 0.38, fp, 0.36, color=col, label=lab)
    ax.set_xticks(x, ["J$_H$4\npolyclonal", "J$_H$4\nlow", "N114\nsubclonal", "N114\nlow"])
    ax.set_ylabel("False calls per library")
    ax.set_title("UMI consensus (simulated)", fontsize=7)
    ax.legend(fontsize=5.5, loc="upper right")
    letter(ax, "C", dx=-0.3)
    save(fig, "FigureS3")


if __name__ == "__main__":
    main()
