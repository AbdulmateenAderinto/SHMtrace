"""Figure S7, page 2: the Sanger pipeline's bookkeeping checks on simulated projects with planted filing errors.

This was Figure 5E; it moved here when panel E became the per-mouse platform comparison.
"""
import json

import numpy as np
from matplotlib import pyplot as plt
from style import DATA, FULL, MID, MM, save, setup

CAL = DATA.parents[1]
NAMES = {"control_contaminant": "mutated clone filed under a control", "file_copy": "copied trace file",
         "retyped_copy": "re-typed copy under another animal", "misfiled_labelled": "misfiled clone, label intact",
         "lineage_misfiled": "whole lineage misfiled", "misfiled_relabelled": "misfiled clone, label changed too"}


def main():
    setup()
    det = json.load(open(CAL / "sanger_validation_summary.json"))["summary"]["detection"]
    fig, ax = plt.subplots(figsize=(FULL * 0.62, 48 * MM))
    fig.subplots_adjust(left=0.45, right=0.95, top=0.86, bottom=0.27)
    keys = list(NAMES)
    vals = [int(det[k].split("/")[0]) / int(det[k].split("/")[1]) for k in keys]
    y = np.arange(len(keys))
    ax.barh(y, [100 * v for v in vals], color="#a3621b", height=0.62)
    for yy, k, v in zip(y, keys, vals, strict=True):
        ax.text(100 * v + 1.5, yy, det[k], va="center", fontsize=6)
    ax.set_yticks(y, [NAMES[k] for k in keys])
    ax.invert_yaxis()
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Planted errors flagged (%)")
    ax.set_title("Sanger bookkeeping checks, 40 simulated projects", fontsize=7, loc="right")
    ax.text(1.0, -0.42, "no flags on 6,018 correctly filed clones", transform=ax.transAxes, ha="right",
            fontsize=6, color=MID)
    save(fig, "FigureS7_bookkeeping")


if __name__ == "__main__":
    main()
