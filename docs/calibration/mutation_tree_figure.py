"""Mutation trees of J꜀4 B-cell clones from the nanopore reads (shmtrace.mutation_tree).

    python docs/calibration/mutation_tree_figure.py

The full-length J꜀4 run (positions 90-1,230, deposited 500-read subsets, recommended mode,
as jh4_shmtrace.py --full-length). For each Wt mouse, the clone (shared V(D)J junction)
with the most tree nodes is drawn: the root is the clone's unmutated ancestor, each node
a set of mutations carried by the same reads, its depth the mutations accumulated since
the root, its area the reads carrying it. Writes figures/mutation_trees/
{mutation_trees.png, .pdf, trees.json}.
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from shmtrace import mutation_tree, pipeline, project  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "src/shmtrace/references"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Processed_FASTQ_Depth_Subsets")
LIBS = {"WT2": "DHJGKB_2_Wildtype_2", "WT3": "DHJGKB_3_Wildtype_3", "WT4": "DHJGKB_4_Wildtype_4",
        "KO1": "DHJGKB_5_AID_KO_1", "KO2": "DHJGKB_6_AID_KO_2", "KO3": "DHJGKB_7_AID_KO_3"}
OUT = Path(__file__).parent / "figures/mutation_trees"
MIN_CLONE_READS = 10
# mutation classes: AID deaminates C (G on the other strand); A:T mutations come from
# error-prone repair of the resulting lesion (MMR with polymerase eta)
CLASS = {"GC_ts": ("G:C transition", "#C0392B"), "GC_tv": ("G:C transversion", "#E67E22"),
         "AT": ("A:T mutation", "#2471A3")}


def mclass(ref, alt):
    if ref in "GC":
        return "GC_ts" if ref + alt in ("GA", "CT") else "GC_tv"
    return "AT"


def layout(tree):
    """y by leaf order (larger branches first), x by mutations accumulated since the root."""
    nodes = {n.id: n for n in tree["nodes"]}
    kids = {0: [n.id for n in tree["nodes"] if n.parent is None]}
    for n in tree["nodes"]:
        kids[n.id] = list(n.children)
    for k in kids:
        kids[k].sort(key=lambda i: -len(nodes[i].reads))
    x, y, nxt = {0: 0}, {}, [0]

    def place(i):
        if i:
            x[i] = mutation_tree.depth(tree, i)
        if not kids[i]:
            y[i] = nxt[0]
            nxt[0] += 1
            return
        for c in kids[i]:
            place(c)
        y[i] = (y[kids[i][0]] + y[kids[i][-1]]) / 2
    place(0)
    return x, y, kids


def draw(ax, tree, ref, title):
    """Horizontal phylogeny: the ancestor at left; each mutation a coloured tick on its branch."""
    nodes = {n.id: n for n in tree["nodes"]}
    x, y, kids = layout(tree)
    total = tree["reads"]
    span = max(x.values()) or 1
    for i, cs in kids.items():
        if cs:
            ax.plot([x[i], x[i]], [min(y[c] for c in cs), max(y[c] for c in cs)], color="#9AA3AD", lw=1.0, zorder=1)
        for c in cs:
            n = nodes[c]
            x0 = x[i]
            ax.plot([x0, x[c]], [y[c], y[c]], color="#4A5560", lw=1.4, zorder=2)
            k = len(n.mutations)
            for m, (pos, alt) in enumerate(n.mutations):
                xm = x0 + (m + 0.5) * (x[c] - x0) / k
                col = CLASS[mclass(ref[pos - 1], alt)][1]
                ax.plot([xm, xm], [y[c] - 0.16, y[c] + 0.16], color=col, lw=2.2, zorder=3)
                ax.text(xm, y[c] - 0.2, f"{ref[pos - 1]}{pos}{alt}", rotation=60, fontsize=5.6, color=col,
                        ha="left", va="bottom", rotation_mode="anchor", zorder=5)
    ax.scatter([0], [y[0]], s=90, marker="s", color="#2E3A46", zorder=4)
    ax.text(-0.02 * span, y[0] + 0.45, f"ancestor\n{tree['reads_at'].get(0, 0)} reads", ha="right", va="top",
            fontsize=6.5, color="#2E3A46")
    for i, n in nodes.items():
        share = len(n.reads) / total
        ax.scatter([x[i]], [y[i]], s=25 + 700 * share, color="white", edgecolor="#2E3A46", lw=1.1, zorder=4)
        ax.text(x[i] + 0.012 * span + 0.05 * (share ** 0.5) * span, y[i], f"{len(n.reads)}", fontsize=6.3,
                color="#2E3A46", ha="left", va="center", zorder=5)
    ax.set_ylim(max(y.values()) + 0.8, min(y.values()) - 1.6)
    ax.set_xlim(-0.12 * span, span * 1.12)
    ax.set_xlabel("Mutations accumulated since the clone's unmutated ancestor")
    ax.set_yticks([])
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    ax.set_title(title, fontsize=9.5, loc="left", fontweight="bold")


def main():
    proj = {"name": "JH4 mutation trees", "region": "jh4_intron",
            "reference": str(REF / "mouse_jh4_intron.fasta"),
            "contaminants": [str(REF / "mouse_jh4_contaminant.fasta")], "mode": "recommended", "depth": 500,
            "overrides": {"analysis_window": [90, 1230]},
            "samples": [{"name": n, "path": str(DEP / f"{f}_C57BL6_sub500.fastq"),
                         "group": "Wt" if n.startswith("WT") else "Aicda-/-", "control": n.startswith("KO")}
                        for n, f in LIBS.items()]}
    res = project.run_project(proj)
    p = project.build_params(proj)
    ref = res["reference"]
    OUT.mkdir(parents=True, exist_ok=True)
    candidates, summary = [], {}
    for m in ("WT2", "WT3", "WT4"):
        d = res["_details"][m]
        per_clone = mutation_tree.clone_reads(d["_aligned"], d["per_read"], d["variants"], p)
        rows = []
        for cid, reads in per_clone.items():
            if len(reads) < MIN_CLONE_READS:
                continue
            t = mutation_tree.build(reads)
            depth_max = max((mutation_tree.depth(t, n.id) for n in t["nodes"]), default=0)
            rows.append({"clone": cid, "reads": len(reads), "nodes": len(t["nodes"]), "max_depth": depth_max,
                         "branch_points": sum(1 for n in t["nodes"] if len(n.children) > 1)
                         + (sum(1 for n in t["nodes"] if n.parent is None) > 1),
                         "fit": t["fit"], "left_out_single_read": t["left_out"]})
            # informative: branches, and more than one mutation deep
            candidates.append((rows[-1]["branch_points"], len(t["nodes"]), m, cid, t))
        summary[m] = rows
    candidates.sort(key=lambda c: (-c[0], -c[1]))
    # informative trees: branched, at least three mutations deep, and most reads consistent
    chosen = [c for c in candidates if c[4]["nodes"] and c[0] and (c[4]["fit"] or 0) >= 0.8
              and max(mutation_tree.depth(c[4], n.id) for n in c[4]["nodes"]) >= 3][:3]
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1], hspace=0.45, wspace=0.18)
    axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    picked = {}
    for ax, (bp, nn, m, cid, t), letter in zip(axes, chosen, "ABC", strict=False):
        draw(ax, t, ref, f"{letter}   Wt{m[2]} clone {cid} · {t['reads']} reads · "
                         f"{round(100 * t['fit'])}% of mutated reads fit")
        picked[f"{m}:{cid}"] = t
    handles = [Line2D([], [], color=c, marker="s", ls="", label=lab) for lab, c in CLASS.values()]
    handles.append(Line2D([], [], marker="o", ls="", mfc="white", mec="#2E3A46", ms=9,
                          label="node: mutations on one branch; area, reads carrying it"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.subplots_adjust(bottom=0.1, top=0.95, left=0.04, right=0.98)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"mutation_trees.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    json.dump({"trees": {k: [{"id": n.id, "parent": n.parent, "reads": len(n.reads),
                              "mutations": [f"{ref[q - 1]}{q}{a}" for q, a in n.mutations]} for n in t["nodes"]]
                         for k, t in picked.items()},
               "clones": summary}, open(OUT / "trees.json", "w"), indent=1)
    for m, rows in summary.items():
        print(m, rows[:5])


if __name__ == "__main__":
    pipeline  # imported for its side effects on presets
    main()
