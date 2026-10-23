"""Figure S17, page 2: the read-level reconstruction of Figure 6D, repeated in every Wt mouse.

Replaces mutation_trees_by_mouse.png from docs/calibration/mutation_tree_figure.py, which was drawn 13.5 in wide and
printed at half size (labels near 3.7 pt). Same trees, same choice of clone, drawn at page width in Figure 6D's
style, with each mutation named.

    ../../../../.venv/bin/python figS17_trees.py --compute   # rebuilds the trees from the delivered reads (SHMtrace)
    python3 figS17_trees.py                                  # draws from the saved trees

The clone per mouse is the one mutation_tree_figure.py picks: most branch points, then deepest branch. The saved
WT2 and WT3 trees are checked node for node against figures/mutation_trees/trees.json.
"""
import json
import sys

from matplotlib import pyplot as plt

from fig6_longreads import CLASS, draw_tree, klass
from style import DATA, FULL, INK, MID, MM, SCRATCH, save, setup

CAL = DATA.parents[1]
TREES = SCRATCH / "s17_trees.json"
MICE = ("WT2", "WT3", "WT4")


def compute():
    sys.path.insert(0, str(CAL))
    import mutation_tree_figure as mtf
    from shmtrace import mutation_tree, project
    proj = {"name": "JH4 mutation trees", "region": "jh4_intron",
            "reference": str(mtf.REF / "mouse_jh4_intron.fasta"),
            "contaminants": [str(mtf.REF / "mouse_jh4_contaminant.fasta")], "mode": "recommended", "depth": 500,
            "overrides": {"analysis_window": [90, 1230]},
            "samples": [{"name": n, "path": str(mtf.DEP / f"{f}.fastq"),
                         "group": "Wt" if n.startswith("WT") else "Aicda-/-", "control": n.startswith("KO")}
                        for n, f in mtf.LIBS.items()]}
    res = project.run_project(proj)
    p = project.build_params(proj)
    ref = res["reference"]
    out = {}
    for m in MICE:
        d = res["_details"][m]
        best = None
        for cid, reads in mutation_tree.clone_reads(d["_aligned"], d["per_read"], d["variants"], p).items():
            if len(reads) < mtf.MIN_CLONE_READS:
                continue
            t = mutation_tree.build(reads)
            if not t["nodes"]:
                continue
            bp = (sum(1 for n in t["nodes"] if len(n.children) > 1)
                  + (sum(1 for n in t["nodes"] if n.parent is None) > 1))
            key = (bp, max(mutation_tree.depth(t, n.id) for n in t["nodes"]))
            if best is None or key > best[0]:
                best = (key, cid, t)
        (bp, deep), cid, t = best
        out[m] = {"clone": cid, "reads": t["reads"], "fit": t["fit"], "depth": deep,
                  "ancestor_reads": t["reads_at"].get(0, 0),
                  "nodes": [{"id": n.id, "parent": n.parent, "reads": len(n.reads),
                             "mutations": [f"{ref[q - 1]}{q}{a}" for q, a in n.mutations]} for n in t["nodes"]]}
    saved = json.load(open(CAL / "figures/mutation_trees/trees.json"))
    assert out["WT2"]["nodes"] == saved["trees"]["WT2:0"], "WT2 tree differs from trees.json"
    assert out["WT3"]["nodes"] == saved["wt3_figure_clone"]["nodes"], "WT3 tree differs from Figure 6D"
    TREES.write_text(json.dumps(out, indent=1))
    print({m: (v["clone"], v["reads"], len(v["nodes"]), v["depth"], v["fit"]) for m, v in out.items()})


def name_mutations(ax, nodes, y_of, x_end):
    """Each mutation's name under its tick, in its class colour, at draw_tree's positions (the read counts sit
    above the line, so the two never meet)."""
    for t in list(ax.texts):                       # draw_tree's read counts, centred over each node, are redrawn
        t.remove()                                 # above-right of it, clear of the names on the row above
    for n in nodes:
        r = (6 + 1.5 * n["reads"]) ** 0.5 / 2      # the node's radius in points (draw_tree sizes it so)
        ax.annotate(str(n["reads"]), (x_end[n["id"]], y_of[n["id"]]), xytext=(0.6 * r + 0.5, 0.6 * r + 0.5),
                    textcoords="offset points", fontsize=6, color=MID, ha="left", va="bottom")
        x0 = x_end[n["id"]] - len(n["mutations"])
        for k, m in enumerate(n["mutations"]):
            ax.text(x0 + k + 0.5, y_of[n["id"]] + 0.27, m, fontsize=6, color=CLASS[klass(m)][1],
                    ha="center", va="top")


def main():
    setup()
    trees = json.loads(TREES.read_text())
    leaves = [max(1, sum(1 for n in trees[m]["nodes"] if not any(c["parent"] == n["id"] for c in trees[m]["nodes"])))
              for m in MICE]
    fig = plt.figure(figsize=(FULL, 222 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[x + 1.2 for x in leaves], left=0.03, right=0.985, top=0.955,
                          bottom=0.1, hspace=0.42)
    for i, m in enumerate(MICE):
        t = trees[m]
        ax = fig.add_subplot(gs[i])
        name_mutations(ax, t["nodes"], *draw_tree(ax, t["nodes"]))
        ax.set_xlim(-3.2 if t["ancestor_reads"] else -1.4, t["depth"] + 0.8)
        ax.set_xticks(range(0, t["depth"] + 1, 2))
        if t["ancestor_reads"]:                    # reads that stay at the ancestor, beside its square
            ax.text(-1.1, 0, f"{t['ancestor_reads']} reads\nat the ancestor", fontsize=6, color=MID, ha="right",
                    va="center")
        ax.set_title(f"Wt{m[2]}: its most branched clone, {t['reads']} reads; deepest branch {t['depth']} mutations; "
                     f"{round(100 * t['fit'])}% of mutated reads fit one path", fontsize=7, pad=6)
        fig.text(0.008, ax.get_position().y1 + 0.012, "ABC"[i], fontsize=8, fontweight="bold", va="bottom")
    hs = [plt.Line2D([], [], color=c, lw=1.5, label=f"mutation: {lab}") for lab, c in CLASS.values()]
    hs += [plt.Line2D([], [], marker="s", ls="", ms=4, color=INK, label="inferred unmutated ancestor"),
           plt.Line2D([], [], marker="o", ls="", ms=5, mfc="white", mec=INK,
                      label="node: number = reads carrying every mutation up to it"),
           plt.Line2D([], [], color="#8a96a3", lw=0.8, label="branch: reads diverge")]
    fig.legend(handles=hs, loc="lower center", ncol=3, fontsize=6, handlelength=1.2, columnspacing=1.2,
               bbox_to_anchor=(0.5, 0.005))
    save(fig, "FigureS17_mutation_trees")


if __name__ == "__main__":
    compute() if "--compute" in sys.argv else main()
