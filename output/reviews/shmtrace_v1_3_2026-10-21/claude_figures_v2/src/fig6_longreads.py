"""Figure 6: what long reads add.

A  unique variants along the whole J_H4 amplicon (positions 90-1,230), per Wt mouse and the AID-null controls
B  AID targeting: calls per position by 5-mer mutability quintile (SHazaM HH_S5F), Wt calls against
   every Aicda-/- mismatch (sequencing error)
C  clones from the V(D)J junction: each mouse's 500 reads by clone, the reads whose junction could not be read,
   and the read- against clone-weighted rate
(Claude revision: A direct labels and the nanopore-only interval; B defines the index; C splits clone composition
from the read- vs clone-weighted comparison; D carries a reading guide.)
D  a read-level reconstruction: the reads of one Wt3 clone ordered into a mutation tree, with the share of
   reads that fit no single path
"""
import csv
import json
from collections import defaultdict

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt

from style import CASE, CTRL, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]
FULLRUN = CAL / "figures/jh4_shmtrace_full"
CLASS = {"GC_ts": ("G:C transition", "#c0392b"), "GC_tv": ("G:C transversion", "#e08a2c"),
         "AT": ("A:T", "#2c5d8a")}


def klass(m):
    f, t = m[0], m[-1]
    if f in "GC":
        return "GC_ts" if f + t in ("GA", "CT") else "GC_tv"
    return "AT"


def draw_tree(ax, nodes):
    kids = defaultdict(list)
    for n in nodes:
        kids[n["parent"]].append(n)
    for v in kids.values():
        v.sort(key=lambda n: -n["reads"])
    y_of, x_end = {}, {}
    leaf = [0]

    def place(n, x0):
        x1 = x0 + len(n["mutations"])
        x_end[n["id"]] = x1
        ch = kids.get(n["id"], [])
        if not ch:
            y_of[n["id"]] = leaf[0]
            leaf[0] += 1
        else:
            ys = [place(c, x1) for c in ch]
            y_of[n["id"]] = ys[0] if len(ch) == 1 else ys[0]
        return y_of[n["id"]]

    roots = kids[None]
    # every root hangs off the clone's unmutated ancestor at x = 0
    for r in roots:
        place(r, 0)
    ax.plot([-0.8], [y_of[roots[0]["id"]]], marker="s", ms=4, color=INK)
    for n in nodes:
        x0 = x_end[n["id"]] - len(n["mutations"])
        y = y_of[n["id"]]
        py = y_of[n["parent"]] if n["parent"] is not None else y_of[roots[0]["id"]]
        ax.plot([x0, x0], [py, y], color=SOFT, lw=0.6)
        ax.plot([x0 if n["parent"] is not None else -0.8, x_end[n["id"]]], [y, y], color=INK, lw=0.8)
        for k, m in enumerate(n["mutations"]):
            ax.plot([x0 + k + 0.5] * 2, [y - 0.22, y + 0.22], color=CLASS[klass(m)][1], lw=1.4)
        ax.scatter(x_end[n["id"]], y, s=6 + 1.5 * n["reads"], facecolor="white", edgecolor=INK, lw=0.6, zorder=3)
        ax.text(x_end[n["id"]], y - 0.42, str(n["reads"]), fontsize=6, color=MID, ha="center")
    ax.set_ylim(leaf[0] - 0.4, -0.95)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Mutations since the clone's unmutated ancestor")
    return y_of, x_end


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 134 * MM))
    gs = gridspec.GridSpec(2, 3, figure=fig, left=0.08, right=0.975, top=0.945, bottom=0.2, wspace=0.42,
                           hspace=0.62, width_ratios=[1, 1, 1])
    bottom = gs[1, :].subgridspec(1, 3, wspace=0.35, width_ratios=[0.95, 0.42, 1.9])
    # A: full-length map
    ax = fig.add_subplot(gs[0, :2])
    proj = json.load(open(FULLRUN / "project.json"))
    depth = defaultdict(lambda: defaultdict(int))
    for r in csv.DictReader(open(FULLRUN / "coverage.csv")):
        depth[r["sample"]][int(r["pos"])] = int(r["depth"])
    edges = np.arange(90, 1231 + 50, 50)
    mids = (edges[:-1] + edges[1:]) / 2
    curves = {}
    for s in ("WT2", "WT3", "WT4", "KO1", "KO2", "KO3"):
        pos = [v["pos"] for v in proj["variants"] if v["sample"] == s]
        h, _ = np.histogram(pos, edges)
        bases = [sum(depth[s][p] for p in range(a, min(b, 1231))) for a, b in zip(edges[:-1], edges[1:])]
        curves[s] = 1e5 * h / np.maximum(bases, 1)
    ax.axvspan(101, 629, color="#f0e2cf", alpha=0.7, lw=0)
    ax.text(365, 0.97, "shared Sanger window (529 bp)", transform=ax.get_xaxis_transform(), ha="center",
            fontsize=6, color="#a3621b", va="top")
    ax.annotate("", xy=(1230, 0.9), xytext=(629, 0.9), xycoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=MID, shrinkA=0, shrinkB=0))
    ax.text(930, 0.93, "additional interval read only by nanopore", transform=ax.get_xaxis_transform(),
            ha="center", fontsize=6, color=MID, va="bottom")
    for s in ("WT2", "WT3", "WT4"):
        ax.plot(mids, curves[s], color="#9fc0d8", lw=0.5, alpha=0.9)
    wt_mean = np.mean([curves[s] for s in ("WT2", "WT3", "WT4")], axis=0)
    ko_mean = np.mean([curves[s] for s in ("KO1", "KO2", "KO3")], axis=0)
    ax.plot(mids, wt_mean, color=CASE, lw=1.5)
    ax.plot(mids, ko_mean, color=CTRL, lw=1.2)
    ax.text(mids[3], wt_mean[3] + 38, "Wt mean", fontsize=6, color=CASE, fontweight="bold", ha="center")
    ax.text(1000, 58, "individual Wt mice", fontsize=6, color="#5d86b3", ha="center", va="bottom")
    ax.text(400, 5, "$\\it{Aicda}^{{-}/{-}}$ mean", fontsize=6, color=CTRL, ha="center", va="bottom")
    ax.set_xlim(90, 1230)
    ax.set_xlabel("J$_H$4 amplicon position (bp); 50-bp bins")
    ax.set_ylabel("Unique variants per 10$^5$ bases")
    ax.set_ylim(-8, 235)
    ax.set_title("The whole amplicon from each read (positions 90–1,230)", fontsize=7)
    letter(ax, "A", dx=-0.08)
    # B: targeting dose-response
    ax = fig.add_subplot(gs[0, 2])
    t = json.load(open(DATA / "targeting_jh4.json"))
    q = t["quintiles"]
    wt = np.array([x["wt_calls_per_position"] for x in q])
    ko = np.array([x["ko_mismatch_events_per_position"] for x in q])
    ax.plot(range(1, 6), wt / wt.mean(), color=CASE, marker="o", ms=3, lw=1.2, label="Wt calls")
    ax.plot(range(1, 6), ko / ko.mean(), color=CTRL, marker="s", ms=3, lw=1.2, label="$\\it{Aicda}^{{-}/{-}}$ mismatches")
    ax.axhline(1, color=SOFT, lw=0.5, ls=":")
    ax.set_xticks(range(1, 6), ["1\nleast", "2", "3", "4", "5\nmost"])
    npos = [x["positions"] for x in q]
    ax.set_xlabel(f"S5F mutability quintile\n(~{round(np.mean(npos))} positions each)")
    ax.set_ylabel("Per position, relative to mean")
    wi = [t["wt"][s]["index"] for s in ("WT2", "WT3", "WT4")]
    ax.text(0.04, 0.97, f"targeting index (mean mutability at\ncalled positions ÷ window mean):\n"
            f"Wt mice {min(wi):.2f}–{max(wi):.2f}; null {t['ko_all_positions']['index']:.2f}",
            transform=ax.transAxes, va="top", fontsize=6, color=MID)
    ax.set_ylim(0.3, 2.05)
    ax.set_title("AID targeting against error", fontsize=7)
    ax.legend(fontsize=6, loc="lower right")
    letter(ax, "B", dx=-0.3)
    # C: clone composition only; the read- against clone-weighted comparison has its own small plot beside it
    ax = fig.add_subplot(bottom[0])
    ex = {s["sample"]: s for s in json.load(open(DATA / "extras_jh4/project.json"))["samples"]}
    clones = defaultdict(list)
    for r in csv.DictReader(open(CAL / "figures/jh4_shmtrace/clones.csv")):
        clones[r["sample"]].append(int(r["reads"]))
    for i, s in enumerate(("WT2", "WT3", "WT4")):
        sizes = sorted(clones[s], reverse=True)
        e = ex[s]
        tot = e["reads"]                     # every analyzed read, so the unreadable junctions show
        left = 0
        for k, z in enumerate(sizes):
            ax.barh(i, 100 * z / tot, left=left, height=0.6, color=CASE if k == 0 else "#9fc0d8",
                    edgecolor="white", lw=0.5)
            left += 100 * z / tot
        ax.barh(i, 100 - left, left=left, height=0.6, color="white", edgecolor=SOFT, lw=0.5, hatch="////")
        ax.text(1.5, i - 0.36, f"{len(sizes)} clones; junction read in {100 * e['reads_with_junction'] / tot:.0f}%",
                fontsize=6, color=MID, ha="left", va="bottom")
    ax.set_yticks(range(3), ["Wt2", "Wt3", "Wt4"])
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel(f"Reads (% of {tot})")
    ax.set_ylim(2.45, -0.75)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=CASE, label="largest clone"), Patch(color="#9fc0d8", label="other clones"),
                       Patch(facecolor="white", edgecolor=SOFT, hatch="////", label="junction unreadable")],
              fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3, handlelength=1,
              columnspacing=0.6, handletextpad=0.3)
    ax.set_title("Clones from the V(D)J junction", fontsize=7)
    letter(ax, "C", dx=-0.3)
    # C, right: what the clone structure does to the frequency, per mouse
    ax = fig.add_subplot(bottom[1])
    for i, s in enumerate(("WT2", "WT3", "WT4")):
        e = ex[s]
        a, b = e["mutation_frequency_x1e3"], e["clone_weighted_frequency_x1e3"]
        ax.annotate("", xy=(b, i), xytext=(a, i), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=MID,
                                                                  mutation_scale=6, shrinkA=2, shrinkB=2))
        ax.scatter([a], [i], s=14, facecolor="white", edgecolor=CASE, linewidth=0.8, zorder=3)
        ax.scatter([b], [i], s=14, color=CASE, zorder=3)
    ax.set_yticks(range(3), [])
    ax.set_ylim(2.45, -0.75)
    ax.set_xlim(0, 28)
    ax.set_xticks([0, 10, 20])
    ax.set_xlabel("per 10$^3$ bp")
    ax.scatter([], [], s=14, facecolor="white", edgecolor=CASE, label="read-weighted")
    ax.scatter([], [], s=14, color=CASE, label="clone-weighted")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.45, -0.32), handletextpad=0.2, borderaxespad=0)
    ax.set_title("Frequency", fontsize=7)
    # D: mutation tree
    ax = fig.add_subplot(bottom[2])
    tr = json.load(open(CAL / "figures/mutation_trees/trees.json"))["wt3_figure_clone"]
    draw_tree(ax, tr["nodes"])
    # every read of this clone carries the founding node, so fit is a share of all its reads
    carriers = max(n["reads"] for n in tr["nodes"])
    n_fit = round(tr["fit"] * carriers)
    ax.set_title(f"Read-level reconstruction of one Wt3 clone ({tr['reads']} reads)", fontsize=7, pad=19)
    ax.text(0.5, 1.01, f"{n_fit} of {carriers} reads ({round(100 * tr['fit'])}%) fit one path; the other {carriers - n_fit} carry\n"
            f"mutations of more than one branch and are counted at each node they carry", transform=ax.transAxes, fontsize=6, color=MID, ha="center", va="bottom")
    # reading guide
    hs = [plt.Line2D([], [], color=c, lw=1.5, label=f"mutation: {lab}") for lab, c in CLASS.values()]
    hs += [plt.Line2D([], [], marker="s", ls="", ms=4, color=INK, label="inferred unmutated ancestor"),
           plt.Line2D([], [], marker="o", ls="", ms=5, mfc="white", mec=INK,
                      label="node: number = reads carrying every\nmutation up to it (size scales with it)"),
           plt.Line2D([], [], color=SOFT, lw=0.8, label="branch: reads diverge")]
    ax.legend(handles=hs, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, handlelength=1.2,
              columnspacing=1.0, borderaxespad=0.1)
    ax.xaxis.set_major_locator(plt.MultipleLocator(5))
    letter(ax, "D", dx=-0.06, dy=1.2)
    save(fig, "Figure6")


if __name__ == "__main__":
    main()
