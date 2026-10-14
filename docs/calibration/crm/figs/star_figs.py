"""STAR Protocols Figures 2-4.

Figure 2  the SHMtrace application on the J_H4 project: set up, results, tables, power (screenshots)
Figure 3  expected outcomes: densities per mouse, read accounting, spectrum, full-length map
Figure 4  troubleshooting: a dominant clone, an error cluster, depth
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from PIL import Image

from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, setup

CRM = DATA.parent
CAL = CRM.parent
STAR = CRM / "star"
SCREENS = STAR / "screens"
OUT = STAR / "figures"
TYPES = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
        fig.savefig(OUT / f"{name}.{ext}", facecolor="white", **kw)
    Image.open(OUT / f"{name}.png").convert("RGB").save(OUT / f"{name}.tif", dpi=(600, 600), compression="tiff_lzw")
    print(name)


def fig2():
    setup()
    fig = plt.figure(figsize=(FULL, 118 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.03, right=0.99, top=0.96, bottom=0.02, wspace=0.08, hspace=0.12)
    shots = [("app_setup.png", "A"), ("app_results.png", "B"), ("app_tables.png", "C"), ("app_power.jpg", "D")]
    for k, (f, lab) in enumerate(shots):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        im = Image.open(SCREENS / f)
        ax.imshow(im, interpolation="lanczos")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_color("#b4bec8")
            s.set_linewidth(0.6)
        ax.text(-0.03, 1.0, lab, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top", ha="right")
    save(fig, "STAR_Figure2")


def fig3():
    setup()
    proj = json.load(open(STAR / "cli_run/project.json"))
    fig = plt.figure(figsize=(FULL, 104 * MM))
    gs = gridspec.GridSpec(2, 3, figure=fig, left=0.09, right=0.99, top=0.93, bottom=0.11, wspace=0.5, hspace=0.7,
                           width_ratios=[0.8, 1.2, 1.2])
    # A densities
    ax = fig.add_subplot(gs[0, 0])
    s = {x["sample"]: x for x in proj["samples"]}
    for i, (grp, f, e) in enumerate((("WT", CASEF, CASE), ("Aicda-/-", CTRLF, CTRL))):
        rows = [x for x in s.values() if x["group"] == grp]
        vals = [x["unique_density_x1e5"] or x["one_event_density_x1e5"] for x in rows]
        for j, (x, v) in enumerate(zip(rows, vals)):
            ax.scatter(i + (j - 1) * 0.12, v, s=14, marker="o" if x["unique_density_x1e5"] else "v",
                       facecolor=f, edgecolor=e, lw=0.6, zorder=3)
        g = 10 ** np.mean(np.log10(vals))
        ax.plot([i - 0.25, i + 0.25], [g, g], color=e, lw=1.4)
    from scipy import stats
    gv = {grp: [x["unique_density_x1e5"] or x["one_event_density_x1e5"] for x in s.values() if x["group"] == grp]
          for grp in ("WT", "Aicda-/-")}
    fold = 10 ** (np.mean(np.log10(gv["WT"])) - np.mean(np.log10(gv["Aicda-/-"])))
    pw = stats.ttest_ind(np.log10(gv["WT"]), np.log10(gv["Aicda-/-"]), equal_var=False).pvalue
    ax.text(0.5, 400, f"≥{fold:.0f}×, p = {pw:.4f}", ha="center", fontsize=5.8)
    ax.set_yscale("log")
    ax.set_ylim(0.2, 1000)
    ax.set_xticks([0, 1], ["Wt", "Aicda$^{-/-}$"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("Unique variants per 10$^5$ bases")
    ax.set_title("Wt vs AID-deficient", fontsize=7)
    letter(ax, "A", dx=-0.45)
    # B read accounting
    ax = fig.add_subplot(gs[0, 1:])
    acct = json.load(open(CAL / "figures/jh4_shmtrace/table_s1.json"))
    parts = [("too_long", "over-length", "#5f6b75"), ("low_quality", "expected error > 0.01", "#9aa5b1"),
             ("not_target", "unmapped or paralog", "#d9952b"), ("chimeric", "chimeric", "#c0392b"),
             ("short_alignment", "short alignment", "#e6b98a"), ("paralogue_alleles", "paralog alleles", "#6a4c93"),
             ("usable", "usable", "#1f6f8b")]
    names = {"WT2": "Wt2", "WT3": "Wt3", "WT4": "Wt4", "KO1": "KO1", "KO2": "KO2", "KO3": "KO3"}
    for i, (lib, a) in enumerate(acct.items()):
        left = 0
        for k, lab, col in parts:
            ax.barh(i, a[k], left=left, color=col, height=0.7, label=lab if i == 0 else None)
            left += a[k]
        ax.text(left + 60, i, f"{a['usable']:,} usable", va="center", fontsize=5.5, color=MID)
    ax.axvline(500, color=INK, lw=0.6, ls=":")
    ax.set_yticks(range(len(acct)), [names[x] for x in acct])
    ax.invert_yaxis()
    ax.set_xlim(0, 6400)
    ax.set_xlabel("Delivered reads")
    ax.set_title("Read accounting per library (dotted line: 500 analyzed)", fontsize=7)
    ax.legend(fontsize=5.2, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.28), handlelength=1, columnspacing=0.8)
    letter(ax, "B", dx=-0.1)
    # C spectrum
    ax = fig.add_subplot(gs[1, :2])
    sp = proj["spectrum_group"]["WT"] if "spectrum_group" in proj else None
    counts = {t_: 0 for t_ in TYPES}
    for v in proj["variants"]:
        if v["sample"].startswith("Wt"):
            counts[f"{v['ref']}>{v['alt']}"] += 1
    tot = sum(counts.values())
    cols = ["#d9952b" if t_[0] in "AT" else "#2c5d8a" for t_ in TYPES]
    ax.bar(range(12), [100 * counts[t_] / tot for t_ in TYPES], color=cols)
    ax.set_xticks(range(12), [t_.replace(">", "→") for t_ in TYPES], rotation=0, fontsize=5.5)
    ax.set_ylabel("Share of Wt calls (%)")
    ax.set_title(f"Spectrum, three Wt mice pooled (n = {tot} unique variants)", fontsize=7)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#2c5d8a", label="from G:C"), Patch(color="#d9952b", label="from A:T")],
              fontsize=5.5, loc="upper right")
    letter(ax, "C", dx=-0.09)
    # D full length
    ax = fig.add_subplot(gs[1, 2])
    full = CAL / "figures/jh4_shmtrace_full"
    p = json.load(open(full / "project.json"))
    depth = defaultdict(lambda: defaultdict(int))
    for r in csv.DictReader(open(full / "coverage.csv")):
        depth[r["sample"]][int(r["pos"])] = int(r["depth"])
    edges = np.arange(90, 1231 + 50, 50)
    mids = (edges[:-1] + edges[1:]) / 2
    curves = {}
    for smp in ("WT2", "WT3", "WT4", "KO1", "KO2", "KO3"):
        h, _ = np.histogram([v["pos"] for v in p["variants"] if v["sample"] == smp], edges)
        b = [sum(depth[smp][q] for q in range(a, min(c, 1231))) for a, c in zip(edges[:-1], edges[1:])]
        curves[smp] = 1e5 * h / np.maximum(b, 1)
    for smp in ("WT2", "WT3", "WT4"):
        ax.plot(mids, curves[smp], color=CASEF, lw=0.5)
    ax.plot(mids, np.mean([curves[x] for x in ("WT2", "WT3", "WT4")], axis=0), color=CASE, lw=1.4, label="Wt")
    ax.plot(mids, np.mean([curves[x] for x in ("KO1", "KO2", "KO3")], axis=0), color=CTRL, lw=1.1, label="Aicda$^{-/-}$")
    ax.set_xlabel("J$_H$4 position (bp)")
    ax.set_ylabel("Calls per 10$^5$ bases")
    ax.set_title("Whole amplicon", fontsize=7)
    ax.legend(fontsize=5.5)
    letter(ax, "D", dx=-0.32)
    save(fig, "STAR_Figure3")


def fig4():
    setup()
    fig = plt.figure(figsize=(FULL, 60 * MM))
    gs = gridspec.GridSpec(1, 3, figure=fig, left=0.08, right=0.99, top=0.86, bottom=0.2, wspace=0.42)
    # A dominant clone
    ax = fig.add_subplot(gs[0, 0])
    ex = {s["sample"]: s for s in json.load(open(DATA / "extras_jh4/project.json"))["samples"]}
    clones = defaultdict(list)
    for r in csv.DictReader(open(CAL / "figures/jh4_shmtrace/clones.csv")):
        clones[r["sample"]].append(int(r["reads"]))
    for i, s in enumerate(("WT2", "WT3", "WT4")):
        sizes = sorted(clones[s], reverse=True)
        tot, left = sum(sizes), 0
        for k, z in enumerate(sizes):
            ax.barh(i, 100 * z / tot, left=left, height=0.6, color=CASE if k == 0 else "#9fc0d8", edgecolor="white", lw=0.3)
            left += 100 * z / tot
        e = ex[s]
        ax.text(102, i, f"{e['mutation_frequency_x1e3']:.1f} → {e['clone_weighted_frequency_x1e3']:.1f}",
                va="center", fontsize=5.5, color=MID)
    ax.set_ylim(2.5, -1.0)
    ax.text(102, -0.72, "per 10$^3$ bp, read-\n→ clone-weighted", fontsize=5, color=SOFT, va="center")
    ax.set_yticks(range(3), ["Wt2", "Wt3", "Wt4"])
    ax.set_xlim(0, 150)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Reads by clone (%)")
    ax.set_title("A dominant clone", fontsize=7)
    letter(ax, "A", dx=-0.3)
    # B error cluster
    ax = fig.add_subplot(gs[0, 1])
    p = json.load(open(DATA / "public/ch12/project.json"))
    for g, col, lab in (("AID", CASE, "AID"), ("R174ER178E", "#a3621b", "R174E/R178E")):
        pos = [v["pos"] for v in p["variants"] if v["group"] == g]
        ax.hist(pos, bins=np.arange(21, 251, 5), color=col, alpha=0.85, label=f"{lab} (n = {len(pos)})", histtype="stepfilled")
    ax.axvspan(164, 172, color="#f3c78e", alpha=0.5, lw=0)
    ax.set_xlabel("Amplicon position (bp)")
    ax.set_ylabel("Calls per 5 bp")
    ax.set_ylim(0, 17)
    ax.set_title("An error cluster (CH12F3)", fontsize=7)
    ax.legend(fontsize=5.3, loc="upper left")
    letter(ax, "B", dx=-0.3)
    # C depth
    ax = fig.add_subplot(gs[0, 2])
    d = json.load(open(CAL / "figures/depth_shmtrace/jh4.json"))
    depths = sorted(d, key=int)
    for grp, col, lab in ((("WT2", "WT3", "WT4"), CASE, "Wt"), (("KO1", "KO2", "KO3"), CTRL, "Aicda$^{-/-}$")):
        m = [10 ** np.mean(np.log10([max(d[x][s]["calls"], 1) for s in grp])) for x in depths]
        ax.plot([int(x) for x in depths], m, marker="o", ms=3, color=col, lw=1.1, label=lab)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([30, 125, 500, 1000], ["30", "125", "500", "1,000"])
    ax.minorticks_off()
    ax.axvline(500, color=SOFT, lw=0.6, ls=":")
    ax.set_xlabel("Reads analyzed")
    ax.set_ylabel("Unique variants (geometric mean)")
    ax.set_title("Depth", fontsize=7)
    ax.legend(fontsize=5.5, loc="center right")
    letter(ax, "C", dx=-0.3)
    save(fig, "STAR_Figure4")


if __name__ == "__main__":
    fig2()
    fig3()
    fig4()
