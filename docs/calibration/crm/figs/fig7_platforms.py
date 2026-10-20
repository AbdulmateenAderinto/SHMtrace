"""Figure 7: the same engine across platforms, loci and questions.

A  Illumina J_H4 (PRJNA973112): mutation frequency per mouse, Wt / Ldha-deleted / AID-null
B  CH12F3 Ung-/-Msh2-/-Aicda-/- (PRJNA918596): mutation frequency per replicate by condition
C  spectrum and targeting of AID and R174E/R178E calls
D  Ramos (PRJNA820976): deletions and insertions per clone
E  UMI consensus on R10 reads (Karst et al.): error per base by reads per molecule
F  ectopic AID in non-B cells (PRJNA170048): transitions among C:G mutations with and without the
   catalytically dead AID control, in three repair contexts
"""
import json
from collections import defaultdict

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt

from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]
PUB = DATA / "public"


def titled(ax, title, question, size=7):
    """Panel title, with the validation question it answers on a smaller line beneath."""
    ax.set_title(title, fontsize=size, pad=10)
    ax.text(0.5, 1.012, question, transform=ax.transAxes, ha="center", va="bottom", fontsize=6, color=MID,
            style="italic")


def strip(ax, groups, values, colors, floor):
    for i, (g, (f, e)) in enumerate(zip(groups, colors)):
        v = np.array(values[g], dtype=float)
        shown = np.where(v > 0, v, floor)
        x = i + np.linspace(-0.15, 0.15, len(v))
        for xx, yy, raw in zip(x, shown, v):
            ax.scatter(xx, yy, s=12, facecolor=f, edgecolor=e, lw=0.6, zorder=3, marker="o" if raw > 0 else "v")
        pos = v[v > 0]
        if len(pos):
            g_ = 10 ** np.mean(np.log10(np.where(v > 0, v, floor)))
            ax.plot([i - 0.28, i + 0.28], [g_, g_], color=e, lw=1.4)
    ax.axhline(floor * 1.5, color=SOFT, lw=0.5, ls=":")


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 112 * MM))
    gs = gridspec.GridSpec(2, 6, figure=fig, left=0.075, right=0.975, top=0.93, bottom=0.17, wspace=1.35, hspace=0.95)
    # A
    ax = fig.add_subplot(gs[0, 0:2])
    p = json.load(open(PUB / "ldha/project.json"))
    groups = ["WT", "LDHA_KO", "AID_KO"]
    vals = {g: [s["mutation_frequency_x1e3"] for s in p["samples"] if s["group"] == g] for g in groups}
    strip(ax, groups, vals, [(CASEF, CASE), ("#b9a3d6", "#6a4c93"), (CTRLF, CTRL)], 0.01)
    ax.set_yscale("log")
    ax.set_xticks(range(3), ["Wt", "$Ldha$ deleted", "AID-null"])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    titled(ax, "J$_H$4, Illumina (Sharma et al.)", "a genotype effect across mice, another platform")
    letter(ax, "A", dx=-0.32)
    # B
    ax = fig.add_subplot(gs[0, 2:4])
    p = json.load(open(PUB / "ch12/project.json"))
    groups = ["AID", "R50E", "R174ER178E", "EV"]
    vals = {g: [s["mutation_frequency_x1e3"] for s in p["samples"] if s["group"] == g] for g in groups}
    strip(ax, groups, vals, [(CASEF, CASE), ("#d9dde2", MID), ("#f3c78e", "#a3621b"), (CTRLF, CTRL)], 0.005)
    ax.set_yscale("log")
    ax.set_xticks(range(4), ["AID", "R50E", "R174E/\nR178E", "empty"])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    titled(ax, "CH12F3, Illumina (Wang et al.)", "activity of AID mutants in one cell line")
    letter(ax, "B", dx=-0.32)
    # C
    ax = fig.add_subplot(gs[0, 4:6])
    s5 = json.load(open(DATA / "targeting_ch12.json"))
    cats = [("C→T or G→A", ("CT", "GA"), "#2c5d8a"), ("other G:C", None, "#9fc0d8"), ("A:T", "AT", "#d9952b")]
    for i, g in enumerate(("AID", "R174ER178E")):
        vs = [v for v in p["variants"] if v["group"] == g]
        n = len(vs)
        ctga = sum(v["ref"] + v["alt"] in ("CT", "GA") for v in vs)
        gc_other = sum(v["ref"] in "GC" for v in vs) - ctga
        at = n - ctga - gc_other
        bottom = 0
        for (lab, _k, col), c in zip(cats, (ctga, gc_other, at)):
            ax.bar(i, 100 * c / n, 0.6, bottom=bottom, color=col, label=lab if i == 0 else None)
            bottom += 100 * c / n
        ax.text(i, 103, f"n = {n}\nS5F {s5[g]['index']:.2f}", ha="center", fontsize=6, color=MID)
    ax.set_xticks([0, 1], ["AID", "R174E/\nR178E"])
    ax.set_ylim(0, 125)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Share of calls (%)")
    titled(ax, "Spectrum and targeting", "do the mutant's calls look like AID's?")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, handlelength=1, columnspacing=0.7)
    letter(ax, "C", dx=-0.42)
    # D
    ax = fig.add_subplot(gs[1, 0:2])
    p = json.load(open(PUB / "hmces/project.json"))
    groups = [("WT", "Wt", CASEF, CASE), ("HMCES_KO", "HMCES-KO", "#b9a3d6", "#6a4c93"), ("AID_KO", "AID-KO", CTRLF, CTRL)]
    for j, kind in enumerate(("deletions", "insertions")):
        for i, (g, lab, f, e) in enumerate(groups):
            v = [s[kind] for s in p["samples"] if s["group"] == g]
            x = j * 3.9 + i * 1.15      # spaced so the 6-pt group labels clear each other
            ax.bar(x, np.mean(v), 0.75, color=f, edgecolor=e, lw=0.6)
            ax.scatter(x + np.linspace(-0.2, 0.2, len(v)), v, s=7, color=e, zorder=3)
    ax.set_xticks([0, 1.15, 2.3, 3.9, 5.05, 6.2], ["Wt", "HMCES\nKO", "AID\nKO"] * 2, fontsize=6)
    ax.text(1.15, -0.36, "deletions", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, fontweight="bold")
    ax.text(5.05, -0.36, "insertions", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, fontweight="bold")
    ax.set_ylabel("Events per clone (20,000 reads)")
    titled(ax, "Indels in Ramos, Illumina (Wu et al.)", "indels, not only substitutions")
    letter(ax, "D", dx=-0.3, dy=1.1)
    # E
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
                      for lab, b in zip(["3–4", "5–9", "10–19", "≥20"], bins)])
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
    titled(ax, "UMI consensus, R10 (Karst et al.)", f"error after consensus, {u['molecules_3plus_reads']} molecules")
    ax.legend(fontsize=6, loc="lower left", title="consensus", title_fontsize=6, handlelength=1.4)
    json.dump({"consensus_error_pct": round(cons, 4), "ci": [round(lo, 4), round(hi, 4)],
               "fold": round(raw / cons, 1), "fold_ci": [round(raw / hi, 1), round(raw / lo, 1)]},
              open(DATA / "fig7_umi.json", "w"), indent=1)
    letter(ax, "E", dx=-0.26, dy=1.1)
    # F: the dead-AID control reveals the UNG-dependent spectrum
    ax = fig.add_subplot(gs[1, 4:6])
    mo = json.load(open(PUB / "morange/project.json"))
    au = json.load(open(PUB / "morange/authors_counts.json"))["libraries"]
    arms = ("Ugi", "Control", "UNG")
    V = defaultdict(list)
    for v in mo["variants"]:
        V[v["sample"]].append(v)
    def ts_share(vs):
        gc = [v for v in vs if v["ref"] in "CG"]
        return 100 * sum((v["ref"], v["alt"]) in (("C", "T"), ("G", "A")) for v in gc) / len(gc) if gc else np.nan

    # one value per experiment and arm; lines join an experiment's three arms (the pairing)
    E = (1, 2, 3)
    ours = np.array([[ts_share(V[f"Exp{e}_{arm}_AID"]) for arm in arms] for e in E])
    theirs = np.array([[au[f"Exp{e}_{arm}_AID"]["pct_transition_at_GC"] for arm in arms] for e in E])
    dead = np.array([[au[f"Exp{e}_{arm}_AIDE58Q"]["pct_transition_at_GC"] for arm in arms] for e in E])
    x = np.arange(3)
    for k in range(3):
        ax.plot(x - 0.08, theirs[k], color=MID, lw=0.5, alpha=0.6)
        ax.plot(x + 0.08, ours[k], color=CASE, lw=0.5, alpha=0.6)
    ax.scatter(np.repeat(x, 3) + 0.22, dead.T.ravel(), marker="x", s=9, color=CTRL, lw=0.7, zorder=3,
               label="dead AID (E58Q), every mismatch, Q > 20")
    ax.scatter(np.tile(x - 0.08, 3), theirs.ravel(), marker="s", s=9, facecolor="white", edgecolor=MID, lw=0.6,
               zorder=3, label="AID, every mismatch, Q > 20")
    ax.scatter(np.tile(x + 0.08, 3), ours.ravel(), marker="o", s=9, color=CASE, zorder=3, label="AID, SHMtrace calls")
    json.dump({"arms": arms, "shmtrace": ours.round(1).tolist(), "authors": theirs.round(1).tolist(),
               "dead": dead.round(1).tolist()}, open(DATA / "fig7_morange.json", "w"), indent=1)
    ax.set_xticks(x, ["Ugi\n(UNG off)", "empty", "UNG\n(extra)"], fontsize=6)
    ax.set_ylim(-4, 104)
    ax.set_ylabel("Transitions at C:G (%)")
    titled(ax, "Ectopic AID, non-B cells (Pérez-Durán)", "does the spectrum follow UNG?", size=6.6)
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=1, handlelength=1.2,
              handletextpad=0.3, labelspacing=0.3)
    letter(ax, "F", dx=-0.3, dy=1.1)
    save(fig, "Figure7")


if __name__ == "__main__":
    main()
