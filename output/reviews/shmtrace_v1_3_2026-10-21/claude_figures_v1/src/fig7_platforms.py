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

from style import SCRATCH  # noqa: E402
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

CAL = DATA.parents[1]
PUB = DATA / "public"


def titled(ax, title, question, size=6.6):
    """Panel title, with the dataset it comes from on a smaller line beneath."""
    ax.set_title(title, fontsize=size, pad=10)
    ax.text(0.0, 1.012, question, transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=MID)


def strip(ax, groups, samples, colors):
    """Libraries as points; zero-call libraries in a separate '0' strip below an axis break, never at a small
    positive value. The bar is the geometric mean with zeros entered at one event on that library's own exposure,
    the paper's zero-count rule for its tests; a group with no call at all gets no bar."""
    pos = [s["mutation_frequency_x1e3"] for g in groups for s in samples[g] if s["mutation_frequency_x1e3"] > 0]
    floors = [1e3 / s["frequency_bases"] for g in groups for s in samples[g]]
    lo = min(pos + floors)               # the one-event values sit inside the measured range, above the break
    zero_y = lo / 5
    for i, (g, (f, e)) in enumerate(zip(groups, colors)):
        v = np.array([s["mutation_frequency_x1e3"] for s in samples[g]], dtype=float)
        one = np.array([1e3 / s["frequency_bases"] for s in samples[g]])
        x = i + np.linspace(-0.15, 0.15, len(v))
        ax.scatter(x, np.where(v > 0, v, zero_y), s=16, facecolor=f, edgecolor=e, lw=0.6, zorder=3)
        if (v > 0).any():
            g_ = 10 ** np.mean(np.log10(np.where(v > 0, v, one)))
            ax.plot([i - 0.28, i + 0.28], [g_, g_], color=e, lw=1.2)
    ax.set_yscale("log")
    top = max(pos) * 3
    ax.set_ylim(zero_y / 2, top)
    ax.axhspan(zero_y / 2, zero_y * 2, color=SOFT, alpha=0.15, lw=0)
    for dy in (2.25, 2.6):
        ax.plot([-0.012, 0.012], [zero_y * dy] * 2, transform=ax.get_yaxis_transform(), color=INK, lw=0.6,
                clip_on=False)
    decades = [10.0 ** k for k in range(int(np.ceil(np.log10(lo))), int(np.floor(np.log10(top))) + 1)]
    ax.set_yticks([zero_y] + decades, ["0"] + [f"$10^{{{int(np.log10(d))}}}$" for d in decades])
    ax.minorticks_off()


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 124 * MM))
    gs = gridspec.GridSpec(2, 6, figure=fig, left=0.075, right=0.975, top=0.87, bottom=0.12, wspace=1.35, hspace=1.05)
    axes = {}
    # A
    ax = fig.add_subplot(gs[0, 0:2])
    p = json.load(open(PUB / "ldha/project.json"))
    groups = ["WT", "LDHA_KO", "AID_KO"]
    samples = {g: [s for s in p["samples"] if s["group"] == g] for g in groups}
    strip(ax, groups, samples, [(CASEF, CASE), ("#b9a3d6", "#6a4c93"), (CTRLF, CTRL)])
    ax.set_xticks(range(3), [f"Wt\n{len(samples['WT'])} mice", f"$Ldha$-deleted\n{len(samples['LDHA_KO'])} mice",
                             f"AID-null\n{len(samples['AID_KO'])} mice"])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    titled(ax, "J$_H$4 in mice", "Illumina, PRJNA973112")
    axes["A"] = ax
    letter(ax, "A", dx=-0.32)
    # B
    ax = fig.add_subplot(gs[0, 2:4])
    p = json.load(open(PUB / "ch12/project.json"))
    groups = ["AID", "R50E", "R174ER178E", "EV"]
    samples = {g: [s for s in p["samples"] if s["group"] == g] for g in groups}
    strip(ax, groups, samples, [(CASEF, CASE), ("#d9dde2", MID), ("#f3c78e", "#a3621b"), (CTRLF, CTRL)])
    ax.set_xticks(range(4), ["AID", "R50E", "R174E/\nR178E", "empty\nvector"])
    ax.set_ylabel("Mutations per 10$^3$ bp")
    titled(ax, "AID variants in CH12F3", "Illumina, PRJNA918596; 4 replicates each")
    axes["B"] = ax
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
            if c:
                ax.text(i, bottom + 50 * c / n, f"{c}", ha="center", va="center", fontsize=6.5,
                        color="white" if col != "#9fc0d8" else INK, fontweight="bold")
            bottom += 100 * c / n
        ax.text(i, 103, f"{n} calls\nS5F index {s5[g]['index']:.2f}", ha="center", fontsize=6, color=MID)
    ax.set_xticks([0, 1], ["AID", "R174E/\nR178E"])
    ax.set_ylim(0, 125)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Share of calls (%)")
    titled(ax, "AID against R174E/R178E", "CH12F3 calls; numbers in bars are calls")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, handlelength=1, columnspacing=0.7)
    axes["C"] = ax
    letter(ax, "C", dx=-0.3)
    # D
    ax = fig.add_subplot(gs[1, 0:3])
    p = json.load(open(PUB / "hmces/project.json"))
    groups = [("WT", "Wt", CASEF, CASE), ("HMCES_KO", "HMCES-KO", "#b9a3d6", "#6a4c93"), ("AID_KO", "AID-KO", CTRLF, CTRL)]
    for j, kind in enumerate(("deletions", "insertions")):
        for i, (g, lab, f, e) in enumerate(groups):
            v = [s[kind] for s in p["samples"] if s["group"] == g]
            x = j * 3.9 + i * 1.15      # spaced so the 6-pt group labels clear each other
            ax.bar(x, np.mean(v), 0.75, color=f, alpha=0.35, edgecolor="none")
            ax.plot([x - 0.37, x + 0.37], [np.mean(v)] * 2, color=e, lw=1.0)
            ax.scatter(x + np.linspace(-0.22, 0.22, len(v)), v, s=16, facecolor=f, edgecolor=e, lw=0.6, zorder=3)
    ax.set_xticks([0, 1.15, 2.3, 3.9, 5.05, 6.2], ["Wt", "HMCES\nKO", "AID\nKO"] * 2, fontsize=6)
    ax.text(1.15, -0.36, "deletions", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, fontweight="bold")
    ax.text(5.05, -0.36, "insertions", transform=ax.get_xaxis_transform(), ha="center", fontsize=6, fontweight="bold")
    ax.set_ylabel("Indels called per clone library")
    titled(ax, "Indels in Ramos clones", "Illumina, PRJNA820976; every library at 20,000 reads")
    axes["D"] = ax
    letter(ax, "D", dx=-0.12, dy=1.1)
    # E: the dead-AID control reveals the UNG-dependent spectrum
    ax = fig.add_subplot(gs[1, 3:6])
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
    series = [(ours, CASE, "o", 0.0, "AID,\nSHMtrace calls", 4),
              (theirs, MID, "s", -0.12, "AID, all Q>20\nmismatches", -4),
              (dead, CTRL, "X", 0.12, "dead AID (E58Q),\nall Q>20 mismatches", 5)]
    for vals, col, mk, dx, lab, dy in series:
        for k in range(3):
            ax.plot(x + dx, vals[k], color=col, lw=0.6, alpha=0.7, zorder=2)
        ax.scatter(np.tile(x + dx, 3), vals.ravel(), marker=mk, s=14, facecolor="white" if mk == "s" else col,
                   edgecolor=col, lw=0.6, zorder=3)
        ax.text(2.25, np.nanmean(vals[:, 2]) + dy, lab, fontsize=6, color=col, va="center", ha="left",
                linespacing=1.1)
    ax.text(0.0, -0.36, "each line joins one experiment's three arms (3 experiments)", transform=ax.transAxes,
            fontsize=6, color=MID, ha="left")
    json.dump({"arms": arms, "shmtrace": ours.round(1).tolist(), "authors": theirs.round(1).tolist(),
               "dead": dead.round(1).tolist()}, open(SCRATCH / "fig7_morange.json", "w"), indent=1)
    ax.set_xticks(x, ["Ugi\n(UNG blocked)", "empty", "UNG\n(overexpressed)"], fontsize=6)
    ax.set_xlim(-0.35, 3.25)
    ax.set_ylim(-4, 104)
    ax.set_ylabel("Transitions among C:G changes (%)")
    titled(ax, "Repair dependence of ectopic AID", "mOrange reporter, PRJNA170048")
    axes["E"] = ax
    letter(ax, "E", dx=-0.12, dy=1.1)
    # one header per question, so the five panels read as three tests rather than five mini-figures
    for keys, text in ((("A", "B"), "Mutation burden: does the AID-dependent signal transfer?"),
                       (("C",), "Mutational signature"),
                       (("D",), "Indels: beyond substitutions"), (("E",), "Repair: does the spectrum follow UNG?")):
        x0 = axes[keys[0]].get_position().x0 - 0.04
        x1 = axes[keys[-1]].get_position().x1
        y = axes[keys[0]].get_position().y1 + (0.085 if keys[0] in "ABC" else 0.075)
        fig.add_artist(plt.Line2D([x0, x1], [y - 0.006, y - 0.006], color=SOFT, lw=0.6, transform=fig.transFigure))
        fig.text(x0, y, text, fontsize=7, fontweight="bold", color=INK, ha="left", va="bottom")
    save(fig, "Figure7")


if __name__ == "__main__":
    main()
