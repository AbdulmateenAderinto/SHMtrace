"""Figure S8, the matplotlib draft: a published Sanger study re-analysed automatically
(Barbulescu et al., 2024, J_H4, cohort 2).

The published Figure S8 is built in SciCanvas (scicanvas/figS8_cohort2.js), from the numbers that
scicanvas/cohort2_data.py takes from this module; running this file draws the draft into out/.

A    mutation frequency per mouse: SHMtrace from the chromatograms against the published values
B    the same per mouse, SHMtrace against published (log-log), with the identity line
C    mutations per clone, binned as in the published Fig. 2g
D    G:C transitions and transversions per mouse (% of each mouse's mutations), as in Fig. 2i
E    pooled 12-type substitution spectrum, SHMtrace against published, per genotype
F    the published comparisons recovered: fold, one-tailed t-test as published, and run time
Published values are replotted from the paper's Supplementary Fig. 9b tables and Fig. 2g percentages
(published/barbulescu2024_jh4.json); SHMtrace's come from out/sanger_cohort2.json (sanger_cohort2.py).
"""
import json
from pathlib import Path

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy import stats
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, letter, setup

PUB = json.load(open(Path(__file__).resolve().parents[1] / "published/barbulescu2024_jh4.json"))
RES = json.load(open(DATA / "sanger_cohort2.json"))
SH_C, PUB_C = "#1f6f8b", "#a3621b"                    # SHMtrace and the published (manual) analysis
MAEA, MAEAF = "#4f7a3a", "#9cc28a"
GCOL = {"WT": (CASE, CASEF), "Maea_KO": (MAEA, MAEAF), "Aicda_KO": (CTRL, CTRLF)}
GLAB = {"WT": "Wt", "Maea_KO": "Maea$^{fl/fl}$\nMb1$^{cre}$", "Aicda_KO": "Aicda$^{-/-}$"}
B = "GCAT"
TYPES = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
TS, TV = ("C>T", "G>A"), ("C>A", "C>G", "G>C", "G>T")
MICE = ["WT_2w1", "WT_2w2", "WT_2w3", "WT_2w4", "MaeaKO_2k1", "MaeaKO_2k2", "MaeaKO_2k3", "MaeaKO_2k4"]


def pub_spec(m):
    mat = PUB["mice"][m]["matrix"]
    return {f"{B[i]}>{B[j]}": mat[i][j] for i in range(4) for j in range(4) if i != j}


def pub_freq(m):
    x = PUB["mice"][m]
    return 1e4 * x["unique_mutations"] / x["nt"]


def gc(sp):
    tot = sum(sp.values())
    return 100 * sum(sp.get(k, 0) for k in TS) / tot, 100 * sum(sp.get(k, 0) for k in TV) / tot


def one_tailed(a, b):
    return stats.ttest_ind(a, b, alternative="greater").pvalue


def _bg_share(group):
    """The share of a group's calls that the AID-null clones' own error rate accounts for."""
    v = [RES["per_mouse"][m]["expected_false_pct"] for m in RES["per_mouse"]
         if RES["per_mouse"][m]["group"] == group and RES["per_mouse"][m]["expected_false_pct"] is not None]
    return f"{min(v):.0f}-{max(v):.0f}%" if v else "-"


def main():
    setup()
    per = RES["per_mouse"]
    fig = plt.figure(figsize=(FULL, 118 * MM))
    gs = gridspec.GridSpec(2, 3, figure=fig, left=0.07, right=0.985, top=0.93, bottom=0.1, wspace=0.42, hspace=0.62,
                           width_ratios=[1.25, 1, 1.15])

    # A: frequency per mouse, SHMtrace and published side by side
    ax = fig.add_subplot(gs[0, 0])
    groups = ["WT", "Maea_KO", "Aicda_KO"]
    for i, g in enumerate(groups):
        e, f = GCOL[g]
        sh = [per[m]["freq_x1e4"] for m in per if per[m]["group"] == g]
        x0 = i * 1.0
        ax.scatter(x0 - 0.17 + np.linspace(-0.06, 0.06, len(sh)), sh, s=13, facecolor=f, edgecolor=e, lw=0.6, zorder=3)
        ax.plot([x0 - 0.3, x0 - 0.04], [np.mean(sh)] * 2, color=SH_C, lw=1.3)
        if g != "Aicda_KO":
            pb = [pub_freq(m) for m in MICE if PUB["mice"][m]["group"] == g]
            ax.scatter(x0 + 0.17 + np.linspace(-0.06, 0.06, len(pb)), pb, s=13, facecolor="white", edgecolor=e, lw=0.6,
                       marker="s", zorder=3)
            ax.plot([x0 + 0.04, x0 + 0.3], [np.mean(pb)] * 2, color=PUB_C, lw=1.3)
    ax.set_xticks(range(3), [GLAB[g] for g in groups])
    ax.set_ylabel("Mutations per 10$^4$ bp")
    ax.set_ylim(0, 230)
    ax.set_title("Frequency per mouse", fontsize=7)
    ax.scatter([], [], s=13, facecolor=SOFT, edgecolor=INK, lw=0.6, label="SHMtrace, chromatograms")
    ax.scatter([], [], s=13, facecolor="white", edgecolor=INK, lw=0.6, marker="s", label="Published")
    ax.legend(loc="upper right", fontsize=5.3, handletextpad=0.2, borderaxespad=0.1)
    letter(ax, "A", dx=-0.22)

    # B: per mouse, SHMtrace against published
    ax = fig.add_subplot(gs[0, 1])
    xs = [pub_freq(m) for m in MICE]
    ys = [per[m]["freq_x1e4"] for m in MICE]
    for m, x, y in zip(MICE, xs, ys, strict=True):
        e, f = GCOL[PUB["mice"][m]["group"]]
        ax.scatter(x, y, s=16, facecolor=f, edgecolor=e, lw=0.6, zorder=3)
    lim = (8, 400)
    ax.plot(lim, lim, color=SOFT, lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    r, _ = stats.spearmanr(xs, ys)
    ratio = 10 ** np.mean(np.log10(np.array(ys) / np.array(xs)))
    ax.text(0.04, 0.96, f"8 mice, Spearman ρ = {r:.2f}\nSHMtrace/published = {ratio:.2f}× (geometric)",
            transform=ax.transAxes, va="top", fontsize=5.3, color=MID)
    ax.set_xlabel("Published (per 10$^4$ bp)")
    ax.set_ylabel("SHMtrace (per 10$^4$ bp)")
    ax.set_title("Agreement by mouse", fontsize=7)
    letter(ax, "B", dx=-0.3)

    # C: mutations per clone, binned as in Fig. 2g
    ax = fig.add_subplot(gs[0, 2])
    bins = PUB["fig2g_bins"]
    shades = ["#f2f2f2", "#d4e1ee", "#a9c3dd", "#6f99c4", "#3f6fa4", "#1d4377"]
    rows = []
    for g in ("WT", "Maea_KO"):
        cnt = np.sum([per[m]["load_bins"] for m in per if per[m]["group"] == g], axis=0)
        rows += [(f"{GLAB[g].replace(chr(10), ' ')}\npublished (n = {PUB['fig2g_n'][g]})", np.array(PUB["fig2g_percent"][g], float)),
                 (f"SHMtrace (n = {int(cnt.sum())})", 100 * cnt / cnt.sum())]
    for k, (_lab, pc) in enumerate(rows):
        left = 0
        y = [4.4, 3.5, 1.6, 0.7][k]
        for _b, v, c in zip(bins, pc, shades, strict=True):
            ax.barh(y, v, left=left, color=c, edgecolor=INK, lw=0.3, height=0.78)
            if v >= 9:
                ax.text(left + v / 2, y, f"{v:.0f}", ha="center", va="center", fontsize=5,
                        color="white" if c in shades[4:] else INK)
            left += v
    ax.set_yticks([4.4, 3.5, 1.6, 0.7], [r[0] for r in rows], fontsize=5.3)
    ax.set_ylim(0.1, 5.0)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of clones")
    for lab_b, c in zip(bins, shades, strict=True):
        ax.bar(0, 0, color=c, edgecolor=INK, lw=0.3, label=lab_b)
    ax.legend(title="Mutations per clone", ncol=6, loc="lower center", bbox_to_anchor=(0.45, 1.0), fontsize=5,
              title_fontsize=5.3, handlelength=0.9, columnspacing=0.6, handletextpad=0.3)
    letter(ax, "C", dx=-0.62)

    # D: G:C transitions and transversions per mouse
    ax = fig.add_subplot(gs[1, 0])
    pos = 0
    ticks = []
    for g in ("WT", "Maea_KO"):
        e, f = GCOL[g]
        for kind, idx in (("Ts", 0), ("Tv", 1)):
            sh = [gc(per[m]["spectrum"])[idx] for m in per if per[m]["group"] == g]
            pb = [gc(pub_spec(m))[idx] for m in MICE if PUB["mice"][m]["group"] == g]
            for dx, vals, mk, fc, col in ((-0.17, sh, "o", f, SH_C), (0.17, pb, "s", "white", PUB_C)):
                ax.scatter(pos + dx + np.linspace(-0.05, 0.05, len(vals)), vals, s=11, marker=mk, facecolor=fc,
                           edgecolor=e, lw=0.6, zorder=3)
                ax.plot([pos + dx - 0.12, pos + dx + 0.12], [np.mean(vals)] * 2, color=col, lw=1.3)
            ticks.append((pos, f"{kind}\n{GLAB[g].split(chr(10))[0]}"))
            pos += 1
    ax.set_xticks([t[0] for t in ticks], [t[1] for t in ticks])
    ax.set_ylim(0, 35)
    ax.set_ylabel("% of mutations at G:C")
    ax.set_title("G:C transitions and transversions", fontsize=7)
    letter(ax, "D", dx=-0.22)

    # E: pooled spectrum, SHMtrace against published
    ax = fig.add_subplot(gs[1, 1])
    allx, ally = [], []
    for g in ("WT", "Maea_KO"):
        e, f = GCOL[g]
        sh = {t: sum(per[m]["spectrum"].get(t, 0) for m in per if per[m]["group"] == g) for t in TYPES}
        pb = {t: sum(pub_spec(m)[t] for m in MICE if PUB["mice"][m]["group"] == g) for t in TYPES}
        x = [100 * pb[t] / sum(pb.values()) for t in TYPES]
        y = [100 * sh[t] / sum(sh.values()) for t in TYPES]
        allx += x
        ally += y
        ax.scatter(x, y, s=13, facecolor=f, edgecolor=e, lw=0.6, zorder=3, label=GLAB[g].replace("\n", " "))
        for t, a, b in zip(TYPES, x, y, strict=True):
            if a > 13 or b > 13:
                off = {"A>T": (-14, 4), "T>C": (3, -7)}.get(t, (3, -3))
                ax.annotate(t.replace(">", "→"), (a, b), xytext=off, textcoords="offset points", fontsize=4.8, color=MID)
    ax.plot([0, 30], [0, 30], color=SOFT, lw=0.6, ls="--")
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 30)
    r = stats.pearsonr(allx, ally)[0]
    ax.text(0.04, 0.96, f"r = {r:.2f}", transform=ax.transAxes, va="top", fontsize=5.5, color=MID)
    ax.set_xlabel("Published (% of mutations)")
    ax.set_ylabel("SHMtrace (% of mutations)")
    ax.set_title("Substitution spectrum", fontsize=7)
    ax.legend(loc="lower right", fontsize=5.2, handletextpad=0.2)
    letter(ax, "E", dx=-0.3)

    # F: the published comparisons recovered
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    wt = [m for m in per if per[m]["group"] == "WT"]
    mk = [m for m in per if per[m]["group"] == "Maea_KO"]
    sh_f = (np.mean([per[m]["freq_x1e4"] for m in wt]) / np.mean([per[m]["freq_x1e4"] for m in mk]),
            one_tailed([per[m]["freq_x1e4"] for m in wt], [per[m]["freq_x1e4"] for m in mk]))
    pw = [pub_freq(m) for m in MICE[:4]]
    pm = [pub_freq(m) for m in MICE[4:]]
    pb_f = (np.mean(pw) / np.mean(pm), PUB["fig2h_p"]["WT_vs_Maea_KO"])
    ts_sh = stats.ttest_ind([gc(per[m]["spectrum"])[0] for m in wt], [gc(per[m]["spectrum"])[0] for m in mk]).pvalue
    tv_sh = stats.ttest_ind([gc(per[m]["spectrum"])[1] for m in wt], [gc(per[m]["spectrum"])[1] for m in mk]).pvalue
    lines = [("", "Published", "SHMtrace"),
             ("Wt / Maea frequency", f"{pb_f[0]:.1f}×", f"{sh_f[0]:.1f}×"),
             ("  P (one-tailed t)", f"{pb_f[1]:.3f}", f"{sh_f[1]:.4f}"),
             ("G:C Ts, Wt vs Maea P", f"{PUB['fig2i_p']['Ts']:.4f}", f"{ts_sh:.4f}"),
             ("G:C Tv, Wt vs Maea P", f"{PUB['fig2i_p']['Tv']:.4f}", f"{tv_sh:.3f}"),
             ("Mice, Wt / Maea / null", "4 / 4 / 3", "4 / 4 / 4"),
             ("Clones (Wt and Maea)", str(sum(PUB['mice'][m]['sequences'] for m in MICE)),
              str(sum(per[m]['clones'] for m in per if per[m]['group'] != 'Aicda_KO'))),
             ("Background explains, Wt", "–", _bg_share("WT")),
             ("  ...   Maea", "–", _bg_share("Maea_KO")),
             ("Shared-clone flags", "–", f"{sum(1 for f in RES['project']['qc']['flags'] if f['check'] == 'duplicate_clone')}"),
             ("Analysis", "hand-edited", f"{RES['seconds']:.0f} s")]
    for k, (a, b, c) in enumerate(lines):
        y = 0.96 - k * 0.096
        w = "bold" if k == 0 else "normal"
        kw = dict(fontsize=5.6, fontweight=w, va="top", transform=ax.transAxes, clip_on=False)
        ax.text(0.0, y, a, color=INK, **kw)
        ax.text(0.66, y, b, color=PUB_C if k else INK, ha="center", **kw)
        ax.text(0.93, y, c, color=SH_C if k else INK, ha="center", **kw)
    ax.plot([0, 1.0], [0.875, 0.875], color=SOFT, lw=0.5, transform=ax.transAxes)
    ax.set_title("Published tests, recovered", fontsize=7)
    letter(ax, "F", dx=-0.05)
    fig.savefig(DATA / "FigureS8_sanger_cohort2_matplotlib.png", dpi=300)
    return {"spearman": float(stats.spearmanr(xs, ys)[0]), "ratio": float(ratio), "spec_r": float(r),
            "sh_fold": float(sh_f[0]), "sh_p": float(sh_f[1]), "pub_fold": float(pb_f[0]), "ts_p": float(ts_sh),
            "tv_p": float(tv_sh)}


if __name__ == "__main__":
    print(main())
