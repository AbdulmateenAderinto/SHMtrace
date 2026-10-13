"""Sanger clones against nanopore reads in the same wild-type mice: do they agree?

Runs both platforms through shmtrace, then draws the manuscript's comparison
figures (Figure 2A/B, Figure 3A-E, Supplementary S2/S3) from those calls, plus
tests the manuscript did not have. Writes figures, stats.json and stats.md.

    python docs/calibration/sanger_vs_nanopore.py [--ont-mode recommended|geneious]
        [--sanger-dir ~/Nanoshm/sanger_ab1_by_mouse] [--out docs/calibration/figures/sanger_vs_nanopore]

Sanger input: one folder per mouse of .ab1 traces (plus a FASTA of clones whose
traces are missing), built by organize_jh4_ab1.py. Nanopore: the six bundled libraries.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Circle, Patch  # noqa: E402
from scipy import stats  # noqa: E402

from shmtrace import motifs, project  # noqa: E402
from shmtrace.io_utils import read_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REF = str(ROOT / "src/shmtrace/references/mouse_jh4_intron.fasta")
CON = [str(ROOT / "src/shmtrace/references/mouse_jh4_contaminant.fasta")]
EX = ROOT / "examples_compare"
WT, KO_ONT, KO_SANGER = ["WT2", "WT3", "WT4"], ["KO1", "KO2", "KO3"], ["A1", "A2", "A3"]
LO, HI = 101, 629                      # the 529-bp comparison window, 1-based
BP = HI - LO + 1

# the manuscript's palette and conventions (figure_scripts/figure3_concordance.py)
SANGER = ("#0072B2", "#005A8D")
ONT = ("#CC79A7", "#A85C86")
WTC, KOC = ("#6D9AD0", "#274B7E"), ("#D07A7E", "#9A3439")
SEQ = LinearSegmentedColormap.from_list("seqblue", ["#F5F8FC", "#C4D8EE", "#6FA0CE", "#2F6BA6", "#173D63"])
PIE = ["#FFFFFF", "#DCEAF7", "#A8CBE8", "#5A9BD5", "#1F4E79"]
PIE_LAB = ["0", "1–3", "4–10", "11–20", "20+"]
GRID = "#C9C9C9"
BASES = ["A", "T", "C", "G"]
MUT = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": "black",
                     "savefig.dpi": 300, "figure.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})


# ------------------------------------------------------------------ data

def sanger_path(sanger_dir: Path, n: str) -> Path | None:
    """A mouse's clones: a folder of .ab1 traces (organize_jh4_ab1.py) or one FASTA."""
    for cand in (sanger_dir / n, sanger_dir / f"Aicda_{n}", sanger_dir / f"{n}.fasta"):
        if cand.exists():
            return cand
    return None


# the deposit's realised 500-read subsets (--ont-dir), named by library
DEPOSIT_NAMES = {"WT2": "DHJGKB_2_Wildtype_2", "WT3": "DHJGKB_3_Wildtype_3", "WT4": "DHJGKB_4_Wildtype_4",
                 "KO1": "DHJGKB_5_AID_KO_1", "KO2": "DHJGKB_6_AID_KO_2", "KO3": "DHJGKB_7_AID_KO_3"}


def ont_path(n: str, ont_dir: Path | None) -> Path:
    """examples_compare/, the deposit's delivered reads (Original_FASTQ/), or its 500-read subsets."""
    if ont_dir is None:
        return EX / f"{n}.fastq"
    if ont_dir.name == "Original_FASTQ":
        return ont_dir / f"{DEPOSIT_NAMES[n]}.fastq"
    return ont_dir / f"{DEPOSIT_NAMES[n]}_C57BL6_sub500.fastq"


def run_platforms(sanger_dir: Path, ont_mode: str, ont_dir: Path | None = None) -> dict:
    """Both platforms through shmtrace. Returns per-sample calls and loads."""
    out = {}
    sp = {"name": "sanger", "region": "jh4_intron", "platform": "sanger", "reference": REF,
          "contaminants": CON,
          "samples": [{"name": n, "path": str(sanger_path(sanger_dir, n)),
                       "group": "WT" if n in WT else "KO", "control": n not in WT}
                      for n in WT + KO_SANGER]}
    op = {"name": "ont", "region": "jh4_intron", "mode": ont_mode, "reference": REF, "depth": 500,
          "contaminants": CON,
          "samples": [{"name": n, "path": str(ont_path(n, ont_dir)),
                       "group": "WT" if n in WT else "KO", "control": n not in WT}
                      for n in WT + KO_ONT]}
    for plat, proj in (("sanger", sp), ("ont", op)):
        res = project.run_project(proj)
        for s in res["samples"]:
            d = res["_details"][s["sample"]]
            vs = [v for v in d["variants"] if LO <= v.pos <= HI]
            out[(plat, s["sample"])] = {
                "n": s["reads"],
                "variants": [(v.pos, v.ref, v.alt) for v in vs],
                "per_seq": [r["mutations"] for r in d["per_read"]],
                # bases each sequence read inside the window at good quality: a short but
                # clean Sanger clone counts what it read, not the full 529 bp
                "bp": sum(r["bases"] for r in d["per_read"]),
                "warnings": s.get("warnings") or "",
                "clone_weighted_x1e3": s.get("clone_weighted_frequency_x1e3"),
            }
    return out


def density(d: dict) -> tuple[float, bool]:
    """Unique mutations per base read, x 1e5: the manuscript's density, with the bases each
    sequence actually read in place of sequences x 529 (the same for full-length reads).

    A sample with none is set to one mutation's worth, as the manuscript does,
    which can only make it look more mutated. Returns (value, censored).
    """
    k = len(d["variants"])
    bp = d.get("bp") or d["n"] * BP
    return 1e5 * max(k, 1) / bp, k == 0


def matrix(variants) -> np.ndarray:
    """From (rows) x To (cols), unique mutations."""
    m = np.zeros((4, 4))
    for _p, r, a in variants:
        if r in BASES and a in BASES:
            m[BASES.index(r), BASES.index(a)] += 1
    return m


def props12(variants) -> np.ndarray:
    c = Counter(f"{r}>{a}" for _p, r, a in variants)
    tot = sum(c[t] for t in MUT)
    return np.array([100 * c[t] / tot if tot else np.nan for t in MUT])


def gc_ts_tv(variants) -> tuple[float, float]:
    """G:C transitions and transversions, % of all unique mutations."""
    n = len(variants)
    ts = sum(1 for _p, r, a in variants if (r, a) in (("C", "T"), ("G", "A")))
    tv = sum(1 for _p, r, a in variants if r in "CG" and (r, a) not in (("C", "T"), ("G", "A")))
    return 100 * ts / n, 100 * tv / n


def hotspot(variants, ref: str) -> tuple[float, float]:
    gc = [p for p, r, _a in variants if r in "GC"]
    at = [p for p, r, _a in variants if r in "AT"]
    h = sum(1 for p in gc if motifs.is_site("WRC", ref, p - 1) or motifs.is_site("GYW", ref, p - 1))
    w = sum(1 for p in at if motifs.is_site("WA", ref, p - 1) or motifs.is_site("TW", ref, p - 1))
    return 100 * h / len(gc) if gc else np.nan, 100 * w / len(at) if at else np.nan


def load_bins(per_seq) -> np.ndarray:
    edges = [(0, 0), (1, 3), (4, 10), (11, 20), (21, 10 ** 6)]
    return np.array([sum(1 for m in per_seq if a <= m <= b) for a, b in edges], float)


# ------------------------------------------------------------------ statistics

def geo(xs) -> float:
    return float(10 ** np.mean(np.log10(xs)))


def boot_ci(xs, n=20000, seed=1) -> tuple[float, float]:
    rng = random.Random(seed)
    vals = sorted(geo([rng.choice(xs) for _ in xs]) for _ in range(n))
    return vals[int(0.025 * n)], vals[int(0.975 * n) - 1]


def welch_log(a, b) -> float:
    return float(stats.ttest_ind(np.log10(a), np.log10(b), equal_var=False).pvalue)


def rm_anova_interaction(y: np.ndarray) -> tuple[float, float, int, int]:
    """Platform x type interaction, two-way repeated measures (both within mouse).

    y[s, p, t]: subjects x platforms x types. Uncorrected (sphericity assumed).
    """
    S, P, T = y.shape
    g = y.mean()
    m_pt, m_ps, m_ts = y.mean(0), y.mean(2), y.mean(1)
    m_p, m_t, m_s = y.mean((0, 2)), y.mean((0, 1)), y.mean((1, 2))
    ss_pt = S * ((m_pt - m_p[:, None] - m_t[None, :] + g) ** 2).sum()
    resid = (y - m_pt[None] - m_ps[:, :, None] - m_ts[:, None, :]
             + m_p[None, :, None] + m_t[None, None, :] + m_s[:, None, None] - g)
    ss_err = (resid ** 2).sum()
    df1, df2 = (P - 1) * (T - 1), (P - 1) * (T - 1) * (S - 1)
    f = (ss_pt / df1) / (ss_err / df2)
    return float(f), float(stats.f.sf(f, df1, df2)), df1, df2


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson interval for a proportion k/n, in percent."""
    if n == 0:
        return (0.0, 0.0)
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def sidak(p: float, k: int) -> float:
    return 1 - (1 - p) ** k


def pearson_ci(x, y) -> dict:
    r, p = stats.pearsonr(x, y)
    z, se = math.atanh(r), 1 / math.sqrt(len(x) - 3)
    slope = float(np.polyfit(x, y, 1)[0])
    return {"r": float(r), "p": float(p), "ci": (math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)),
            "slope": slope}


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if a | b else float("nan")


# ------------------------------------------------------------------ drawing helpers

def dots(ax, x, vals, edge, censored=(), s=60, r=0.16):
    for i, v in enumerate(vals):
        cen = i in censored
        ax.scatter(x + (i - 1) * 0.07, v, marker="o", s=s, zorder=6, linewidths=1.6,
                   facecolors="white" if not cen else "#eeeeee", edgecolors=edge,
                   linestyle="--" if cen else "-")
    gm = geo(vals)
    gsd = 10 ** np.std(np.log10(vals), ddof=1)
    ax.plot([x - r, x + r], [gm, gm], color=edge, lw=2.2, zorder=7)
    ax.plot([x, x], [gm / gsd, gm * gsd], color=edge, lw=1.3, zorder=5)
    for yv in (gm / gsd, gm * gsd):
        ax.plot([x - r / 2, x + r / 2], [yv, yv], color=edge, lw=1.3, zorder=5)


def bracket(ax, x1, x2, y, text):
    ax.plot([x1, x1, x2, x2], [y, y * 1.15, y * 1.15, y], color="black", lw=1.1)
    ax.text((x1 + x2) / 2, y * 1.25, text, ha="center", va="bottom", fontsize=10)


def heat(ax, m, title, vmax, show_from=True):
    tot = m.sum()
    for i in range(4):
        for j in range(4):
            if i == j:
                ax.add_patch(plt.Rectangle((j, 3 - i), 1, 1, fc="#e9e9e9", ec="white", lw=2))
                ax.text(j + .5, 3 - i + .5, "—", ha="center", va="center", color="#888")
                continue
            pct = 100 * m[i, j] / tot if tot else 0
            fc = SEQ(min(pct / vmax, 1) if vmax else 0)
            ax.add_patch(plt.Rectangle((j, 3 - i), 1, 1, fc=fc, ec="white", lw=2))
            lum = sum(c * w for c, w in zip(to_rgb(fc), (.2126, .7152, .0722), strict=True))
            col = "white" if lum < .5 else "#1a1a1a"
            ax.text(j + .5, 3 - i + .62, f"{int(m[i, j])}", ha="center", va="center", fontsize=10,
                    fontweight="bold", color=col)
            if tot:      # undefined for a sample with no mutations
                ax.text(j + .5, 3 - i + .28, f"{pct:.1f}%", ha="center", va="center", fontsize=6.8,
                        color=col)
    for k, b in enumerate(BASES):
        ax.text(k + .5, 4.18, b, ha="center", va="center", fontweight="bold")
        ax.text(-.25, 3 - k + .5, b, ha="center", va="center", fontweight="bold")
    ax.text(2, 4.55, "To", ha="center", fontsize=9.5, fontweight="bold")
    if show_from:
        ax.text(-.7, 2, "From", ha="center", va="center", rotation=90, fontsize=9.5, fontweight="bold")
    ax.text(2, 4.85, title, ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_xlim(-.9, 4.1)
    ax.set_ylim(-.1, 5.9)
    ax.set_aspect("equal")
    ax.axis("off")


def save(fig, out: Path, name: str):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ figures

def fig_density(D, out, st):
    """Manuscript Figure 2A/B: Wt against Aicda-/-, each platform on its own axis."""
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.2))
    for ax, plat, kos, title in ((axes[0], "ont", KO_ONT, "A   ONT"), (axes[1], "sanger", KO_SANGER, "B   Sanger")):
        wt = [density(D[(plat, n)])[0] for n in WT]
        ko_d = [density(D[(plat, n)]) for n in kos]
        ko = [v for v, _c in ko_d]
        cens = {i for i, (_v, c) in enumerate(ko_d) if c}
        dots(ax, 0, wt, WTC[1])
        dots(ax, 1, ko, KOC[1], censored=cens)
        p = welch_log(wt, ko)
        fold = geo(wt) / geo(ko)
        top = max(wt) * 2
        bracket(ax, 0, 1, top, f"P = {p:.4f} (~{fold:.0f}×)")
        ax.set_yscale("log")
        ax.set_xticks([0, 1], ["$\\mathit{Wt}$", "$\\mathit{Aicda}^{-/-}$"])
        ax.set_xlim(-.6, 1.6)
        ax.set_ylabel("Mutations / bp (× 10$^{-5}$)")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.yaxis.grid(True, ls=(0, (5, 4)), color=GRID)
        ax.set_axisbelow(True)
        lo = min(wt + ko) / 4
        ax.set_ylim(lo, top * 6)
        st[f"density_{plat}"] = {
            "wt": dict(zip(WT, wt, strict=True)), "ko": dict(zip(kos, ko, strict=True)),
            "ko_censored": [kos[i] for i in cens],
            "wt_geomean": geo(wt), "wt_ci": boot_ci(wt), "ko_geomean": geo(ko), "ko_ci": boot_ci(ko),
            "fold": fold, "p_welch_log": p,
            "sequences": {n: D[(plat, n)]["n"] for n in WT + kos}}
    fig.text(0.5, -0.04, "Each dot one mouse; bar, geometric mean ± geometric SD; grey dashed dot: no "
             "mutations, shown at one mutation's worth. Welch's t on log10.", ha="center", fontsize=8.5,
             color="#444")
    fig.tight_layout()
    save(fig, out, "fig1_density_wt_vs_ko")


def fig_spectrum(D, out, st, ref):
    """Manuscript Figure 3A-E."""
    sv = {n: D[("sanger", n)]["variants"] for n in WT}
    ov = {n: D[("ont", n)]["variants"] for n in WT}
    MS, MO = sum(matrix(sv[n]) for n in WT), sum(matrix(ov[n]) for n in WT)
    vmax = max(100 * MS.max() / MS.sum(), 100 * MO.max() / MO.sum())

    fig = plt.figure(figsize=(16, 11.5))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=.32)
    gs_top = outer[0].subgridspec(1, 2, width_ratios=[1.25, 2], wspace=.12)
    gs_bot = outer[1].subgridspec(1, 3, width_ratios=[1, 1.05, 1.05], wspace=.35)
    # A matrices
    gA = gs_top[0, 0].subgridspec(2, 2, height_ratios=[1, .05], hspace=.08, wspace=.02)
    heat(fig.add_subplot(gA[0, 0]), MS, f"Wt (Sanger)\nn = {int(MS.sum())}", vmax)
    heat(fig.add_subplot(gA[0, 1]), MO, f"Wt (ONT)\nn = {int(MO.sum())}", vmax, show_from=False)
    cax = fig.add_subplot(gA[1, :])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, vmax)), cax=cax,
                      orientation="horizontal")
    cb.set_label("% of unique mutations", fontsize=9)
    fig.text(0.11, 0.89, "A", fontsize=15, fontweight="bold")

    # B 12 types
    axB = fig.add_subplot(gs_top[0, 1])
    Ps = np.array([props12(sv[n]) for n in WT])      # mice x types
    Po = np.array([props12(ov[n]) for n in WT])
    x = np.arange(12)
    pairs = []
    for i in range(12):
        p = stats.ttest_rel(Ps[:, i], Po[:, i]).pvalue
        pairs.append(sidak(float(p), 12) if p == p else float("nan"))
    # Bars: the platform's pooled proportion (every unique mutation of the three mice)
    # with a 95% Wilson interval. A mean of three mice would let a mouse with a handful
    # of mutations (Sanger Wt2 has 6) swing a type by tens of points. Points: each
    # mouse, area proportional to its number of mutations.
    for vv, arr, o, (fill, edge), lab in ((sv, Ps, -.2, SANGER, "Sanger"), (ov, Po, .2, ONT, "ONT")):
        allv = [m for n in WT for m in vv[n]]
        cnt = Counter(f"{r}>{a}" for _p, r, a in allv)
        N = sum(cnt[t] for t in MUT)
        pooled = np.array([100 * cnt[t] / N for t in MUT])
        ci = np.array([wilson(cnt[t], N) for t in MUT])
        axB.bar(x + o, pooled, .38, yerr=[pooled - ci[:, 0], ci[:, 1] - pooled], capsize=3,
                color=fill, edgecolor=edge, lw=1.4, error_kw=dict(ecolor=edge, lw=1.2),
                label=f"{lab} (n = {N})", zorder=2)
        sizes = [max(8, 70 * math.sqrt(len(vv[n]) / max(len(vv[m]) for m in WT))) for n in WT]
        for k in range(3):
            axB.scatter(x + o + (k - 1) * .07, arr[k], s=sizes[k], facecolor="white",
                        edgecolor=edge, lw=1.0, zorder=4, alpha=.9)
    top = min(np.nanmax(np.concatenate([Ps, Po])), 40) + 4
    for i, p in enumerate(pairs):
        axB.text(x[i], top, "adj. p\n" + (">0.99" if p > .99 else f"{p:.2f}"), ha="center", fontsize=7.5,
                 color="#333")
    yall = np.stack([Ps, Po], axis=1)
    F, p_int, df1, df2 = rm_anova_interaction(yall)
    axB.set_xticks(x, MUT)
    axB.set_ylabel("% of unique mutations")
    axB.set_ylim(0, top + 7)
    small = [n for n in WT if len(sv[n]) < 20]
    if small:
        axB.text(.99, .62, "Points: each mouse, area ∝ its mutations.\n" + ", ".join(
            f"Sanger {n.replace('WT', 'Wt')}: {len(sv[n])} mutations" for n in small)
            + "\n(the smallest points)", transform=axB.transAxes, ha="right",
            fontsize=8, color="#555")
    axB.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    axB.set_title(f"Platform × type interaction (repeated measures): F({df1},{df2}) = {F:.2f}, "
                  f"p = {p_int:.3f}", fontsize=10, loc="right", pad=26)
    axB.yaxis.grid(True, ls=(0, (5, 4)), color=GRID)
    axB.set_axisbelow(True)
    axB.text(-0.07, 1.06, "B", transform=axB.transAxes, fontsize=15, fontweight="bold")

    # C pooled scatter
    axC = fig.add_subplot(gs_bot[0, 0])
    allS = [v for n in WT for v in sv[n]]
    allO = [v for n in WT for v in ov[n]]
    ps, po = props12(allS), props12(allO)
    pc = pearson_ci(ps, po)
    lim = max(ps.max(), po.max()) * 1.12
    axC.plot([0, lim], [0, lim], ls="--", color="#999", lw=1.2)
    axC.scatter(ps, po, s=70, facecolor=ONT[0], edgecolor="white", lw=.9, zorder=3)
    # label the two dominant types and any that differ by more than 2 points
    for t, a, b in zip(MUT, ps, po, strict=True):
        if max(a, b) > 14 or abs(a - b) > 2:
            axC.annotate(t, (a, b), xytext=(6, -3), textcoords="offset points", fontsize=8.5)
    axC.set_xlim(0, lim)
    axC.set_ylim(0, lim)
    axC.set_aspect("equal")
    axC.set_xlabel("Sanger, % of unique mutations")
    axC.set_ylabel("ONT, % of unique mutations")
    axC.text(.04, .96, f"r = {pc['r']:.2f} [{pc['ci'][0]:.2f}–{pc['ci'][1]:.2f}]\nslope = {pc['slope']:.2f}",
             transform=axC.transAxes, va="top", fontsize=10)
    axC.text(-0.2, 1.06, "C", transform=axC.transAxes, fontsize=15, fontweight="bold")

    # D G:C Ts / Tv
    axD = fig.add_subplot(gs_bot[0, 1])
    tsv = {"sanger": np.array([gc_ts_tv(sv[n]) for n in WT]), "ont": np.array([gc_ts_tv(ov[n]) for n in WT])}
    ptt = [sidak(float(stats.ttest_rel(tsv["sanger"][:, k], tsv["ont"][:, k]).pvalue), 2) for k in (0, 1)]
    xx = np.arange(2)
    for plat, o, (fill, edge), lab in (("sanger", -.18, SANGER, "Sanger"), ("ont", .18, ONT, "ONT")):
        a = tsv[plat]
        allv = [m for n in WT for m in (sv if plat == "sanger" else ov)[n]]
        N = len(allv)
        ks = (sum(1 for _p, r, al in allv if (r, al) in (("C", "T"), ("G", "A"))),
              sum(1 for _p, r, al in allv if r in "CG" and (r, al) not in (("C", "T"), ("G", "A"))))
        pooled = np.array([100 * k / N for k in ks])
        ci = np.array([wilson(k, N) for k in ks])
        axD.bar(xx + o, pooled, .34, yerr=[pooled - ci[:, 0], ci[:, 1] - pooled], capsize=4,
                color=fill, edgecolor=edge, lw=1.4, error_kw=dict(ecolor=edge, lw=1.2), label=lab,
                zorder=2)
        for k in range(3):
            axD.scatter(xx + o + (k - 1) * .06, a[k], s=22, facecolor="white", edgecolor=edge, lw=1.2, zorder=4)
    for k in (0, 1):
        axD.text(xx[k], max(tsv["sanger"][:, k].max(), tsv["ont"][:, k].max()) + 3,
                 f"adj. p = {ptt[k]:.2f}", ha="center", fontsize=9)
    axD.set_xticks(xx, ["G:C transitions", "G:C transversions"])
    axD.set_ylabel("% of unique mutations")
    axD.legend(frameon=False, loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=2)
    axD.set_ylim(0, max(tsv["sanger"].max(), tsv["ont"].max()) + 8)
    axD.yaxis.grid(True, ls=(0, (5, 4)), color=GRID)
    axD.set_axisbelow(True)
    axD.text(-0.15, 1.06, "D", transform=axD.transAxes, fontsize=15, fontweight="bold")

    # E donuts
    gE = gs_bot[0, 2].subgridspec(2, 2, width_ratios=[1, .55], hspace=.12, wspace=0)
    bins = {}
    for k, plat in enumerate(("sanger", "ont")):
        a = fig.add_subplot(gE[k, 0])
        v = sum(load_bins(D[(plat, n)]["per_seq"]) for n in WT)
        bins[plat] = v
        N = int(v.sum())
        w, _ = a.pie(v, colors=PIE, startangle=90, counterclock=False, radius=1,
                     wedgeprops=dict(width=.42, edgecolor="#333", lw=.8))
        for wedge, val, fc in zip(w, v, PIE, strict=True):
            if val / N < .04:
                continue
            ang = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)
            lum = sum(c * ww for c, ww in zip(to_rgb(fc), (.2126, .7152, .0722), strict=True))
            a.text(.79 * np.cos(ang), .79 * np.sin(ang), f"{100 * val / N:.0f}%", ha="center",
                   va="center", fontsize=8.5, color="white" if lum < .5 else "#1a1a1a")
        a.text(0, .1, "Sanger" if plat == "sanger" else "ONT", ha="center", va="center",
               fontsize=9.5, fontweight="bold")
        a.text(0, -.17, f"n = {N:,}", ha="center", va="center", fontsize=8)
        a.set_aspect("equal")
        if k == 0:
            a.text(-1.6, 1.15, "E", fontsize=15, fontweight="bold")
    axl = fig.add_subplot(gE[:, 1])
    axl.axis("off")
    axl.legend([Patch(fc=c, ec="#333") for c in PIE], PIE_LAB, title="Mutations\nper sequence",
               loc="center left", frameon=False, fontsize=9, title_fontsize=9.5)
    save(fig, out, "fig2_spectrum_concordance")
    st["spectrum"] = {
        "sanger_unique": int(MS.sum()), "ont_unique": int(MO.sum()),
        "pooled_props": {"sanger": dict(zip(MUT, ps.round(2).tolist(), strict=True)),
                         "ont": dict(zip(MUT, po.round(2).tolist(), strict=True))},
        "pearson": pc, "mean_abs_diff_points": float(np.mean(np.abs(ps - po))),
        "max_abs_diff_points": float(np.max(np.abs(ps - po))),
        "rm_anova_interaction": {"F": F, "p": p_int, "df": (df1, df2)},
        "paired_sidak_min": float(np.nanmin(pairs)),
        "paired_sidak": dict(zip(MUT, pairs, strict=True)),
        "gc_ts_tv": {"sanger_mean": tsv["sanger"].mean(0).tolist(), "ont_mean": tsv["ont"].mean(0).tolist(),
                     "sidak_p": ptt},
        "per_sequence_bins": {k: v.tolist() for k, v in bins.items()},
        "per_sequence_mean": {plat: float(np.mean([m for n in WT for m in D[(plat, n)]["per_seq"]]))
                              for plat in ("sanger", "ont")},
    }


def fig_beyond(D, out, st, ref):
    """What the manuscript did not test: is each mouse's Sanger closer to its own
    nanopore than to another mouse's? Do the platforms agree on where along the
    intron, on hotspot preference, and on which mouse is most mutated?"""
    sv = {n: D[("sanger", n)]["variants"] for n in WT}
    ov = {n: D[("ont", n)]["variants"] for n in WT}
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # A: spectrum correlation, every Sanger mouse against every ONT mouse
    R = np.array([[stats.pearsonr(props12(sv[a]), props12(ov[b]))[0] for b in WT] for a in WT])
    ax = axes[0, 0]
    im = ax.imshow(R, cmap=SEQ, vmin=min(.5, R.min()), vmax=1)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=12,
                    fontweight="bold" if i == j else "normal", color="white" if R[i, j] > .85 else "#1a1a1a")
            if i == j:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, ec="#1a1a1a", lw=2.2))
    ax.set_xticks(range(3), [f"{n}\nONT" for n in WT])
    ax.set_yticks(range(3), [f"{n} Sanger" for n in WT])
    ax.set_title("Spectrum correlation, each mouse against each", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=.046, pad=.04, label="Pearson r, 12 types")
    ax.text(-0.3, 1.05, "A", transform=ax.transAxes, fontsize=15, fontweight="bold")

    def panel_where(ax, letter):
        # where along the window, pooled wild type
        edges = np.linspace(LO, HI + 1, 22)
        mids = (edges[:-1] + edges[1:]) / 2
        hs = np.histogram([p for n in WT for p, _r, _a in sv[n]], edges)[0]
        ho = np.histogram([p for n in WT for p, _r, _a in ov[n]], edges)[0]
        fs, fo = 100 * hs / hs.sum(), 100 * ho / ho.sum()
        rp = stats.pearsonr(fs, fo)
        ax.plot(mids, fs, "-o", color=SANGER[0], ms=4, label="Sanger")
        ax.plot(mids, fo, "-o", color=ONT[0], ms=4, label="ONT")
        ax.set_xlabel("Position in the JH4 reference (bp)")
        ax.set_ylabel("% of unique mutations, 25-bp bins")
        ax.legend(frameon=False)
        ax.set_title(f"Along the window: r = {rp[0]:.2f} (p = {rp[1]:.1e})", fontsize=11)
        ax.yaxis.grid(True, ls=(0, (5, 4)), color=GRID)
        ax.set_axisbelow(True)
        ax.text(-0.13, 1.05, letter, transform=ax.transAxes, fontsize=15, fontweight="bold")
        return rp

    def panel_hotspot(ax, letter):
        # hotspot preference per mouse
        H = {plat: np.array([hotspot(D[(plat, n)]["variants"], ref) for n in WT]) for plat in ("sanger", "ont")}
        xx = np.arange(2)
        for plat, o, (fill, edge), lab in (("sanger", -.18, SANGER, "Sanger"), ("ont", .18, ONT, "ONT")):
            a = H[plat]
            ax.bar(xx + o, a.mean(0), .34, color=fill, edgecolor=edge, lw=1.4, label=lab, zorder=2)
            for k in range(3):
                ax.scatter(xx + o + (k - 1) * .06, a[k], s=22, facecolor="white", edgecolor=edge, lw=1.2, zorder=4)
        # the share expected if mutations ignored motifs: the window's own share
        gc = [i for i in range(LO - 1, HI) if ref[i] in "GC"]
        at = [i for i in range(LO - 1, HI) if ref[i] in "AT"]
        exp = (100 * sum(1 for i in gc if motifs.is_site("WRC", ref, i) or motifs.is_site("GYW", ref, i)) / len(gc),
               100 * sum(1 for i in at if motifs.is_site("WA", ref, i) or motifs.is_site("TW", ref, i)) / len(at))
        for k in (0, 1):
            ax.plot([k - .4, k + .4], [exp[k], exp[k]], ls="--", color="#555", lw=1.2)
        ax.text(1.42, exp[1], "share of\nsites in\nwindow", fontsize=8, va="center", color="#555")
        ax.set_xticks(xx, ["G:C in WRC/GYW", "A:T in WA/TW"])
        ax.set_ylabel("% of that base's mutations")
        ax.legend(frameon=False, loc="upper left")
        ax.set_ylim(0, 100)
        ax.yaxis.grid(True, ls=(0, (5, 4)), color=GRID)
        ax.set_axisbelow(True)
        ax.text(-0.13, 1.05, letter, transform=ax.transAxes, fontsize=15, fontweight="bold")
        return H, exp

    rp = panel_where(axes[0, 1], "B")
    H, exp = panel_hotspot(axes[1, 0], "C")

    # D: mutation load per mouse, both platforms
    ax = axes[1, 1]
    fs_ = [1e3 * sum(D[("sanger", n)]["per_seq"]) / D[("sanger", n)]["bp"] for n in WT]
    fo_ = [1e3 * sum(D[("ont", n)]["per_seq"]) / D[("ont", n)]["bp"] for n in WT]
    fc_ = [D[("ont", n)]["clone_weighted_x1e3"] for n in WT]
    for n, a, b, c in zip(WT, fs_, fo_, fc_, strict=True):
        ax.scatter(a, b, s=90, color=ONT[0], edgecolor="white", zorder=3)
        ax.annotate(n, (a, b), xytext=(6, -12), textcoords="offset points")
        if c:
            ax.scatter(a, c, s=90, facecolor="white", edgecolor=ONT[1], lw=1.6, zorder=3)
            ax.annotate("", (a, c), (a, b), arrowprops=dict(arrowstyle="->", color=ONT[1], lw=1))
    ax.scatter([], [], s=90, color=ONT[0], label="ONT per read")
    ax.scatter([], [], s=90, facecolor="white", edgecolor=ONT[1], lw=1.6, label="ONT per clone (junction)")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    lim = (min(fs_ + fo_ + [c for c in fc_ if c]) / 2, max(fs_ + fo_ + [c for c in fc_ if c]) * 2)
    ax.plot(lim, lim, ls="--", color="#999")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Sanger: mutations per 10³ bp sequenced")
    ax.set_ylabel("ONT: mutations per 10³ bp sequenced")
    ax.set_title("Mutation load, mouse by mouse", fontsize=11)
    ax.text(-0.2, 1.05, "D", transform=ax.transAxes, fontsize=15, fontweight="bold")

    fig.tight_layout()
    save(fig, out, "fig3_agreement_beyond_spectrum")

    # the two panels the manuscript uses (Supplementary Figure S9)
    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4.8))
    panel_where(ax2[0], "A")
    panel_hotspot(ax2[1], "B")
    fig2.tight_layout()
    save(fig2, out, "fig3_agreement_position_hotspot")

    # site sharing between platforms: same mouse against different mice
    J = np.array([[jaccard({(p, a) for p, _r, a in sv[x]}, {(p, a) for p, _r, a in ov[y]}) for y in WT]
                  for x in WT])
    st["beyond"] = {
        "spectrum_r_matrix": R.round(3).tolist(),
        "spectrum_r_matched_mean": float(np.mean(np.diag(R))),
        "spectrum_r_mismatched_mean": float(R[~np.eye(3, dtype=bool)].mean()),
        "positional_r": float(rp[0]), "positional_p": float(rp[1]),
        "hotspot": {plat: H[plat].round(1).tolist() for plat in H}, "hotspot_expected": exp,
        "load_x1e3": {"sanger": dict(zip(WT, fs_, strict=True)), "ont": dict(zip(WT, fo_, strict=True)),
                      "ont_clone_weighted": dict(zip(WT, fc_, strict=True))},
        "site_jaccard": J.round(3).tolist(),
        "site_jaccard_matched_mean": float(np.mean(np.diag(J))),
        "site_jaccard_mismatched_mean": float(J[~np.eye(3, dtype=bool)].mean()),
    }


def fig_per_mouse(D, out):
    """Supplementary S2/S3: one matrix per mouse and platform."""
    fig, axes = plt.subplots(2, 6, figsize=(21, 7.4))
    for row, (plat, kos, lab) in enumerate((("sanger", KO_SANGER, "Sanger"), ("ont", KO_ONT, "ONT"))):
        ms = {n: matrix(D[(plat, n)]["variants"]) for n in WT + kos}
        # shading scaled within wild type and within knockouts, as in the
        # manuscript: one knockout call is 100% of its sample
        vm = {g: max([100 * ms[n].max() / ms[n].sum() for n in names if ms[n].sum()] or [1])
              for g, names in (("wt", WT), ("ko", kos))}
        for k, n in enumerate(WT + kos):
            m = ms[n]
            heat(axes[row, k], m, f"{n} ({lab})  n = {int(m.sum())}",
                 vm["wt" if n in WT else "ko"], show_from=(k == 0))
    fig.tight_layout()
    save(fig, out, "figS_per_mouse_matrices")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ont-mode", default="recommended")
    ap.add_argument("--suffix", default="", help="appended to the output folder name")
    # one folder per mouse, built from the AB1 traces by organize_jh4_ab1.py
    # (SANGER_VS_NANOPORE.md, "Which clone belongs to which mouse")
    ap.add_argument("--sanger-dir", default=str(Path.home() / "Nanoshm/sanger_ab1_by_mouse"))
    ap.add_argument("--ont-dir", default=None,
                    help="the deposit's Original_FASTQ folder (delivered reads; the manuscript's input, SHMtrace "
                         "drawing 500 usable reads) or Processed_FASTQ_Depth_Subsets (the old *_sub500 files)")
    ap.add_argument("--out", default=str(ROOT / "docs/calibration/figures/sanger_vs_nanopore"))
    a = ap.parse_args()
    out = Path(a.out) / (a.ont_mode + a.suffix)
    out.mkdir(parents=True, exist_ok=True)
    ref = read_fasta(REF)[0][1].upper()
    sd = Path(a.sanger_dir).expanduser()
    # knockout mice with no Sanger clones (A3 has no usable trace) are left out
    KO_SANGER[:] = [n for n in KO_SANGER if sanger_path(sd, n)]
    D = run_platforms(sd, a.ont_mode, Path(a.ont_dir).expanduser() if a.ont_dir else None)
    st = {"ont_mode": a.ont_mode,
          "unique_mutations": {f"{p}:{n}": len(d["variants"]) for (p, n), d in D.items()},
          "sequences": {f"{p}:{n}": d["n"] for (p, n), d in D.items()},
          "warnings": {f"{p}:{n}": d["warnings"] for (p, n), d in D.items() if d["warnings"]}}
    fig_density(D, out, st)
    fig_spectrum(D, out, st, ref)
    fig_beyond(D, out, st, ref)
    fig_per_mouse(D, out)
    with open(out / "stats.json", "w") as fh:
        json.dump(st, fh, indent=1, default=float)
    print(json.dumps(st, indent=1, default=float))


if __name__ == "__main__":
    main()
