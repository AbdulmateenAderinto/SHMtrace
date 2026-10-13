#!/usr/bin/env python3
"""Generate the Ver13 graphical abstract from audited manuscript values.

The artwork is deterministic Matplotlib vector drawing. It deliberately separates
the two biological validation systems, the shared ONT workflow, and the results
that are supported directly by the manuscript. Run without arguments to replace
the PDF, PNG, and TIFF files in ``5_Graphical_Abstract``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


# Pooled Wt spectra (% of calls, A>T ... G>C) from the analysis behind Figure 3: Sanger
# from the chromatograms, ONT from SHMtrace (sanger_vs_nanopore.py, recommended_manuscript)
import json as _json
_ST = _json.load(open(Path(__file__).resolve().parents[2] /
                      "figures/sanger_vs_nanopore/recommended_manuscript/stats.json"))["spectrum"]["pooled_props"]
_ORDER = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
JH4_SANGER = np.array([_ST["sanger"][k] for k in _ORDER])
JH4_ONT = np.array([_ST["ont"][k] for k in _ORDER])
R_SPECTRUM = float(np.corrcoef(JH4_SANGER, JH4_ONT)[0, 1])

# Headline numbers from the SHMtrace results (never typed in): fold differences with a library
# without calls at its one-event density (a lower bound), recovery at 500 reads, control background
_CAL = Path(__file__).resolve().parents[2]


def _geo(x):
    return float(10 ** np.mean(np.log10(x)))


def _fold(run, case, ctrl):
    S = {s["sample"]: s for s in _json.load(open(_CAL / f"figures/{run}/project.json"))["samples"]}
    d = lambda ks: [S[k]["unique_density_x1e5"] or S[k]["one_event_density_x1e5"] for k in ks]  # noqa: E731
    censored = any(not S[k]["unique_density_x1e5"] for k in ctrl)
    return _geo(d(case)) / _geo(d(ctrl)), censored, _geo(d(ctrl))


def _recovery(assay, case):
    dd = _json.load(open(_CAL / f"figures/depth_shmtrace/{assay}.json"))
    return float(np.mean([dd["500"][k]["calls"] / max(dd[x][k]["calls"] for x in dd) for k in case])) * 100


JH4_FOLD, JH4_CENS, JH4_BG = _fold("jh4_shmtrace", ["WT2", "WT3", "WT4"], ["KO1", "KO2", "KO3"])
N114_FOLD, N114_CENS, _ = _fold("n114_shmtrace", ["C1", "D3", "D5"], ["C3", "C4", "D1"])
JH4_REC = _recovery("jh4", ["WT2", "WT3", "WT4"])
N114_REC = _recovery("n114", ["C1", "D3", "D5"])

# Okabe-Ito-derived, colour-blind-safe palette.
INK = "#17212B"
MUTED = "#687684"
LINE = "#D7E0E8"
PANEL = "#F6F8FA"
BLUE = "#0072B2"
BLUE_D = "#174E78"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#D55E00"
AMBER = "#E69F00"
PURPLE = "#CC79A7"
GREY = "#9AA7B2"
PALE_BLUE = "#E9F3FA"
PALE_GREEN = "#E8F6F1"
PALE_ORANGE = "#FFF2E9"
PALE_PURPLE = "#F7ECF3"


mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "mathtext.fontset": "dejavusans",
})


def rounded(ax, x, y, w, h, *, fc="white", ec=LINE, lw=1.2, radius=0.8, z=1):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.18,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, x1, y1, x2, y2, *, color=MUTED, lw=1.6, z=4):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
        linewidth=lw, color=color, shrinkA=0, shrinkB=0, zorder=z,
    ))


def section_heading(ax, x, y, number, text, color):
    ax.add_patch(Circle((x, y), 0.78, facecolor=color, edgecolor="none", zorder=5))
    ax.text(x, y, str(number), ha="center", va="center", color="white",
            fontsize=10.5, fontweight="bold", zorder=6)
    ax.text(x + 1.35, y, text, ha="left", va="center", fontsize=14,
            fontweight="bold", color=INK, zorder=6)


def cell_cluster(ax, cx, cy, *, color, scale=1.0):
    for dx, dy in [(-0.9, 0.0), (0.0, 0.55), (0.85, 0.0), (0.0, -0.7)]:
        ax.add_patch(Circle((cx + dx * scale, cy + dy * scale), 0.68 * scale,
                            facecolor="white", edgecolor=color, linewidth=1.7, zorder=4))
        ax.add_patch(Circle((cx + dx * scale, cy + dy * scale), 0.18 * scale,
                            facecolor=color, edgecolor="none", zorder=5))


def mouse_icon(ax, x, y, *, color=GREY):
    ax.add_patch(Circle((x, y), 0.85, facecolor="#E2E7EB", edgecolor=color, linewidth=1.2, zorder=4))
    ax.add_patch(Circle((x + 0.82, y + 0.26), 0.42, facecolor="#E2E7EB", edgecolor=color, linewidth=1.2, zorder=4))
    ax.add_patch(Circle((x + 0.77, y + 0.68), 0.20, facecolor="#E2E7EB", edgecolor=color, linewidth=1.0, zorder=4))
    ax.add_patch(Circle((x + 0.96, y + 0.34), 0.06, facecolor=INK, edgecolor="none", zorder=5))
    ax.plot([x - 0.75, x - 1.35, x - 1.75], [y - 0.15, y + 0.25, y - 0.05],
            color=color, linewidth=1.2, solid_capstyle="round", zorder=3)


def amplicon(ax, x, y, w, *, color, label, start_trim=0.0, end_trim=0.0):
    ax.plot([x, x + w], [y, y], color=MUTED, linewidth=1.0, zorder=2)
    ax.add_patch(Rectangle((x, y - 0.35), w, 0.70, facecolor="white",
                           edgecolor=color, linewidth=1.4, zorder=3))
    if start_trim:
        ax.add_patch(Rectangle((x, y - 0.35), w * start_trim, 0.70,
                               facecolor="#E2E7EB", edgecolor="none", zorder=4))
    if end_trim:
        ax.add_patch(Rectangle((x + w * (1 - end_trim), y - 0.35), w * end_trim, 0.70,
                               facecolor="#E2E7EB", edgecolor="none", zorder=4))
    ax.text(x + w / 2, y + 0.85, label, ha="center", va="bottom",
            fontsize=9.2, fontweight="bold", color=color, zorder=5)


def pipeline_box(
    ax, x, y, w, title, subtitle, color, fill, *,
    subtitle_fontsize=7.05, subtitle_y=1.28, subtitle_linespacing=1.08,
):
    rounded(ax, x, y, w, 5.1, fc=fill, ec=color, lw=1.1, radius=0.55, z=2)
    ax.text(x + w / 2, y + 3.45, title, ha="center", va="center",
            fontsize=8.9, fontweight="bold", color=color, zorder=4)
    ax.text(x + w / 2, y + subtitle_y, subtitle, ha="center", va="center",
            fontsize=subtitle_fontsize, color=INK,
            linespacing=subtitle_linespacing, zorder=4)


def result_card(ax, x, y, w, h, *, title, color, fill):
    rounded(ax, x, y, w, h, fc=fill, ec=color, lw=1.1, radius=0.55, z=2)
    ax.text(x + 0.8, y + h - 1.1, title, ha="left", va="center",
            fontsize=10.3, fontweight="bold", color=color, zorder=5)


def main(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    # Elsevier minimum 1328 x 531 px. This renders at 3984 x 1594 px, 300 dpi.
    fig = plt.figure(figsize=(13.28, 5.313), dpi=300, facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    ax.axis("off")

    ax.text(2.1, 38.4, "Cloning-free long-read profiling of somatic hypermutation",
            fontsize=17.5, fontweight="bold", ha="left", va="center", color=INK)
    ax.text(2.1, 36.6, "Two biological systems test one control-based Oxford Nanopore workflow",
            fontsize=10.2, ha="left", va="center", color=MUTED)

    # Three primary columns.
    rounded(ax, 1.6, 3.4, 31.0, 31.8, fc=PANEL, ec=LINE, lw=1.25, radius=0.8)
    rounded(ax, 34.2, 3.4, 31.2, 31.8, fc=PANEL, ec=LINE, lw=1.25, radius=0.8)
    rounded(ax, 67.0, 3.4, 31.4, 31.8, fc=PANEL, ec=LINE, lw=1.25, radius=0.8)
    section_heading(ax, 3.3, 33.2, 1, "Complementary test systems", BLUE_D)
    section_heading(ax, 35.9, 33.2, 2, "Shared ONT workflow", GREEN)
    section_heading(ax, 68.7, 33.2, 3, "Validated outcomes", ORANGE)

    # ------------------------------------------------------------------ left
    rounded(ax, 3.0, 19.4, 28.2, 11.7, fc=PALE_BLUE, ec="#B9D7EA", lw=1.0, radius=0.65)
    ax.text(4.1, 29.7, r"Mouse GC B cells: $J_{\mathrm{H}4}$ downstream intron",
            fontsize=10.6, fontweight="bold", color=BLUE_D, ha="left", va="center")
    mouse_icon(ax, 5.3, 25.0)
    cell_cluster(ax, 10.0, 25.0, color=BLUE, scale=0.85)
    arrow(ax, 6.8, 25.0, 8.0, 25.0, color=BLUE_D, lw=1.3)
    amplicon(ax, 13.2, 25.0, 7.3, color=BLUE_D, label="~1.23-kb amplicon")
    arrow(ax, 11.7, 25.0, 12.6, 25.0, color=BLUE_D, lw=1.3)
    rounded(ax, 22.0, 26.0, 7.9, 2.2, fc="white", ec=BLUE, lw=0.9, radius=0.35)
    rounded(ax, 22.0, 22.7, 7.9, 2.2, fc="white", ec=ORANGE, lw=0.9, radius=0.35)
    ax.text(25.95, 27.1, r"Wt  ($n=3$)", ha="center", va="center", fontsize=8.8,
            fontweight="bold", color=BLUE_D)
    ax.text(25.95, 23.8, r"$\it{Aicda}^{-/-}$  ($n=3$)", ha="center", va="center",
            fontsize=8.6, fontweight="bold", color=ORANGE)
    ax.text(17.1, 20.6, "Same mice benchmarked by\ncloned-amplicon Sanger and ONT",
            fontsize=8.0, color=MUTED, ha="center", va="center", linespacing=1.2)

    rounded(ax, 3.0, 6.0, 28.2, 11.7, fc=PALE_GREEN, ec="#B5DECf", lw=1.0, radius=0.65)
    # Use the same left margin as the mouse-system heading, but separate the
    # system name from its descriptive subtitle for a cleaner two-line hierarchy.
    ax.text(4.1, 16.55, "N114 hybridoma",
            fontsize=10.0, fontweight="bold", color=GREEN, ha="left", va="center")
    ax.text(4.1, 15.45, "Fixed immunoglobulin V region",
            fontsize=8.5, fontweight="bold", color=GREEN, ha="left", va="center")
    cell_cluster(ax, 5.8, 11.2, color=GREEN, scale=0.9)
    ax.text(5.8, 8.2, "N114", ha="center", va="center", fontsize=8.4, color=MUTED)
    arrow(ax, 8.2, 11.2, 10.2, 11.2, color=GREEN, lw=1.3)
    rounded(ax, 9.7, 12.05, 9.2, 2.2, fc="white", ec=GREEN, lw=0.9, radius=0.35)
    rounded(ax, 9.7, 8.65, 9.2, 2.2, fc="white", ec=ORANGE, lw=0.9, radius=0.35)
    ax.text(14.3, 13.15, r"pMX-$\it{Aicda}$", ha="center", va="center", fontsize=8.4,
            fontweight="bold", color=GREEN)
    ax.text(14.3, 9.75, "pMX-Empty vector", ha="center", va="center", fontsize=8.1,
            fontweight="bold", color=ORANGE)
    arrow(ax, 19.2, 11.2, 20.3, 11.2, color=GREEN, lw=1.3)
    amplicon(ax, 20.8, 11.2, 8.2, color=GREEN, label="~0.65-kb amplicon",
             start_trim=20/655, end_trim=20/655)
    ax.text(19.4, 7.7, r"$n=3$ independent subclones per condition",
            ha="center", va="center", fontsize=7.9, color=MUTED)

    # ---------------------------------------------------------------- center
    ax.text(36.2, 29.5, "R10.4.1  |  Kit 14  |  SUP basecalling",
            fontsize=10.0, fontweight="bold", color=GREEN, ha="left", va="center")
    ax.text(36.2, 27.8, "Intact amplicons are retained as single-molecule reads",
            fontsize=8.8, color=MUTED, ha="left", va="center")
    for y in [25.6, 24.4, 23.2]:
        ax.add_patch(Rectangle((36.4, y), 26.8, 0.55, facecolor=GREEN,
                               edgecolor="none", alpha=0.86, zorder=3))
    for x in [44.0, 54.5, 60.6]:
        ax.add_patch(Rectangle((x, 23.0), 0.42, 3.45, facecolor=PURPLE,
                               edgecolor="none", zorder=4))
    ax.text(49.8, 21.6, "substitutions remain phased within each read",
            ha="center", va="center", fontsize=8.2, color=MUTED)

    # Give the information-dense final step extra width rather than forcing its
    # definition into the same narrow footprint as the three simpler steps.
    # The process row, calling-rule strip and interpretation card share the
    # same 36.0–63.5 alignment grid.
    box_widths = [6.25, 6.25, 5.90, 8.05]
    xs = [36.0, 42.60, 49.20, 55.45]
    # SHMtrace's own steps, in its order: quality gates, competitive alignment, depth, calling
    titles = ["Read QC", "Align", "Depth", "Call SNVs"]
    subtitles = [
        "length ≤ max\nEE ≤ 0.01/base",
        "vs. paralogue;\nchimeras\nremoved",
        "random\n500 reads",
        "tested against\ncontrol errors",
    ]
    colors = [GREEN, BLUE_D, AMBER, PURPLE]
    fills = [PALE_GREEN, PALE_BLUE, "#FFF7DD", PALE_PURPLE]
    for i, (x, box_w, title, subtitle, color, fill) in enumerate(
            zip(xs, box_widths, titles, subtitles, colors, fills)):
        pipeline_box(ax, x, 13.7, box_w, title, subtitle, color, fill)
        if i < len(xs) - 1:
            arrow(ax, x + box_w + 0.08, 16.72, xs[i+1] - 0.18, 16.72,
                  color=MUTED, lw=1.2)

    # Keep the detailed calling definition out of the narrow process box. This
    # dedicated two-line annotation is robust to PDF font rendering and scaling.
    rounded(ax, 36.0, 11.05, 27.5, 2.0, fc=PALE_PURPLE, ec="#E7BBD4",
            lw=0.75, radius=0.38)
    ax.text(49.75, 12.05,
            "SHMtrace: Benjamini–Hochberg FDR ≤ 0.05; sites at ≥ 50% removed\n"
            "Each position–substitution counted once per sample",
            ha="center", va="center", fontsize=6.45, color=PURPLE,
            fontweight="bold", linespacing=1.08)

    rounded(ax, 36.0, 5.9, 27.5, 4.65, fc="white", ec=LINE, lw=0.9, radius=0.55)
    ax.text(37.0, 9.05, "Control-based interpretation", fontsize=9.5,
            fontweight="bold", color=BLUE_D, ha="left", va="center")
    ax.text(37.0, 7.35,
            "$\\it{Aicda}^{-/-}$ or empty-vector controls\nmeasure residual background",
            fontsize=8.0, color=INK, ha="left", va="center", linespacing=1.12)

    # --------------------------------------------------------------- outcomes
    result_card(ax, 68.4, 20.0, 13.8, 10.8,
                title="Sanger concordance", color=BLUE_D, fill=PALE_BLUE)
    sax = fig.add_axes([0.696, 0.555, 0.097, 0.155])
    sax.plot([0, 20.5], [0, 20.5], linestyle="--", color=GREY, linewidth=1.0)
    sax.scatter(JH4_SANGER, JH4_ONT, s=23, facecolor=PURPLE,
                edgecolor="white", linewidth=0.55, zorder=3)
    sax.set_xlim(0, 20.5)
    sax.set_ylim(0, 20.5)
    sax.set_xticks([0, 10, 20])
    sax.set_yticks([0, 10, 20])
    sax.tick_params(labelsize=6.2, length=2, pad=1)
    sax.set_xlabel("Sanger (%)", fontsize=6.8, labelpad=0.5)
    sax.set_ylabel("ONT (%)", fontsize=6.8, labelpad=0.5)
    for spine in ("top", "right"):
        sax.spines[spine].set_visible(False)
    sax.text(19.3, 1.4, rf"$r={R_SPECTRUM:.2f}$", ha="right", va="bottom",
             fontsize=7.2, fontweight="bold", color=BLUE_D)

    result_card(ax, 83.1, 20.0, 13.8, 10.8,
                title="AID dependence", color=GREEN, fill=PALE_GREEN)
    ax.text(84.2, 27.4, r"Mouse $J_{\mathrm{H}4}$", fontsize=8.4,
            fontweight="bold", color=INK, ha="left", va="center")
    ax.text(95.8, 27.4, f"{'≥' if JH4_CENS else '~'}{JH4_FOLD:.0f}×", fontsize=12.5, fontweight="bold",
            color=BLUE, ha="right", va="center")
    ax.text(84.2, 24.7, "N114 V region", fontsize=8.4,
            fontweight="bold", color=INK, ha="left", va="center")
    ax.text(95.8, 24.7, f"{'≥' if N114_CENS else '~'}{N114_FOLD:.0f}×", fontsize=12.5, fontweight="bold",
            color=GREEN, ha="right", va="center")
    ax.text(84.2, 22.0, "G:C-biased AID signature", fontsize=8.2,
            color=MUTED, ha="left", va="center")

    result_card(ax, 68.4, 7.0, 13.8, 10.8,
                title="Longer region", color=PURPLE, fill=PALE_PURPLE)
    ax.text(69.4, 14.4, r"$J_{\mathrm{H}4}$: 1.23-kb intact reads", fontsize=8.2,
            color=INK, ha="left", va="center")
    # 1,141-bp reliable window with the 529-bp shared Sanger region overlaid.
    ax.add_patch(Rectangle((69.5, 11.0), 11.2, 0.95, facecolor=SKY,
                           edgecolor=BLUE_D, linewidth=0.8, zorder=3))
    ax.add_patch(Rectangle((69.8, 11.0), 5.2, 0.95, facecolor=PURPLE,
                           edgecolor="none", alpha=0.92, zorder=4))
    ax.text(75.1, 9.8, "1,141-bp reliable analysis", ha="center", va="center",
            fontsize=7.9, fontweight="bold", color=BLUE_D)
    ax.text(72.4, 12.2, "529-bp Sanger window", ha="center", va="bottom",
            fontsize=6.8, color=PURPLE)
    ax.text(75.1, 8.3, "long-range linkage retained",
            ha="center", va="center", fontsize=7.6, color=MUTED)

    result_card(ax, 83.1, 7.0, 13.8, 10.8,
                title="Depth + background", color=ORANGE, fill=PALE_ORANGE)
    rounded(ax, 84.0, 12.1, 5.8, 3.2, fc="white", ec=BLUE, lw=0.9, radius=0.35)
    rounded(ax, 90.3, 12.1, 5.8, 3.2, fc="white", ec=GREEN, lw=0.9, radius=0.35)
    ax.text(86.9, 13.8, f"{JH4_REC:.0f}%", fontsize=12.0, fontweight="bold", color=BLUE,
            ha="center", va="center")
    ax.text(86.9, 12.7, r"$J_{\mathrm{H}4}$", fontsize=7.5, color=MUTED,
            ha="center", va="center")
    ax.text(93.2, 13.8, f"{N114_REC:.0f}%", fontsize=12.0, fontweight="bold", color=GREEN,
            ha="center", va="center")
    ax.text(93.2, 12.7, "N114", fontsize=7.5, color=MUTED,
            ha="center", va="center")
    ax.text(90.0, 10.9, "variants recovered at 500 reads", fontsize=7.6,
            color=INK, ha="center", va="center")
    ax.text(90.0, 9.0, r"$Aicda^{-/-}$ control background", fontsize=7.3,
            color=MUTED, ha="center", va="center")
    ax.text(90.0, 8.0, rf"≈{JH4_BG:.2f} × 10$^{{-5}}$ variants/bp", fontsize=7.5,
            color=ORANGE, ha="center", va="center", fontweight="bold")

    # Cross-column arrows and final message.
    arrow(ax, 32.7, 19.2, 33.8, 19.2, color=BLUE_D, lw=1.7)
    arrow(ax, 65.5, 19.2, 66.6, 19.2, color=BLUE_D, lw=1.7)
    ax.text(50.0, 1.35,
            "Cloning-free  •  intact long molecules  •  depth-standardised, substitution-focused profiling",
            fontsize=10.3, fontweight="bold", ha="center", va="center", color=BLUE_D)

    # Deterministic visual-fit assertions for the canvas.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in ax.texts:
        bbox = text.get_window_extent(renderer=renderer)
        inv = ax.transData.inverted()
        (x0, y0), (x1, y1) = inv.transform([(bbox.x0, bbox.y0), (bbox.x1, bbox.y1)])
        if x0 < -0.2 or x1 > 100.2 or y0 < -0.2 or y1 > 40.2:
            raise RuntimeError(f"Text outside canvas: {text.get_text()!r} -> {(x0, y0, x1, y1)}")

    stem = output_dir / "ONT_SHM_GraphicalAbstract"
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".tif"), dpi=300, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    from PIL import Image

    for suffix in (".png", ".tif"):
        path = stem.with_suffix(suffix)
        with Image.open(path) as image:
            print(f"{path.name}: {image.width} x {image.height} px; dpi={image.info.get('dpi')}")
            if image.width < 1328 or image.height < 531:
                raise RuntimeError(f"Elsevier size minimum not met by {path.name}")
    print(f"Wrote {stem.with_suffix('.pdf').name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    default_output = Path(__file__).resolve().parents[1] / "figures" / "graphical_abstract"
    parser.add_argument("--out-dir", type=Path, default=default_output)
    args = parser.parse_args()
    main(args.out_dir)
