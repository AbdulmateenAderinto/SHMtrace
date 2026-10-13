#!/usr/bin/env python3
"""Publication-size Figure 1: experimental design and ONT read processing.

The canvas is 7.1 inches wide so text remains legible when placed at the
manuscript's full text width. Complete filter-attrition counts are not available
for every sample, so the figure does not infer or display them.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle


BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "figures"
PDF_OUT = BASE / "figures" / "pdf"
OUT.mkdir(exist_ok=True)
PDF_OUT.mkdir(parents=True, exist_ok=True)

# Colour-blind-safe status colours with one shared step colour.
NAVY = "#1F4E79"
BLUE = "#0072B2"
TEAL = "#009E73"
VERMILION = "#D55E00"
ORANGE = "#E69F00"
PURPLE = "#7B4FA3"
STEP_GREEN = "#4F9F56"
STEP_AMBER = "#C88424"
STEP_RUST = "#A9472B"
INK = "#202428"
MID = "#626C75"
LIGHT = "#DCE2E7"
PALE_BLUE = "#EDF4F8"
PALE_TEAL = "#EEF7F4"
PALE_ORANGE = "#FFF0E8"
PALE_YELLOW = "#FFF1C7"
PALE_PURPLE = "#EFE7F5"
WHITE = "#FFFFFF"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": WHITE,
        "savefig.edgecolor": "none",
    }
)

W, H = 7.10, 8.00
fig = plt.figure(figsize=(W, H), facecolor=WHITE)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")


def box(x, y, w, h, fc=WHITE, ec=LIGHT, lw=0.8, radius=0.07, z=1):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    )
    ax.add_patch(p)
    return p


def arrow(x1, y1, x2, y2, color=MID, lw=1.0, ms=8, connection="arc3"):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=ms,
            linewidth=lw, color=color, connectionstyle=connection,
            shrinkA=0, shrinkB=0, zorder=8,
        )
    )


def pill(x, y, w, h, text, fc, fontsize=6.7, weight="bold"):
    box(x, y, w, h, fc=fc, ec="none", radius=h / 2, z=5)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=weight, color=INK, zorder=6)


def status_mark(x, y, ok, r=0.055):
    col = TEAL if ok else VERMILION
    ax.add_patch(Circle((x, y), r, facecolor=col, edgecolor=WHITE, linewidth=0.4, zorder=9))
    if ok:
        ax.plot([x - 0.026, x - 0.006, x + 0.030],
                [y - 0.002, y - 0.022, y + 0.022],
                color=WHITE, linewidth=0.8, solid_capstyle="round", zorder=10)
    else:
        ax.plot([x - 0.021, x + 0.021], [y - 0.021, y + 0.021],
                color=WHITE, linewidth=0.8, solid_capstyle="round", zorder=10)
        ax.plot([x - 0.021, x + 0.021], [y + 0.021, y - 0.021],
                color=WHITE, linewidth=0.8, solid_capstyle="round", zorder=10)


def read_bar(x, y, w, color=TEAL, h=0.055, dashed=False, chimera=False, z=4):
    if chimera:
        box(x, y, w / 2, h, fc=TEAL, ec="none", radius=h / 2, z=z)
        box(x + w / 2, y, w / 2, h, fc=VERMILION, ec="none", radius=h / 2, z=z)
    elif dashed:
        p = box(x, y, w, h, fc=WHITE, ec=MID, lw=0.55, radius=h / 2, z=z)
        p.set_linestyle((0, (2.2, 1.6)))
    else:
        box(x, y, w, h, fc=color, ec="none", radius=h / 2, z=z)


def workflow_legend(y=4.345):
    """Compact read-state key placed directly above the numbered workflow."""
    items = [
        (0.34, "target read", "target"),
        (1.56, "contaminant", "contaminant"),
        (2.88, "chimera", "chimera"),
        (4.00, "rejected / wrong length", "dashed"),
        (5.78, "error marks", "errors"),
    ]
    for x, label, kind in items:
        if kind == "target":
            read_bar(x, y - 0.018, 0.19, TEAL, h=0.037)
        elif kind == "contaminant":
            read_bar(x, y - 0.018, 0.19, VERMILION, h=0.037)
        elif kind == "chimera":
            read_bar(x, y - 0.018, 0.19, chimera=True, h=0.037)
        elif kind == "dashed":
            read_bar(x, y - 0.018, 0.24, dashed=True, h=0.037)
        else:
            read_bar(x, y - 0.018, 0.19, TEAL, h=0.037)
            for dx in (0.055, 0.105, 0.155):
                ax.add_patch(Circle((x + dx, y), 0.012,
                                    facecolor=ORANGE, edgecolor="none", zorder=7))
        offset = 0.30 if kind == "dashed" else 0.25
        ax.text(x + offset, y, label, fontsize=4.55, color=MID,
                ha="left", va="center")


def mouse(x, y, s=1.0):
    ax.add_patch(Circle((x, y), 0.105 * s, facecolor=LIGHT, edgecolor=MID, linewidth=0.65, zorder=4))
    ax.add_patch(Circle((x + 0.115 * s, y + 0.030 * s), 0.070 * s, facecolor=LIGHT, edgecolor=MID, linewidth=0.65, zorder=4))
    ax.add_patch(Circle((x + 0.108 * s, y + 0.098 * s), 0.033 * s, facecolor=LIGHT, edgecolor=MID, linewidth=0.5, zorder=4))
    ax.add_patch(Circle((x + 0.145 * s, y + 0.043 * s), 0.010 * s, facecolor=INK, edgecolor="none", zorder=5))
    ax.plot([x - 0.09 * s, x - 0.18 * s, x - 0.23 * s],
            [y, y + 0.07 * s, y + 0.03 * s], color=MID, linewidth=0.65, zorder=3)


def cells(x, y, color, s=1.0):
    for dx, dy in [(0, 0), (0.12, 0.03), (-0.10, 0.05), (0.04, 0.13), (-0.05, -0.09)]:
        ax.add_patch(Circle((x + dx * s, y + dy * s), 0.075 * s,
                            facecolor=WHITE, edgecolor=color, linewidth=0.8, zorder=4))
        ax.add_patch(Circle((x + dx * s, y + dy * s), 0.022 * s,
                            facecolor=color, edgecolor="none", zorder=5))


def amplicon(x, y, w, label):
    ax.plot([x, x + w], [y, y], color=MID, linewidth=0.8, zorder=3)
    ax.add_patch(Rectangle((x + 0.12 * w, y - 0.045), 0.76 * w, 0.09,
                           facecolor=PALE_BLUE, edgecolor=NAVY, linewidth=0.65, zorder=4))
    ax.add_patch(Polygon([[x + 0.12 * w, y], [x + 0.18 * w, y + 0.04], [x + 0.18 * w, y - 0.04]],
                         closed=True, facecolor=NAVY, edgecolor="none", zorder=5))
    ax.add_patch(Polygon([[x + 0.88 * w, y], [x + 0.82 * w, y + 0.04], [x + 0.82 * w, y - 0.04]],
                         closed=True, facecolor=NAVY, edgecolor="none", zorder=5))
    ax.text(x + w / 2, y + 0.15, label, ha="center", va="bottom", fontsize=6.8,
            fontweight="bold", color=INK)


def sequencer(x, y, s=1.0):
    box(x - 0.13 * s, y - 0.12 * s, 0.26 * s, 0.24 * s,
        fc=WHITE, ec=NAVY, lw=0.75, radius=0.025, z=4)
    ax.add_patch(Rectangle((x - 0.085 * s, y + 0.015 * s), 0.17 * s, 0.055 * s,
                           facecolor=PALE_BLUE, edgecolor=NAVY, linewidth=0.45, zorder=5))
    for i in range(3):
        ax.add_patch(Circle((x - 0.05 * s + i * 0.05 * s, y - 0.06 * s),
                            0.009 * s, facecolor=TEAL, edgecolor="none", zorder=5))


def sanger(x, y, s=1.0):
    ax.plot([x - 0.12 * s, x - 0.06 * s, x, x + 0.06 * s, x + 0.12 * s],
            [y - 0.02 * s, y + 0.07 * s, y - 0.06 * s, y + 0.05 * s, y - 0.01 * s],
            color=PURPLE, linewidth=1.0, zorder=4)


def gfp_high_gate(x, y, s=1.0):
    """Compact flow-cytometry plot with the GFP-high population gated."""
    w, h = 0.36 * s, 0.28 * s
    left, bottom = x - w / 2, y - h / 2
    ax.plot([left, left, left + w], [bottom + h, bottom, bottom],
            color=MID, linewidth=0.55, zorder=4)

    # Fixed point coordinates keep the schematic reproducible. The muted
    # lower-left cloud represents GFP-low/background cells; the teal upper-right
    # cloud is the population selected by the rectangular gate.
    low_points = [
        (0.08, 0.10), (0.14, 0.16), (0.19, 0.12), (0.22, 0.21),
        (0.28, 0.18), (0.32, 0.27), (0.37, 0.24), (0.42, 0.33),
        (0.47, 0.28), (0.52, 0.39), (0.34, 0.38), (0.24, 0.32),
    ]
    high_points = [
        (0.66, 0.62), (0.72, 0.70), (0.78, 0.66), (0.82, 0.78),
        (0.88, 0.73), (0.90, 0.86), (0.75, 0.84), (0.64, 0.76),
    ]
    for px, py in low_points:
        ax.add_patch(Circle((left + px * w, bottom + py * h), 0.010 * s,
                            facecolor="#AFC0CA", edgecolor="none", alpha=0.82, zorder=5))
    for px, py in high_points:
        ax.add_patch(Circle((left + px * w, bottom + py * h), 0.012 * s,
                            facecolor="#B7E36B", edgecolor=TEAL,
                            linewidth=0.25, zorder=6))
    ax.add_patch(Rectangle((left + 0.57 * w, bottom + 0.52 * h),
                           0.39 * w, 0.43 * h, fill=False,
                           edgecolor=TEAL, linewidth=0.75, zorder=7))
    ax.text(left + w, bottom + h + 0.025 * s, "GFP-high gate",
            fontsize=3.8 * s, fontweight="bold", color=TEAL,
            ha="right", va="bottom", zorder=7)


def downsampling_funnel(x, y, s=1.0):
    """A small funnel that explicitly encodes random depth reduction."""
    ax.add_patch(Polygon([
        (x - 0.105 * s, y + 0.095 * s),
        (x + 0.105 * s, y + 0.095 * s),
        (x + 0.035 * s, y - 0.005 * s),
        (x + 0.035 * s, y - 0.105 * s),
        (x - 0.035 * s, y - 0.105 * s),
        (x - 0.035 * s, y - 0.005 * s),
    ], closed=True, facecolor=PALE_BLUE, edgecolor=NAVY,
       linewidth=0.70, zorder=6))


def plate_96(x, y, s=1.0):
    """Small 8 × 12 well-plate icon."""
    w, h = 0.32 * s, 0.22 * s
    box(x - w / 2, y - h / 2, w, h, fc=WHITE, ec=NAVY,
        lw=0.65, radius=0.025, z=4)
    for row in range(4):
        for col in range(6):
            ax.add_patch(Circle((x - 0.125 * s + col * 0.050 * s,
                                 y - 0.070 * s + row * 0.047 * s),
                                0.010 * s, facecolor=PALE_BLUE,
                                edgecolor=NAVY, linewidth=0.30, zorder=5))


def culture_flask(x, y, s=1.0):
    """Culture flask used for puromycin maintenance."""
    ax.add_patch(Polygon([
        (x - 0.12 * s, y - 0.09 * s),
        (x + 0.12 * s, y - 0.09 * s),
        (x + 0.08 * s, y + 0.07 * s),
        (x - 0.08 * s, y + 0.07 * s),
    ], closed=True, facecolor=WHITE, edgecolor=NAVY, linewidth=0.65, zorder=4))
    ax.add_patch(Rectangle((x - 0.035 * s, y + 0.07 * s), 0.07 * s, 0.055 * s,
                           facecolor=WHITE, edgecolor=NAVY, linewidth=0.60, zorder=4))
    ax.add_patch(Rectangle((x - 0.105 * s, y - 0.075 * s), 0.21 * s, 0.050 * s,
                           facecolor="#DDF2EB", edgecolor="none", zorder=5))
    ax.text(x, y - 0.051 * s, "Puro", ha="center", va="center",
            fontsize=4.8 * s, fontweight="bold", color=TEAL, zorder=6)


def puromycin_selection(x, y, s=1.0):
    """Show positive selection: resistant cells expand; sensitive cells die."""
    culture_flask(x - 0.27 * s, y, 0.72 * s)

    # A short fork makes clear that puromycin exposure has two cell fates.
    branch_x = x - 0.10 * s
    ax.plot([x - 0.17 * s, branch_x], [y, y], color=MID,
            linewidth=0.65, zorder=4)
    arrow(branch_x, y, x + 0.07 * s, y + 0.14 * s,
          color=TEAL, lw=0.65, ms=5)
    arrow(branch_x, y, x + 0.07 * s, y - 0.14 * s,
          color=MID, lw=0.65, ms=5)

    # Puro-resistant cells: several healthy cells indicate survival and
    # expansion, rather than merely persistence.
    for dx, dy in [(0.11, 0.14), (0.18, 0.18),
                   (0.19, 0.11), (0.26, 0.16)]:
        ax.add_patch(Circle((x + dx * s, y + dy * s), 0.035 * s,
                            facecolor="#B7E36B", edgecolor=TEAL,
                            linewidth=0.50, zorder=5))
        ax.add_patch(Circle((x + dx * s, y + dy * s), 0.009 * s,
                            facecolor=TEAL, edgecolor="none", zorder=6))
    ax.text(x + 0.19 * s, y + 0.32 * s,
            "Puro-resistant\ncells expand", fontsize=4.4 * s,
            fontweight="bold", color=TEAL, ha="center", va="top",
            linespacing=0.92, zorder=7)

    # Cells lacking resistance are shown as a faded cell with an explicit
    # vermilion death mark. The paired label means colour is not the sole cue.
    dead_x, dead_y = x + 0.18 * s, y - 0.15 * s
    ax.add_patch(Circle((dead_x, dead_y), 0.041 * s,
                        facecolor=LIGHT, edgecolor=MID,
                        linewidth=0.50, alpha=0.78, zorder=5))
    ax.add_patch(Circle((dead_x, dead_y), 0.010 * s,
                        facecolor=MID, edgecolor="none", alpha=0.75, zorder=6))
    ax.plot([dead_x - 0.032 * s, dead_x + 0.032 * s],
            [dead_y - 0.032 * s, dead_y + 0.032 * s],
            color=VERMILION, linewidth=0.85, zorder=7)
    ax.plot([dead_x - 0.032 * s, dead_x + 0.032 * s],
            [dead_y + 0.032 * s, dead_y - 0.032 * s],
            color=VERMILION, linewidth=0.85, zorder=7)
    ax.text(x + 0.28 * s, y - 0.38 * s,
            "Puro-sensitive\ncells die", fontsize=4.4 * s,
            fontweight="bold", color=VERMILION, ha="center", va="top",
            linespacing=0.92, zorder=7)


# ---------------------------------------------------------------------------
# A - Experimental design. Stacked lanes preserve legibility at 7.1-inch width.
# ---------------------------------------------------------------------------
ax.text(0.24, 7.85, "A", fontsize=12.5, fontweight="bold", color=INK, va="top")
ax.text(0.52, 7.83, "Experimental design", fontsize=10.2, fontweight="bold", color=INK, va="top")
ax.plot([0.24, 6.87], [7.61, 7.61], color=LIGHT, linewidth=0.7)

# JH4 lane
# The lane box is extended 0.12 in downward (6.56 -> 6.44) purely to open a
# footer strip for the analysis intervals. No other element moves; the lowest
# existing label ("JH4 downstream intron") descends only to about y = 6.62.
box(0.27, 6.44, 6.58, 0.97, fc=PALE_BLUE, ec="none", radius=0.08, z=0)
ax.text(0.47, 7.17, r"$\mathbfit{J_{H}4}$ intron", fontsize=8.3, fontweight="bold", color=NAVY, va="center")
ax.text(0.47, 6.84, "4-month-old\nC57BL/6 mice", fontsize=6.4, color=INK, ha="left", va="center")
mouse(1.28, 6.96, 0.95)
pill(1.51, 7.08, 0.73, 0.20, "Wt  n = 3", "#DDECF4", fontsize=6.2)
pill(1.51, 6.74, 0.73, 0.20,
     r"$\mathbfit{Aicda}^{\mathbf{-/-}}$  n = 3",
     PALE_ORANGE, fontsize=5.8)
arrow(2.34, 6.98, 2.62, 6.98, ms=7)
cells(2.88, 6.99, BLUE, 0.92)
ax.text(2.88, 6.78, "Peyer's patch GC B cells\nB220$^{+}$ GL7$^{+}$ Fas$^{+}$", fontsize=6.1,
        color=INK, ha="center", va="top", linespacing=1.03)
arrow(3.23, 6.98, 3.51, 6.98, ms=7)
amplicon(3.56, 6.97, 0.96, "~1.23 kb")
ax.text(4.04, 6.70, r"$J_{H}4$ downstream intron", fontsize=5.9, color=INK, ha="center", va="top")
arrow(4.60, 6.98, 4.86, 6.98, ms=7)
sequencer(5.10, 7.13, 0.90)
sanger(5.10, 6.80, 0.90)
ax.text(5.31, 7.13, "ONT · R10.4.1 / Kit 14 SUP\n~500 reads per sample", fontsize=6.0, color=INK, ha="left", va="center")
ax.text(5.31, 6.80, "Sanger · cloned amplicons", fontsize=6.0, color=INK, ha="left", va="center")
# Footer strip. The two platforms share the animals but not the amplification,
# which is what makes the spectrum concordance in Figure 3 a genuine
# cross-platform result rather than a shared-PCR artefact. The two intervals are
# the denominators every density in the paper is computed over.
ax.text(3.56, 6.50,
        "Same mice, separate amplicon reactions · SuperFi II high-fidelity PCR"
        " · ONT reliable interval 1,141 bp · shared Sanger window 529 bp",
        fontsize=4.8, color=MID, ha="center", va="center")

# N114 lane - explicit selection, limiting-dilution cloning and maintenance.
box(0.27, 5.02, 6.58, 1.30, fc=PALE_TEAL, ec="none", radius=0.08, z=0)
ax.text(0.46, 6.17, "N114 V region", fontsize=8.1, fontweight="bold", color="#14735C", va="center")
cells(0.78, 5.70, TEAL, 0.84)
ax.text(0.78, 5.34, "N114\nhybridoma", fontsize=5.9, color=INK, ha="center", va="top")
pill(1.19, 5.86, 1.14, 0.20,
     r"pMX-$\mathbfit{Aicda}$ · GFP · Puro$^{\mathbf{R}}$",
     "#DDF2EB", fontsize=4.75)
pill(1.19, 5.51, 1.14, 0.20,
     r"pMX-Empty vector · GFP · Puro$^{\mathbf{R}}$",
     "#F2E8E4", fontsize=4.25)
ax.text(1.76, 5.27, "retroviral constructs", fontsize=5.5, color=MID, ha="center", va="top")

arrow(2.36, 5.70, 2.52, 5.70, ms=7)
gfp_high_gate(2.70, 5.73, 0.94)
ax.text(2.70, 5.37, "GFP-high\nselection", fontsize=5.7, fontweight="bold", color=INK, ha="center", va="top")

arrow(2.92, 5.70, 3.18, 5.70, ms=7)
plate_96(3.42, 5.72, 1.05)
ax.text(3.42, 6.12, "limiting dilution\n96-well subcloning", fontsize=5.5, fontweight="bold", color=INK, ha="center", va="top")

# Stop the incoming connector at the flask boundary so it does not cross the
# vessel.
arrow(3.63, 5.70, 3.79, 5.70, ms=5)
puromycin_selection(4.16, 5.72, 1.00)

# Only the resistant, expanding branch proceeds to amplicon analysis.
arrow(4.48, 5.86, 4.63, 5.72, color=TEAL, lw=0.85, ms=6)
amplicon(4.66, 5.70, 0.61, "~0.65 kb")
ax.text(4.97, 5.35, "fixed V region", fontsize=5.3, color=INK,
        ha="center", va="top")

arrow(5.30, 5.70, 5.39, 5.70, ms=7)
sequencer(5.59, 5.72, 0.82)
ax.text(6.30, 5.70,
        "ONT · R10.4.1\nKit 14 · SUP\nn = 3; ~500 reads/subclone",
        fontsize=5.5, color=INK, ha="center", va="center", linespacing=1.25)
# Footer strip. ~0.65 kb is the amplicon; 615 bp is the analysed denominator.
ax.text(3.56, 5.10,
        "SuperFi II high-fidelity PCR · analysed positions 21–635 of the"
        " 655 bp reference (615 bp)",
        fontsize=4.8, color=MID, ha="center", va="center")


# ---------------------------------------------------------------------------
# B - Five-step workflow. Three upper cards + two larger lower cards.
# ---------------------------------------------------------------------------
ax.text(0.24, 4.66, "B", fontsize=12.5, fontweight="bold", color=INK, va="top")
ax.text(0.52, 4.64, "Read filtering and variant calling", fontsize=10.2, fontweight="bold", color=INK, va="top")
ax.plot([0.24, 6.87], [4.42, 4.42], color=LIGHT, linewidth=0.7)
workflow_legend()


def stage(x, y, w, h, number, title, accent):
    box(x, y, w, h, fc=WHITE, ec=LIGHT, lw=0.75, radius=0.07, z=0)
    box(x + 0.008, y + h - 0.30, w - 0.016, 0.292,
        fc=accent, ec="none", radius=0.055, z=2)
    ax.add_patch(Circle((x + 0.19, y + h - 0.15), 0.105,
                        facecolor=WHITE, edgecolor="none", zorder=5))
    ax.text(x + 0.19, y + h - 0.15, str(number), ha="center", va="center",
            fontsize=7.0, fontweight="bold", color=accent, zorder=6)
    ax.text(x + 0.38, y + h - 0.085, title, ha="left", va="top",
            fontsize=6.9, fontweight="bold", color=WHITE,
            linespacing=0.95, zorder=6)


top_y, top_h, top_w = 2.55, 1.74, 2.03
top_x = [0.28, 2.54, 4.80]
stage(top_x[0], top_y, top_w, top_h, 1, "Raw reads → competitive\nalignment", NAVY)
stage(top_x[1], top_y, top_w, top_h, 2, "Amplicon-length filter", TEAL)
stage(top_x[2], top_y, top_w, top_h, 3, "Expected-error filter", STEP_GREEN)

bottom_y, bottom_h = 0.27, 1.74
# Centre the combined lower row under the full three-card upper row. The lower
# pair spans x=0.79–6.32, giving nearly identical left and right margins beneath
# the x=0.28–6.83 upper row.
bottom_x4, bottom_x5 = 0.79, 3.50
stage(bottom_x4, bottom_y, 2.50, bottom_h, 4, "Depth standardisation", STEP_AMBER)
stage(bottom_x5, bottom_y, 2.82, bottom_h, 5, "Reference alignment and SNV calling", STEP_RUST)

arrow(2.33, 3.42, 2.51, 3.42, ms=8)
arrow(4.59, 3.42, 4.77, 3.42, ms=8)
# Collision-free elbow from step 3 to step 4.
ax.plot([5.82, 5.82], [2.53, 2.28], color=MID, linewidth=1.0, zorder=7)
ax.plot([5.82, 2.04], [2.28, 2.28], color=MID, linewidth=1.0, zorder=7)
arrow(2.04, 2.28, 2.04, 2.04, ms=8)
arrow(3.31, 1.14, 3.47, 1.14, ms=8)

# Step 1: competitive reference matching. The card is arranged as a compact
# decision table so each input read can be compared across both references
# before the retain/exclude outcome is read at the far right.
x, y = top_x[0], top_y

# Column header band and aligned row guides make the comparison direction
# explicit without adding a separate legend.
box(x + 0.10, 3.70, 1.83, 0.22, fc="#F3F6F8", ec="none", radius=0.035, z=1)
ax.text(x + 0.36, 3.81, "input read", fontsize=4.4, fontweight="bold",
        color=MID, ha="center", va="center")
ax.text(x + 1.02, 3.81, "target ref", fontsize=4.4, fontweight="bold",
        color=INK, ha="center", va="center")
ax.text(x + 1.43, 3.81, "contam. ref", fontsize=4.4, fontweight="bold",
        color=INK, ha="center", va="center")
ax.text(x + 1.79, 3.81, "retain?", fontsize=4.2, fontweight="bold",
        color=MID, ha="center", va="center")
for yy in (3.35, 3.10):
    ax.plot([x + 0.10, x + 1.93], [yy, yy], color=LIGHT,
            linewidth=0.38, zorder=1)

def weak_alignment(bx, by):
    read_bar(bx, by, 0.34, dashed=True, h=0.045)
    for dx in (0.09, 0.17, 0.25):
        ax.add_patch(Circle((bx + dx, by + 0.022), 0.012,
                            facecolor=ORANGE, edgecolor="none", zorder=7))

# Intended read: strongest match to the target reference, therefore retained.
read_bar(x + 0.11, 3.49, 0.15, TEAL, h=0.045)
ax.text(x + 0.30, 3.51, "target", fontsize=4.8, color=INK,
        ha="left", va="center")
read_bar(x + 0.87, 3.49, 0.30, TEAL, h=0.055)
weak_alignment(x + 1.28, 3.49)
status_mark(x + 1.79, 3.52, True, r=0.048)

# Contaminant read: strongest match to its own reference, therefore excluded.
read_bar(x + 0.11, 3.22, 0.15, VERMILION, h=0.045)
ax.text(x + 0.30, 3.24, "contaminant", fontsize=4.5, color=INK,
        ha="left", va="center")
weak_alignment(x + 0.87, 3.22)
read_bar(x + 1.28, 3.22, 0.34, VERMILION, h=0.055)
status_mark(x + 1.79, 3.25, False, r=0.048)

# Split support across the two references identifies a chimera.
read_bar(x + 0.11, 2.85, 0.15, chimera=True, h=0.045)
ax.text(x + 0.30, 2.87, "chimera", fontsize=4.8, color=INK,
        ha="left", va="center")
read_bar(x + 0.87, 2.85, 0.22, TEAL, h=0.050)
read_bar(x + 1.28, 2.85, 0.22, VERMILION, h=0.050)
ax.text(x + 1.25, 2.955, "split match", fontsize=3.8, color=MID,
        ha="center", va="bottom")
status_mark(x + 1.79, 2.88, False, r=0.048)

pill(x + 0.23, 2.60, 1.57, 0.22,
     "target-best reads retained", "#DDF2EB", fontsize=5.35)

# Step 2: two amplicon-specific length windows on one shared size scale. This
# avoids implying that JH4 and N114 use the same lower bound.
x = top_x[1]
scale_left, scale_right = x + 0.48, x + 1.69
scale_w = scale_right - scale_left
ax.text(x + 1.08, 3.82, "accepted size windows", fontsize=5.1,
        fontweight="bold", color=MID, ha="center", va="center")
ax.text(scale_left, 3.68, "0", fontsize=4.2, color=MID, ha="center", va="bottom")
ax.text(scale_right, 3.68, "2,000 bp", fontsize=4.2, color=MID, ha="center", va="bottom")

for label, yy, lower, upper, threshold in [
    (r"$\mathbf{J}_{\mathbf{H}\mathbf{4}}$", 3.41, 800, 2000, "800–2,000 bp"),
    ("N114", 3.02, 400, 900, "400–900 bp"),
]:
    # Apply one visual hierarchy to both amplicon names and both ranges.
    ax.text(x + 0.25, yy + 0.025, label, fontsize=6.0,
            fontweight="bold", color=INK, ha="center", va="center")
    read_bar(scale_left, yy, scale_w, LIGHT, h=0.060)
    start = scale_left + (lower / 2000.0) * scale_w
    end = scale_left + (upper / 2000.0) * scale_w
    read_bar(start, yy, end - start, TEAL, h=0.060)
    ax.plot([start, start], [yy - 0.055, yy + 0.115], color=TEAL,
            linewidth=0.55, zorder=6)
    ax.plot([end, end], [yy - 0.055, yy + 0.115], color=TEAL,
            linewidth=0.55, zorder=6)
    status_mark(end + 0.10, yy + 0.030, True, r=0.043)
    ax.text(x + 1.08, yy - 0.13, threshold, fontsize=5.4,
            fontweight="bold", color=INK,
            ha="center", va="top")

pill(x + 0.41, 2.60, 1.28, 0.22,
     "amplicon-sized reads", "#DDF2EB", fontsize=5.25)

# Step 3: mean expected-error filter.
x = top_x[2]
for i, ok in enumerate([False, True, False, True]):
    yy = 3.67 - i * 0.24
    read_bar(x + 0.25, yy, 1.20, TEAL)
    for j in range(6 if not ok else 1):
        ax.add_patch(Circle((x + 0.39 + j * 0.17, yy + 0.027), 0.022,
                            facecolor=ORANGE, edgecolor="none", zorder=7))
    status_mark(x + 1.69, yy + 0.027, ok)
ax.plot([x + 0.22, x + 1.79], [3.31, 3.31], color=LIGHT, linewidth=0.55,
        linestyle=(0, (3, 2)))
pill(x + 0.39, 2.62, 1.25, 0.27, "EE ≤ 0.01 per base", PALE_YELLOW, fontsize=6.5)

# Step 4: sampled depth; explicitly avoid a deduplication claim.
x, y = bottom_x4, bottom_y
# Centred at 1.58, not 1.64: this label is now two lines, and at 1.64 the upper
# line runs into the amber card header. Bars top out at y = 1.34, so 1.58 sits
# clear of both.
ax.text(x + 0.64, 1.58, "delivered 2,336–5,000\nQC-passing 911–2,903", fontsize=5.2,
        color=MID, ha="center", va="center", linespacing=1.18)
ax.text(x + 2.00, 1.58, "matched depth", fontsize=5.6,
        color=MID, ha="center", va="center")

# Bars are schematic: they show that libraries entering this step differ in
# depth and leave it equal. The earlier version labelled them 1.50x, 1.25x and
# 1.00x, which invented per-sample ratios that are not in the data. The real
# per-sample counts are in Supplementary Table S1 and are now summarised in the
# header above, so the invented multipliers have been removed.
input_centres = [x + 0.34, x + 0.65, x + 0.96]
input_counts = [6, 5, 4]
for centre, count in zip(input_centres, input_counts):
    for i in range(count):
        read_bar(centre - 0.105, 1.01 + i * 0.060, 0.21, TEAL, h=0.032)

# The input block ends at x + 1.073 and the funnel starts at x + 1.270, so the
# first arrow has 0.197 in to sit in and is fine. The outbound side had only
# 0.087 in between the funnel and the first output group, which is less than the
# length of a mutation_scale=6 arrowhead, so the head was landing on the bars.
# The output groups move right by 0.06 in (which still leaves ~0.10 in of card
# margin) and the second arrow is centred in the widened gap.
arrow(x + 1.13, 1.15, x + 1.24, 1.15, color=NAVY, ms=6, lw=0.8)
downsampling_funnel(x + 1.36, 1.15, 0.86)
arrow(x + 1.479, 1.15, x + 1.569, 1.15, color=NAVY, ms=6, lw=0.8)

# All three samples are reduced to the same four-bar schematic height. Shifted
# right by 0.06 in from x + 1.65/1.94/2.23 to clear the inbound arrowhead.
for centre in [x + 1.71, x + 2.00, x + 2.29]:
    for i in range(4):
        read_bar(centre - 0.105, 1.01 + i * 0.060, 0.21, TEAL, h=0.032)

# The pill is lifted from y = 0.53 to 0.60 to make room for a three-line
# caption. Card 4 is 2.50 in wide, so the caption must wrap to about 50
# characters per line; set on two lines it ran 0.11 in past both card edges.
pill(x + 0.49, 0.60, 1.52, 0.27, "~500 reads per sample", "#DCEAF7", fontsize=6.5)
ax.text(x + 1.25, 0.30,
        "random down-sampling standardises sampling depth,"
        "\nnot per-read accuracy; PCR duplicates are not"
        "\nidentified without UMIs",
        fontsize=5.0, color=MID, ha="center", va="bottom", linespacing=1.18)

# Step 5: reference-aligned reads and final call definition.
x, y = bottom_x5, bottom_y
for r in range(4):
    yy = 1.59 - r * 0.20
    read_bar(x + 0.23, yy, 1.47, TEAL, h=0.052)
    for j, sx in enumerate([x + 0.64, x + 1.05, x + 1.48]):
        if (r + j) % 2 == 0:
            ax.add_patch(Rectangle((sx, yy), 0.045, 0.052, facecolor=PURPLE, edgecolor="none", zorder=7))
    if r == 2:
        ax.add_patch(Rectangle((x + 0.87, yy), 0.11, 0.052, facecolor=WHITE,
                               edgecolor=VERMILION, linewidth=0.65, zorder=7))
        status_mark(x + 1.79, yy + 0.026, False, 0.05)
ref_y = 0.80
ax.plot([x + 0.22, x + 1.70], [ref_y, ref_y], color=MID, linewidth=1.15)
ax.text(x + 0.22, ref_y - 0.08, "reference", fontsize=5.5, color=MID, ha="left", va="top")
for sx in [x + 0.64, x + 1.05, x + 1.48]:
    ax.plot([sx + 0.022, sx + 0.022], [ref_y + 0.02, 1.67], color=PURPLE,
            linewidth=0.4, linestyle=(0, (2, 2)), alpha=0.55)
# Right-hand criteria column. Previously centred at x + 2.19 with the longest
# line set on one line, which put its left edge on top of the read bars (which
# end at x + 1.70) and on the exclusion mark. The column is shifted right to
# x + 2.30, the longest label is wrapped, and the sizes are trimmed slightly so
# that every line clears the bars on the left and the card edge on the right.
ax.text(x + 2.30, 1.44, "SNPs only", fontsize=6.5, fontweight="bold", color=INK, ha="center")
# Calling is SHMtrace's quality-calibrated test against the negative controls
# (Benjamini-Hochberg FDR 0.05); previously the Geneious variant P <= 0.05.
ax.text(x + 2.30, 1.20, "vs. controls,\nFDR ≤ 0.05", fontsize=5.5, color=INK, ha="center",
        va="center", linespacing=1.12)
ax.text(x + 2.30, 0.98, "<50% variant\nfrequency", fontsize=5.5, color=MID,
        ha="center", va="center", linespacing=1.16)
ax.text(x + 2.30, 0.76, "indels excluded", fontsize=5.5, color=MID, ha="center")
pill(x + 0.52, 0.37, 1.82, 0.25, "filtered unique SNV set (SHMtrace)", PALE_PURPLE, fontsize=5.8)


basename = "Figure1_study_design_pipeline"
paths = [
    OUT / f"{basename}.png",
    OUT / f"{basename}.pdf",
    OUT / f"{basename}.svg",
    PDF_OUT / f"{basename}.pdf",
]
fig.savefig(paths[0], dpi=600, bbox_inches="tight", pad_inches=0.03)
fig.savefig(paths[1], bbox_inches="tight", pad_inches=0.03)
fig.savefig(paths[2], bbox_inches="tight", pad_inches=0.03)
fig.savefig(paths[3], bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
for path in paths:
    print(path)
