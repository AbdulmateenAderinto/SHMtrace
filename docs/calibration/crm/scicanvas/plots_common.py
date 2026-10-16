# ruff: noqa: F821  (plt, PAGE and FONT_SIZE come from the SciCanvas plot worker and build.js's preamble)
# Shared by the CRM code plots (SciCanvas 2.0). Prepended to each plot's code by the builders, so every plot in a
# .scifig re-runs on its own in SciCanvas. PAGE and the rcParams scaling come from build.js's preamble: sizes below
# are the printed sizes, in points.
import json  # noqa: F401  (each plot reads its numbers with json.loads)

INK, MID, SOFT, FAINT = "#1d2733", "#4a5562", "#8a96a3", "#cfd8df"
SH, PUB = "#1f6f8b", "#a3621b"           # SHMtrace, published
CASE, CASEF = "#274b7e", "#6d9ad0"        # wild type
KO, KOF = "#7a4f9a", "#bda3d4"            # knockout
CTRL, CTRLF = "#9a3439", "#d07a7e"        # AID-null
plt.rcParams.update({"axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK, "axes.labelcolor": MID,
                     "text.color": INK})


def area(d):
    """Scatter size for a marker d printed points across."""
    return (d * PAGE) ** 2


def lw(w):
    """A line w printed points wide."""
    return w * PAGE


def plain_log(axis):
    """Log-axis tick labels as plain numbers (10, 100), not 10^n, whose exponents print below 6 pt."""
    axis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_formatter(plt.NullFormatter())


def strip(ax, x, vals, face, edge, d=5.0, spread=0.10, marker="o", zorder=3):
    """One group's points, spread sideways so they do not hide each other."""
    vals = list(vals)
    if not vals:
        return
    xs = [x] if len(vals) == 1 else [x - spread + 2 * spread * k / (len(vals) - 1) for k in range(len(vals))]
    ax.scatter(xs, vals, s=area(d), facecolor=face, edgecolor=edge, linewidth=lw(0.6), marker=marker, zorder=zorder)


def mean_bar(ax, x, vals, color, half=0.16):
    vals = list(vals)
    if vals:
        m = sum(vals) / len(vals)
        ax.plot([x - half, x + half], [m, m], color=color, linewidth=lw(1.3), solid_capstyle="butt", zorder=4)
