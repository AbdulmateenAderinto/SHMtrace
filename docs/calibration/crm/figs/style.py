"""Shared style for the Cell Reports Methods figures (Cell Press: 174 mm full width, Arial 6-8 pt,
lines 0.5-1.5 pt, RGB, no red-green pairing). Each caller keeps one colour in every figure."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "figures"
DATA = Path(__file__).resolve().parents[1] / "out"
MM = 1 / 25.4
FULL, HALF = 174 * MM, 85 * MM

CALLERS = {   # key: (label, colour, marker)
    "recommended": ("SHMtrace", "#1f6f8b", "o"),
    "geneious": ("Geneious rule", "#a3621b", "s"),
    "lofreq": ("LoFreq, defaults", "#9aa5b1", "v"),
    "lofreq_nobaq": ("LoFreq, BAQ off", "#6a4c93", "^"),
    "lofreq_permissive": ("LoFreq, permissive", "#c08fc7", "D"),
    "v0.1": ("Per-position null", "#5f6b75", "x"),
    "pbinom": ("Quality model only", "#2c5d8a", "P"),
    "clairs": ("ClairS", "#c99a06", "X"),
    "deepsomatic": ("DeepSomatic", "#b5487a", "P"),
}
CASE, CASEF, CTRL, CTRLF = "#274b7e", "#6d9ad0", "#9a3439", "#d07a7e"
INK, MID, SOFT = "#1d2733", "#4a5562", "#8a96a3"


def setup():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "lines.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "text.color": INK, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "legend.frameon": False, "axes.titleweight": "bold", "axes.titlelocation": "left",
    })


def letter(ax, s, dx=-0.18, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600}), ("svg", {})):
        fig.savefig(OUT / f"{name}.{ext}", facecolor="white", **kw)
    from PIL import Image
    Image.open(OUT / f"{name}.png").convert("RGB").save(OUT / f"{name}.tif", dpi=(600, 600), compression="tiff_lzw")
    w, h = fig.get_size_inches()
    print(f"{name}: {w / MM:.0f} x {h / MM:.0f} mm")
