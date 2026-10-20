"""Visualisations for an SHM run — the plots people actually put in papers.

Three generalisable views, computed from data the engine already produces:
  1. substitution spectrum  (the 12 from->to classes
  G:C vs A:T coloured)
  2. spatial mutation map    (unique variants per position along the reference)
  3. mutations-per-read load  (distribution of SHM load
  saturation readout)

Everything is assay-agnostic: it works for any reference/window, so the same
code serves JH4, V regions, or a user's own locus. Palette matches the
manuscript figures. Each function returns PNG bytes (for web embedding) and can
optionally write a file.
"""
from __future__ import annotations

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .variants import BASES, Variant

# manuscript palette (fill, edge)
BLUE = ("#6D9AD0", "#274B7E")     # G:C-origin
RED = ("#D07A7E", "#9A3439")      # A:T-origin
GREY = "#C9C9C9"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.edgecolor": "black", "savefig.dpi": 150})


def _bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def _finish(fig, path: str | None) -> bytes:
    if path:
        fig.savefig(path, bbox_inches="tight", facecolor="white")
    return _bytes(fig)


def spectrum(spec: dict, title: str = "Substitution spectrum",
             path: str | None = None) -> bytes:
    """12-class from->to bar chart
    G:C origins blue, A:T origins red."""
    classes = [(r, a) for r in BASES for a in BASES if r != a]
    vals = [spec.get(r, {}).get(a, 0) for r, a in classes]
    colors = [BLUE if r in ("G", "C") else RED for r, _ in classes]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    x = range(len(classes))
    ax.bar(x, vals, color=[c[0] for c in colors],
           edgecolor=[c[1] for c in colors], linewidth=1.3, zorder=3)
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{r}→{a}" for r, a in classes], fontsize=9)
    ax.set_ylabel("Unique variants")
    ax.set_title(title, fontweight="bold", fontsize=12)
    ax.yaxis.grid(True, linestyle=(0, (4, 3)), color=GREY, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # legend proxies
    ax.bar(0, 0, color=BLUE[0], edgecolor=BLUE[1], label="G:C origin")
    ax.bar(0, 0, color=RED[0], edgecolor=RED[1], label="A:T origin")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    return _finish(fig, path)


def spatial(variants: list[Variant], reflen: int,
            window: tuple | None = None,
            title: str = "Mutation position map",
            path: str | None = None) -> bytes:
    """Stem plot: unique variants per reference position, coloured by ref base."""
    per_pos = {}
    for v in variants:
        per_pos.setdefault(v.pos, [0, v.ref])
        per_pos[v.pos][0] += 1
    fig, ax = plt.subplots(figsize=(8.4, 3.0))
    for pos, (n, refbase) in per_pos.items():
        c = BLUE[1] if refbase in ("G", "C") else RED[1]
        ax.vlines(pos, 0, n, color=c, linewidth=1.4, zorder=3)
        ax.plot(pos, n, "o", color=c, markersize=3, zorder=4)
    if window:
        ax.axvspan(window[0], window[1], color="#EEF2F8", zorder=0)
    ax.set_xlim(1, reflen)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Reference position (bp)")
    ax.set_ylabel("Unique variants")
    ax.set_title(title, fontweight="bold", fontsize=12)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return _finish(fig, path)


def mutation_load(per_read: list[int],
                  title: str = "Mutations per read",
                  path: str | None = None) -> bytes:
    """Histogram of SHM load per read (called-variant sites carried per read)."""
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    if per_read:
        hi = max(per_read)
        bins = range(0, hi + 2)
        ax.hist(per_read, bins=bins, color=BLUE[0], edgecolor=BLUE[1],
                linewidth=1.2, align="left", zorder=3)
    ax.set_xlabel("Called variants per read")
    ax.set_ylabel("Reads")
    ax.set_title(title, fontweight="bold", fontsize=12)
    ax.yaxis.grid(True, linestyle=(0, (4, 3)), color=GREY, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return _finish(fig, path)


def compare_spectrum(res: dict, path: str | None = None) -> bytes:
    """Grouped 12-class spectrum, CONTROL vs CASE (base-composition uncorrected counts)."""
    spec = res["spectrum"]
    classes = [(r, a) for r in BASES for a in BASES if r != a]
    cvals = [spec[r][a]["control"] for r, a in classes]
    avals = [spec[r][a]["case"] for r, a in classes]
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    x = range(len(classes))
    w = 0.4
    ax.bar([i - w/2 for i in x], cvals, w, color=BLUE[0], edgecolor=BLUE[1],
           linewidth=1.1, label="Control", zorder=3)
    ax.bar([i + w/2 for i in x], avals, w, color=RED[0], edgecolor=RED[1],
           linewidth=1.1, label="Case", zorder=3)
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{r}→{a}" for r, a in classes], fontsize=9)
    ax.set_ylabel("Unique variants")
    ax.set_title("Spectrum: control vs case",
                                                    fontweight="bold", fontsize=12)
    ax.yaxis.grid(True, linestyle=(0, (4, 3)), color=GREY, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    return _finish(fig, path)


def frequency_bars(res: dict, path: str | None = None) -> bytes:
    """Base-composition-corrected mutation frequencies, control vs case, with stars."""
    rows = res["frequency_rows"]
    labels = [r["metric"].replace(" mutations", "").replace(" mutation frequency", "")
              for r in rows]
    cvals = [r["control_freq"] for r in rows]
    avals = [r["case_freq"] for r in rows]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    y = range(len(rows))
    h = 0.38
    ax.barh([i + h/2 for i in y], cvals, h, color=BLUE[0], edgecolor=BLUE[1], label="Control", zorder=3)
    ax.barh([i - h/2 for i in y], avals, h, color=RED[0], edgecolor=RED[1], label="Case", zorder=3)
    for i, r in enumerate(rows):
        ax.text(max(r["control_freq"], r["case_freq"]) * 1.02, i, r["sig"],
                va="center", fontsize=9, color="#333")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Mutation frequency (base-composition corrected)")
    ax.set_title("Mutation frequencies: control vs case (pooled, legacy)", fontweight="bold", fontsize=12)
    ax.invert_yaxis()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    return _finish(fig, path)


def motif_enrichment(res: dict, path: str | None = None) -> bytes:
    """Hot/cold-spot motif mutation frequencies, control vs case."""
    rows = res["motif_rows"]
    labels = [r["metric"].replace("_", " ") for r in rows]
    cvals = [r["control_freq"] for r in rows]
    avals = [r["case_freq"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    x = range(len(rows))
    w = 0.38
    ax.bar([i - w/2 for i in x], cvals, w, color=BLUE[0], edgecolor=BLUE[1], label="Control", zorder=3)
    ax.bar([i + w/2 for i in x], avals, w, color=RED[0], edgecolor=RED[1], label="Case", zorder=3)
    for i, r in enumerate(rows):
        ax.text(i, max(r["control_freq"], r["case_freq"]) * 1.02, r["sig"],
                ha="center", fontsize=9, color="#333")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Motif mutation frequency")
    ax.set_title("Hot/cold-spot targeting (pooled, legacy)", fontweight="bold", fontsize=12)
    ax.yaxis.grid(True, linestyle=(0, (4, 3)), color=GREY, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    return _finish(fig, path)


def comparison_plots(res: dict, out_dir: str | None = None) -> dict:
    import os
    def p(n):
        return os.path.join(out_dir, n) if out_dir else None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    return {"spectrum_compare": compare_spectrum(res, p("spectrum_compare.png")),
            "frequency_bars": frequency_bars(res, p("frequency_bars.png")),
            "motif_enrichment": motif_enrichment(res, p("motif_enrichment.png"))}


def render_all(result: dict, out_dir: str | None = None) -> dict:
    """Return {name: png_bytes} for all plots
    also write files if out_dir set."""
    import os
    s = result["summary"]
    win = s["params"].get("analysis_window")
    reflen = s["calling"]["reference_length"]
    variants = result["variants"]
    per_read = result.get("per_read_mutations", [])

    def p(name):
        return os.path.join(out_dir, name) if out_dir else None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    return {
        "spectrum": spectrum(s["spectrum_unique"], path=p("spectrum.png")),
        "spatial": spatial(variants, reflen, tuple(win) if win else None,
                           path=p("spatial.png")),
        "mutation_load": mutation_load(per_read, path=p("mutation_load.png")),
    }




# ---- project figures (many samples, groups) ---------------------------------
# These live in viz_project so this module stays about one sample; they are
# re-exported here because callers (report.py, the web app) import viz.
def project_figures(res: dict, fmt: str = "png") -> dict[str, bytes]:
    from .viz_project import project_figures as _pf
    return _pf(res, fmt=fmt)


def __getattr__(name: str):
    if name in ("group_dotplot", "group_spectrum", "group_positional", "group_load",
                "group_matrix", "GROUP_COLOURS", "SUBS12"):
        from . import viz_project
        return getattr(viz_project, name)
    raise AttributeError(name)
