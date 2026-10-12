"""What the manuscript's Sanger analysis counted for each mouse, against what is counted now.

    python docs/calibration/sanger_inclusion_figure.py

Reads figures/sanger_vs_nanopore/included.json: the composition of the manuscript's
Geneious assemblies (every clone classed by the mouse it really comes from, plus paralogue
sequences counted only in the denominator) and of the organised AB1 set after nanoshm's
filters. Writes figS_sanger_inclusion.{png,pdf,svg} next to it.
"""
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent / "figures/sanger_vs_nanopore"
D = json.load(open(HERE / "included.json"))
SANGER = Path.home() / "Nanoshm/sanger_ab1_by_mouse"


def refresh_now() -> None:
    """The 'now' side, from the organised set's manifest and the comparison's stats.

    Clones counted = the comparison's analysed Sanger sequences. Traces were filtered
    when the set was organised, so the GenBank-only clones counted are the rest.
    """
    import csv
    man = list(csv.DictReader(open(SANGER / "manifest.csv")))
    st = json.load(open(HERE / "recommended/stats.json"))
    folder = {"WT2": "WT2", "WT3": "WT3", "WT4": "WT4", "A1": "Aicda_A1", "A2": "Aicda_A2",
              "A3": "Aicda_A3"}
    for m, f in folder.items():
        kept = sum(1 for r in man if r["mouse"] == f and r["status"] == "kept")
        rel = sum(1 for r in man if r["mouse"] == f and r["status"] == "relabelled")
        n = st["sequences"][f"sanger:{m}"]
        u = st["unique_mutations"][f"sanger:{m}"]
        D["now"][m] = {"composition": {"own traces": kept, "relabelled traces": rel,
                                       "GenBank-only": n - kept - rel},
                       "clones_counted": n, "unique_mutations": u,
                       "density_x1e5": round(1e5 * max(u, 1) / (n * 529), 1)}
    json.dump(D, open(HERE / "included.json", "w"), indent=1)


refresh_now()
MICE = ["WT2", "WT3", "WT4", "A1", "A2", "A3"]
LABEL = {"WT2": "Wt2", "WT3": "Wt3", "WT4": "Wt4", "A1": "Aicda⁻/⁻ 1", "A2": "Aicda⁻/⁻ 2",
         "A3": "Aicda⁻/⁻ 3"}

# segment colours: correct material in blues, everything wrongly counted in warm tones
SEG_MS = [("own", "this mouse's clones", "#2F6BA6", None),
          ("WT2", "another wild-type mouse's clones", "#E69F00", None),
          ("WT3", None, "#E69F00", None),
          ("A1", "another Aicda⁻/⁻ mouse's clones", "#F0C27A", None),
          ("A2", None, "#F0C27A", None),
          ("MKLN", "MKLN-knockout clones", "#CC79A7", None),
          ("JH6", "JH6-amplicon clones", "#9A3439", None),
          ("w1", "unattributable (w1)", "#9A9A9A", None),
          ("paralogue (denominator only)", "paralogue (counted as unmutated JH4)", "#FFFFFF", "////")]
SEG_NOW = [("own traces", "AB1 trace, own label", "#2F6BA6", None),
           ("relabelled traces", "AB1 trace, relabelled by sequence", "#2F6BA6", "xxx"),
           ("GenBank-only", "GenBank record (trace not found)", "#A8CBE8", None)]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300})
fig = plt.figure(figsize=(15.5, 9.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.75, 1], height_ratios=[1, .62], hspace=.42, wspace=.22)

# ---------------- A: composition
ax = fig.add_subplot(gs[:, 0])
y = 0
yt, yl = [], []
for m in MICE:
    for which, segs, key in (("Manuscript", SEG_MS, "manuscript"), ("Now", SEG_NOW, "now")):
        comp = D[key][m]["composition"]
        x = 0
        for k, _lab, col, hatch in segs:
            v = comp.get(k, 0)
            if not v:
                continue
            ax.barh(y, v, left=x, height=.72, color=col, hatch=hatch, edgecolor="#333", lw=.7)
            if v >= 3:
                light = col in ("#FFFFFF", "#A8CBE8", "#F0C27A", "#E69F00")
                ax.text(x + v / 2, y, str(v), ha="center", va="center", fontsize=8.5,
                        color="#1a1a1a" if light else "white", fontweight="bold")
            x += v
        d = D[key][m]
        dens = d["density_x1e5"]
        ax.text(x + 1, y, f"n = {d['clones_counted']} · {d['unique_mutations']} mutations · "
                f"{dens:.0f} × 10⁻⁵", va="center", fontsize=8.5, color="#333")
        yt.append(y)
        yl.append(which)
        y += .85
    ax.text(-0.115, y - 1.27, LABEL[m], ha="right", va="center", fontsize=11, fontweight="bold",
            fontstyle="italic" if m.startswith("A") else "normal",
            transform=matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData))
    y += .75
ax.set_yticks(yt, yl, fontsize=9)
ax.tick_params(axis="y", length=0, pad=4)
ax.invert_yaxis()
ax.set_xlim(0, 75)
ax.set_xlabel("Sequences counted in the denominator")
ax.spines["left"].set_visible(False)
ax.set_title("A   What each mouse's Sanger value was built from", loc="left", fontweight="bold",
             fontsize=12, pad=12)
handles = [Patch(fc=c, hatch=h, ec="#333", label=lab) for _k, lab, c, h in SEG_MS if lab]
handles += [Patch(fc=c, hatch=h, ec="#333", label=lab) for _k, lab, c, h in SEG_NOW]
leg1 = ax.legend(handles=handles[:-3], loc="lower right", frameon=False, fontsize=8.5,
                 title="Manuscript (Geneious assemblies)", title_fontsize=9,
                 bbox_to_anchor=(1.0, 0.25), ncol=1, alignment="left")
ax.add_artist(leg1)
# the 'now' key sits in the empty right of the panel, beside the knockout rows
ax.legend(handles=handles[-3:], loc="lower right", frameon=False, fontsize=8.5,
          title="Now (organised AB1 set)", title_fontsize=9, bbox_to_anchor=(1.0, 0.03),
          alignment="left")

# ---------------- B: density per mouse, before and after
bx = fig.add_subplot(gs[0, 1])
ends = []
for m in MICE:
    a, b = D["manuscript"][m], D["now"][m]
    wt = m.startswith("WT")
    col = "#274B7E" if wt else "#9A3439"
    va, vb = a["density_x1e5"], b["density_x1e5"]
    bx.plot([0, 1], [va, vb], color=col, lw=1.6, alpha=.8)
    for xx, v, d in ((0, va, a), (1, vb, b)):
        cens = d["unique_mutations"] == 0
        bx.scatter(xx, v, s=55, zorder=3, facecolor="white" if not cens else "#eeeeee",
                   edgecolor=col, lw=1.6, linestyle="--" if cens else "-")
    ends.append((vb, LABEL[m], col, wt))
bx.set_yscale("log")
# labels at the line ends, pushed apart in log space so they never overlap
placed = []
for v, lab, col, wt in sorted(ends):
    ly = math.log10(v)
    if placed and ly - placed[-1] < .12:
        ly = placed[-1] + .12
    placed.append(ly)
    bx.annotate(lab, (1, v), xytext=(1.1, 10 ** ly), textcoords="data", va="center", fontsize=8.5,
                color=col, fontstyle="normal" if wt else "italic",
                arrowprops=dict(arrowstyle="-", color=col, lw=.6, alpha=.6))
bx.set_xlim(-.25, 1.75)
bx.set_xticks([0, 1], ["Manuscript", "Now"])
bx.set_ylabel("Unique mutations / bp (× 10⁻⁵)")
bx.yaxis.grid(True, ls=(0, (5, 4)), color="#C9C9C9")
bx.set_axisbelow(True)
bx.set_title("B   Each mouse's Sanger density", loc="left", fontweight="bold", fontsize=12)
bx.text(0, -.17, "grey dashed: no mutation, shown at one mutation's worth",
        transform=bx.transAxes, fontsize=8, color="#555")

# ---------------- C: the group comparison it feeds
cx = fig.add_subplot(gs[1, 1])
cx.axis("off")


def group(key):
    wt = [D[key][m]["density_x1e5"] for m in MICE[:3]]
    ko = [D[key][m]["density_x1e5"] for m in MICE[3:]]
    gw = 10 ** (sum(map(math.log10, wt)) / 3)
    gk = 10 ** (sum(map(math.log10, ko)) / 3)
    p = stats.ttest_ind([math.log10(v) for v in wt], [math.log10(v) for v in ko], equal_var=False).pvalue
    return gw, gk, gw / gk, p


rows = [("", "Wt", "Aicda⁻/⁻", "Fold", "P")]
for key, lab in (("manuscript", "Sanger, manuscript"), ("now", "Sanger, now")):
    gw, gk, f, p = group(key)
    rows.append((lab, f"{gw:.0f}", f"{gk:.0f}", f"{f:.0f}×", f"{p:.4f}" if p >= .001 else f"{p:.1e}"))
rows.append(("Nanopore (unchanged)", "57", "0.55", "103×", "0.0017"))
for i, r in enumerate(rows):
    for j, cell in enumerate(r):
        cx.text([0, .56, .77, .9, 1.06][j], .92 - i * .2, cell, transform=cx.transAxes, fontsize=10,
                fontweight="bold" if i == 0 or (i == 2 and j == 4) else "normal",
                color="#9A3439" if (i == 2 and j == 4) else "#1a1a1a",
                ha="left" if j == 0 else "right")
cx.text(0, .05, "Geometric means, × 10⁻⁵; Welch's t on log10, three mice per group.\n"
        "Correctly assigned, Wt2's own clones sit at the knockout level,\n"
        "so Sanger alone no longer separates the groups; nanopore does.",
        transform=cx.transAxes, fontsize=8.5, color="#444", va="bottom")
cx.set_title("C   The group comparison", loc="left", fontweight="bold", fontsize=12)

for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"figS_sanger_inclusion.{ext}", bbox_inches="tight")
print(HERE / "figS_sanger_inclusion.png")
