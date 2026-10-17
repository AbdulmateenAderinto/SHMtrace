"""Data and equations for the SciCanvas Figure S14 (figS14_engine.js).

    python docs/calibration/crm/scicanvas/engine_data.py

Reads what figs/fig_engine.py reads (out/intuition_jh4.json, out/engine_alignment.json and the J_H4
read accounting), reduces each panel to the numbers it plots, renders the six equations as small SVGs
(mathtext in Arial, drawn as outlines) and writes all of it into figS14_engine.js between the DATA
markers, so the builder stays one file that can also be pasted into SciCanvas.
"""
import base64
import io
import json
import re
import sys
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "figs"))
from fig_engine import ACCT, D, alignment_points  # noqa: E402

BUILDER = HERE / "figS14_engine.js"
PT = 1400 / (174 / 25.4 * 72)             # page units per point at 174 mm
EQ_PT = 7.0

EQUATIONS = [
    r"$\bar{E}=\frac{1}{L}\sum_{i=1}^{L}10^{-Q_i/10}\;\leq\;0.01$",
    r"$t^{*}(r)=\arg\max_{t}\;\mathrm{score}(r,\,t)$",
    # the rate the code uses (null_model.qrate): the pooled rate g for this change, strand and quality is the
    # floor, and a position's own control count is shrunk toward it with m pseudo-observations
    r"$\hat{q}_{x}=\kappa\,\max\!\left(g,\;\frac{k_{x}+m\,g}{n_{x}+m}\right),\quad m=500$",
    r"$\lambda_j=\sum_{s}n_s\,\hat{q}_s(\geq q_{(j)}),\quad"
    r"p=k\,\min_{j\leq k}P\!\left(X\geq j\right),\;X\sim\mathrm{Poisson}(\lambda_j)$",
    r"$\widehat{FP}(t)=\frac{1}{|D|}\sum_{d\in D}\left|\{p^{d}\leq t\}\right|\max\!\left(1,\frac{B}{B_d}\right),"
    r"\quad t^{*}=\max\{t:\widehat{FP}(t)\leq\alpha\left|\{p\leq t\}\right|\}$",
    r"$\rho=\frac{\#\,\mathrm{unique\;calls}}{\#\,\mathrm{analysed\;bases}}\times10^{5},"
    r"\quad \mathrm{Welch\;}t\mathrm{\;on\;}\log_{10}\rho,\;n=\mathrm{animals}$",
]


def equation(tex):
    """One equation as an SVG data URI, with its size in page units."""
    # glyphs as outlines, so the SVG needs no font; Arial, with STIX Sans for the symbols Arial lacks
    plt.rcParams.update({"svg.fonttype": "path", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "mathtext.sf": "Arial",
                         "mathtext.fallback": "stixsans", "font.family": "Arial"})
    fig = plt.figure(figsize=(6, 1))
    fig.text(0, 0.5, tex, fontsize=EQ_PT, color="#1d2733", va="center")
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.01, transparent=True)
    plt.close(fig)
    svg = buf.getvalue()
    w, h = (float(re.search(rf'{k}="([\d.]+)pt"', svg).group(1)) for k in ("width", "height"))
    uri = "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()
    return {"src": uri, "w": round(w * PT, 1), "h": round(h * PT, 1)}


def counts(ps, grid):
    ps = np.sort(np.asarray(ps, float))
    return [int(np.searchsorted(ps, 10 ** g, "right")) for g in grid]


def main():
    parts = [("too_long", "low_quality"), ("not_target", "chimeric", "short_alignment", "paralogue_alleles"),
             ("usable",)]
    p1 = {"libs": list(ACCT), "pct": [[round(100 * sum(a[k] for k in ks) / a["delivered"], 2) for ks in parts]
                                      for a in ACCT.values()]}
    al = alignment_points()
    p2 = {"target": al["target"][::2], "paralogue": al["paralogue"]}
    xs, ys = [], []
    for k, (mm, bases) in sorted(D["quality"].items(), key=lambda kv: int(kv[0])):
        if bases >= 2000:
            xs.append(int(k))
            ys.append(round(-10 * np.log10(max(mm, 0.5) / bases), 2))
    p3 = {"x": xs, "y": ys}
    sites, case = D["anatomy"]["sites"], D["anatomy"]["case"]
    p4 = []
    for key in ("true", "error"):
        s = sites[key]
        libs = s["libraries"]
        p4.append({"label": [f"{s['ref']}{s['pos']}{s['alt']}", s["context"]],
                   "case": libs[case]["alt_plus"] + libs[case]["alt_minus"],
                   "ctl": round(float(np.mean([v["alt_plus"] + v["alt_minus"] for k, v in libs.items() if k != case])), 2),
                   "called": bool(s["called"])})
    grid = list(np.round(np.linspace(-12, 0, 61), 3))
    dec = [d[0] if isinstance(d, list) and d and isinstance(d[0], list) else d for _, d in sorted(D["decoys"].items())]
    p5 = {"grid": grid, "case": counts(D["candidates"][case]["p"], grid), "decoys": [counts(d, grid) for d in dec],
          "t": D["anatomy"]["threshold"]}
    bases = (D["window"][1] - D["window"][0] + 1) * D["depth"]
    p6 = [[round(float(np.log10(max(len(D["calls"][n]), 0.5) / bases * 1e5)), 4) for n in names if n in D["calls"]]
          for names in (D["cases"], D["controls"])]
    data = {"p1": p1, "p2": p2, "p3": p3, "p4": p4, "p5": p5, "p6": p6, "eq": [equation(t) for t in EQUATIONS]}
    src = BUILDER.read_text()
    new = re.sub(r"/\*DATA\*/.*?/\*END\*/", lambda _m: "/*DATA*/" + json.dumps(data, separators=(",", ":")) + "/*END*/",
                 src, flags=re.S)
    BUILDER.write_text(new)
    print(f"{BUILDER.name}: data {len(new) - len(src):+d} bytes; equations "
          + ", ".join(f"{e['w']:.0f}x{e['h']:.0f}" for e in data["eq"]))


if __name__ == "__main__":
    main()
