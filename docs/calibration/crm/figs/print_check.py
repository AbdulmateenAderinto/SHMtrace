"""Cell Press print minimums for the matplotlib figures: text >= 6 pt (super- and subscripts too), lines 0.5-1.5 pt.

    python docs/calibration/crm/figs/print_check.py [Figure2 FigureS7 ...]

Reads each figure's SVG, which matplotlib writes at its printed size (1 SVG unit = 1 pt), and lists every text run
and stroke below the minimum with the text it belongs to. The SciCanvas figures are checked by
scicanvas/build.js on every build; this does the same for the figures drawn by figs/*.py. Exits 1 on any breach.
"""
import re
import sys
from pathlib import Path

from style import OUT

MIN_PT, MIN_LW, MAX_LW = 5.95, 0.495, 1.505
STAR = str(Path(__file__).resolve().parents[1] / "star" / "figures")
FIGURES = ["Figure2", "Figure3", "Figure4", "Figure5", "Figure6", "Figure7", "FigureS2", "FigureS7", "FigureS9",
           "FigureS10", "FigureS23", "FigureS19_switch_pacbio", "FigureS18_public_smu",
           "FigureS17_engineered", "FigureS4_fdr_calibration", "FigureS5_controls", "FigureS11_immcantation", f"{STAR}/STAR_Figure3", f"{STAR}/STAR_Figure4"]


def check(name):
    svg = (Path(name + ".svg") if "/" in name else OUT / f"{name}.svg").read_text()   # a path, or a name in figures/
    bad = []
    # every element or tspan that sets a font size; the text is what follows it up to the next tag
    for m in re.finditer(r'font(?:-size)?:\s*(?:[a-z ]*?)([\d.]+)px[^>]*>([^<]*)', svg):
        size, text = float(m.group(1)), m.group(2).strip()
        if size < MIN_PT and text:
            bad.append(f"text {size:.2f} pt: {text[:40]!r}")
    for m in re.finditer(r'stroke-width:\s*([\d.]+)', svg):
        w = float(m.group(1))
        if 0.01 < w < MIN_LW or w > MAX_LW:
            bad.append(f"stroke {w:.2f} pt")
    return bad


def main(names):
    failed = 0
    for name in names or FIGURES:
        bad = check(name)
        strokes = sum(1 for b in bad if b.startswith("stroke"))
        texts = sorted({b for b in bad if b.startswith("text")})
        print(f"{name}: {'ok' if not bad else f'{len(texts)} text, {strokes} stroke breaches'}")
        for b in texts[:12]:
            print("   ", b)
        if strokes:
            print("    strokes:", sorted({b for b in bad if b.startswith("stroke")}))
        failed += bool(bad)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
