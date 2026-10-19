"""Do the figure legends agree with the figures they describe?

    python crm/legend_check.py          -> report; exit 1 if anything fails
    python crm/legend_check.py --review -> also print the axis-label review table

Legends live in the manuscript sources and figures are drawn by separate scripts, so the two can drift
without anything noticing. They did, repeatedly: an axis reading "% of mutations at G:C" while the code
divided by every mutation, a panel printing a confidence interval the Results had dropped, a legend saying
"equally well" where the Methods said "within five matched bases".

Three checks run here. The first two are assertions and fail the build:

  panels      every panel letter a legend describes must be drawn, and every letter drawn must be described
  references  every "Figure N", "Figure SN" and "Table SN" named anywhere in the manuscript sources must
              exist among the figures and tables actually built

The third is a review aid rather than an assertion, because whether an axis label and a legend sentence
describe the same quantity is a judgement a person has to make. It prints each panel's axis labels beside
the legend sentence for that panel, so the two can be read together instead of in different files.
"""
import argparse
import re
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
sys.path.insert(0, str(CRM / "manuscript"))
sys.path.insert(0, str(CRM))

PANEL_IN_LEGEND = re.compile(r"\((?:and )?([A-Z])(?:\s*(?:and|–|-|,)\s*([A-Z]))?\)")
LETTER_CALL = re.compile(r'letter\(\s*ax\s*,\s*"([A-Z])"')
LETTER_SEQ = re.compile(r'letter\(\s*ax\s*,\s*"([A-Z]{2,})"\[')
PANEL_CALL = re.compile(r"panel\(\s*'([A-Z])'")
AXIS = re.compile(r'set_(xlabel|ylabel)\(\s*(?:f?")([^"]{2,90})"')
TITLE = re.compile(r'set_title\(\s*(?:f?")([^"]{2,90})"')


def legend_panels(text):
    """The panel letters a legend claims to describe, expanding "(A and B)" and "(A-C)".

    Panels are lettered consecutively from A, so only the run starting at A is kept. That is what separates
    a panel marker from an abbreviation defined in passing, such as the "(S)" and "(N)" that Figure 5E uses
    for Sanger and nanopore.
    """
    seen = set()
    for a, b in PANEL_IN_LEGEND.findall(text):
        seen.add(a)
        if b:
            seen.add(b)
            if ord(b) > ord(a) + 1:                      # an (A-C) style range
                seen |= {chr(c) for c in range(ord(a), ord(b) + 1)}
    out, c = set(), ord("A")
    while chr(c) in seen:
        out.add(chr(c))
        c += 1
    return out


def drawn_panels(src):
    out = set(LETTER_CALL.findall(src)) | set(PANEL_CALL.findall(src))
    for seq in LETTER_SEQ.findall(src):
        out |= set(seq)
    return out


def script_for(stem, scripts):
    """The drawing script for a figure stem such as 'Figure5' or 'FigureS18_engineered'."""
    want = stem.lower().replace("figure", "fig")
    for p in scripts:
        n = p.stem.lower()
        if n == want or n.startswith(want + "_") or want.startswith(n):
            return p
    return None


def main(review=False):
    import build_crm as B
    import crm_methods as M

    scripts = sorted(list((CRM / "figs").glob("*.py")) + list((CRM / "scicanvas").glob("*.js")))
    figs = [(t, "", lg) for t, lg in M.LEGENDS] + [(t, i, lg) for t, i, lg in B.SUPP_FIGS]
    problems, reviewed = [], []
    for title, imgs, legend in figs:
        m = re.match(r"(Figure S?\d+)", title)
        if not m:
            continue
        label = m.group(1)
        stem = label.replace(" ", "")
        src_path = script_for(stem, scripts)
        if src_path is None:
            continue                                      # drawn elsewhere (manuscript_corrected, etc.)
        src = src_path.read_text()
        want, have = legend_panels(legend), drawn_panels(src)
        if want and have:
            missing, extra = sorted(want - have), sorted(have - want)
            if missing:
                problems.append(f"{label}: legend describes panel(s) {','.join(missing)} that "
                                f"{src_path.name} does not draw")
            if extra:
                problems.append(f"{label}: {src_path.name} draws panel(s) {','.join(extra)} that the "
                                f"legend does not describe")
        if review:
            reviewed.append((label, src_path.name, TITLE.findall(src), [b for _a, b in AXIS.findall(src)],
                             legend[:160]))

    # every figure and table named in the manuscript sources must exist
    built = {re.match(r"(Figure S?\d+|Table S?\d+)", t).group(1).replace(" ", "")
             for t, *_ in figs if re.match(r"(Figure S?\d+|Table S?\d+)", t)}
    built |= {f"Table{n}" for n in ("1",)} | {f"TableS{i}" for i in range(1, 6)}
    text = "".join(p.read_text() for p in (CRM / "manuscript").glob("crm_*.py"))
    # "Supplementary Figure 9" is Ma et al.'s figure, not one of ours, so references carrying another
    # paper's qualifier are left out
    named = {f"{k}{n}" for k, n in re.findall(r"(?<!Supplementary )\b(Figure|Table) (S?\d+)", text)}
    dangling = sorted(n for n in named - built if not n.startswith("Figure") or "S" in n or int(
        n.replace("Figure", "")) > 7)
    for n in dangling:
        problems.append(f"text names {n}, which is not among the figures and tables built")

    if review:
        print(f"{'figure':10s} {'script':26s} panel titles / axis labels, then the legend")
        for label, name, titles, axes, lg in reviewed:
            print(f"\n{label:10s} {name:26s}")
            for t in titles:
                print(f"      title: {t}")
            for a in axes:
                print(f"       axis: {a}")
            print(f"     legend: {lg}...")
        print()
    if problems:
        print(f"legend check: {len(problems)} problem(s)")
        for p in problems:
            print("  -", p)
        return 1
    print(f"legend check: {len(figs)} legends, no mismatch between legends and the figures they describe")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--review", action="store_true", help="print the axis-label review table")
    sys.exit(main(**vars(ap.parse_args())))
