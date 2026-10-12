"""Write sample sheets for the lab's JH4 Sanger traces from organize_jh4_ab1.py's manifest.

    python make_jh4_samplesheets.py ~/Nanoshm/sanger_ab1_original/manifest.csv

Every trace keeps the mouse it was filed under (GenBank folder; traces with no GenBank
record by their label), as in the manuscript. Traces the pipeline will drop (noisy,
paralogue) are listed anyway so the pipeline records why. Two clones never filed with a
mouse, which carry AID-type mutations under knockout labels, are listed as excluded. A
clone whose tube label names another mouse (the WT2/WT3 swap) stays where it was filed, with a
note saying so: the manuscript keeps the original publication's classification.

  jh4_manuscript.csv   Wt2-4 and Aicda-/- 1-3, as analysed in the manuscript
  jh4_with_mkln.csv    the same plus the four MKLN-knockout mice, for the QC checks
"""
import csv
import re
import sys
from pathlib import Path

HOME = str(Path.home())
STUDY = {"WT2": "WT", "WT3": "WT", "WT4": "WT", "Aicda_A1": "KO", "Aicda_A2": "KO", "Aicda_A3": "KO"}
MKLN = {"MKLN_K1": "MKLN", "MKLN_K2": "MKLN", "MKLN_K3": "MKLN", "MKLN_K4": "MKLN"}
FOLDER = {"W2": "WT2", "W3": "WT3", "W4": "WT4", "A1": "Aicda_A1", "A2": "Aicda_A2",
          "A3": "Aicda_A3", "K1": "MKLN_K1", "K2": "MKLN_K2", "K3": "MKLN_K3", "K4": "MKLN_K4"}
ALIASES = {"WT2": "w2;wt2", "WT3": "w3;wt3", "WT4": "w4;wt4", "Aicda_A1": "a1", "Aicda_A2": "a2",
           "Aicda_A3": "a3", "MKLN_K1": "k1;jh4k1", "MKLN_K2": "k2", "MKLN_K3": "k3", "MKLN_K4": "k4"}
# w4-9 (filed with WT4) is identical, all 7 mutations over 520 bp, to two MKLN K2 clones; the
# author removed it from WT4 on 12 Oct 2026, and the manifest now files it with MKLN_K2
NEVER_FILED = {("432625", "E6"), ("414309", "F7")}


def main(manifest):
    rows = list(csv.DictReader(open(Path(manifest).expanduser())))
    out = []
    for r in rows:
        if not r["original_path"].startswith("/"):
            continue                     # GenBank-only records: no trace
        if r["status"] in ("kept", "relabelled"):
            mouse = r["mouse"]
        else:
            mouse = FOLDER.get(r["genbank_folder"]) or (
                r["mouse_by_label"] if r["mouse_by_label"] in STUDY | MKLN else "")
        if not mouse:
            continue
        exclude = ""
        if (r["order"], r["well"]) in NEVER_FILED:
            exclude = "never filed with a mouse; AID-type mutations under a knockout label"
        # a tube label naming another mouse: kept where the published analysis filed it
        toks = set(re.split(r"[^0-9a-z]+", r["original_label"].lower()))
        named = [m for m, al in ALIASES.items() if m != mouse and toks & set(al.split(";"))]
        note = ""
        if named and not toks & set(ALIASES[mouse].split(";")):
            note = (f"label names {named[0]}; kept with {mouse} as classified in the original "
                    "publication (author's decision, 12 Oct 2026)")
        out.append({"file": r["original_path"].replace(HOME, "~"), "sample": mouse,
                    "group": (STUDY | MKLN)[mouse], "control": "yes" if mouse.startswith("Aicda") else "no",
                    "order": r["order"], "well": r["well"], "label": r["original_label"],
                    "exclude": exclude, "aliases": ALIASES[mouse], "note": note})
    out.sort(key=lambda x: (x["sample"], x["order"], x["well"]))
    here = Path(__file__).parent
    for name, keep in (("jh4_manuscript.csv", STUDY), ("jh4_with_mkln.csv", STUDY | MKLN)):
        sel = [x for x in out if x["sample"] in keep]
        with open(here / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sel[0]))
            w.writeheader()
            w.writerows(sel)
        print(name, len(sel), "rows")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "~/Nanoshm/sanger_ab1_original/manifest.csv")
