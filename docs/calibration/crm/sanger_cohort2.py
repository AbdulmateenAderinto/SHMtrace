"""A published Sanger study re-analysed automatically, from one sample sheet.

    python docs/calibration/crm/sanger_cohort2.py [cohort2|mkln]

cohort2 (default): Barbulescu et al. (2024), J_H4, cohort 2.

Wild-type, Maea^fl/fl Mb1cre and Aicda-/- mice, 375 chromatograms in four 96-well plates and a pilot plate.

mkln: Barbulescu et al. (2025, J Immunol), J_H4 of wild-type, Mkln1-/- and Aicda-/- mice, whose clones are
filed in per-mouse folders of the original analysis project.

In both, each clone is assigned to the mouse it was filed under in that project, through a sample sheet.
A trace filed under two mice is listed twice and analysed as filed, through temporary copies, because
SHMtrace accepts a file in one sample only. Writes out/sanger_<key>.json: the project result, the run time,
and per-mouse values in the units of the published figures, for the figure scripts and the manuscript
numbers.
"""
import csv
import json
import os
import shutil
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from shmtrace import project, samplesheet

CAL = Path(__file__).resolve().parents[1]
REF = CAL.parents[1] / "src/shmtrace/references"
STUDIES = {"cohort2": "samplesheets/cohort2_barbulescu2024.csv",
           "mkln": "samplesheets/jh4_mkln_barbulescu2025.csv"}
BINS = [(0, 0), (1, 1), (2, 5), (6, 10), (11, 20), (21, 10 ** 6)]
TS, TV = {"C>T", "G>A"}, {"C>A", "C>G", "G>C", "G>T"}


def run(key="cohort2"):
    sheet_path = CAL / STUDIES[key]
    out_path = Path(__file__).parent / f"out/sanger_{key}.json"
    with tempfile.TemporaryDirectory() as tmp:
        rows = [r for r in csv.DictReader(open(sheet_path)) if not r["exclude"]]
        seen = Counter()
        for r in rows:
            f = os.path.expanduser(r["file"])
            r["file"] = f
            seen[f] += 1
            if seen[f] > 1:                                  # filed twice: analyse as filed, from a copy
                cp = Path(tmp) / f"{r['clone_id']}_{Path(f).name}"
                shutil.copyfile(f, cp)
                r["file"] = str(cp)
        sheet = Path(tmp) / "sheet.csv"
        with open(sheet, "w", newline="") as fh:            # closed before SHMtrace reads it
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        t0 = time.time()
        sh = samplesheet.load(str(sheet))
        proj = sh.project(str(REF / "mouse_jh4_intron.fasta"), "jh4_intron",
                          contaminants=[str(REF / "mouse_jh4_contaminant.fasta")])
        res = project.run_project(proj)
        seconds = time.time() - t0
    pub = project.public(res)
    per = {}
    for s in pub["samples"]:
        m = s["sample"]
        load = {int(k): v for k, v in pub["mutation_load"].get(m, {}).items()}
        sp = pub["spectrum_sample"].get(m, {})
        tot = sum(sp.values())
        # both papers report unique mutations over nucleotides sequenced, so that is the comparable
        # figure; mutation_frequency_x1e3 counts every clone's mutations, including recurrences
        per[m] = {"group": s["group"], "clones": s["reads"],
                  "freq_x1e4": 1e4 * s["unique_mutations"] / s["frequency_bases"] if s["frequency_bases"] else 0,
                  "freq_all_x1e4": 10 * (s["mutation_frequency_x1e3"] or 0),
                  "unique_mutations": s["unique_mutations"], "bases": s["frequency_bases"],
                  "load_bins": [sum(v for k, v in load.items() if lo <= k <= hi) for lo, hi in BINS],
                  "spectrum": sp,
                  "gc_ts_pct": 100 * sum(sp.get(k, 0) for k in TS) / tot if tot else None,
                  "gc_tv_pct": 100 * sum(sp.get(k, 0) for k in TV) / tot if tot else None}
    qc = Counter(f["check"] for f in pub["qc"]["flags"])
    json.dump({"seconds": seconds, "traces": len(rows), "clones_analysed": pub["qc"]["clones_analysed"],
               "excluded": Counter(e["reason"] for e in pub["qc"]["excluded"]), "qc_flags": qc,
               "per_mouse": per, "project": pub}, out_path.open("w"), default=str)
    print(f"{len(rows)} traces, {pub['qc']['clones_analysed']} clones analysed in {seconds:.0f} s; QC flags {dict(qc)}")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "cohort2")
