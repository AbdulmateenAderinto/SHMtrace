"""A published Sanger study re-analysed automatically: J_H4 clones of Barbulescu et al. (2024), cohort 2.

    python docs/calibration/crm/sanger_cohort2.py

Wild-type, Maea^fl/fl Mb1cre and Aicda-/- mice, 375 chromatograms in four 96-well plates and a pilot plate.
Each clone is assigned to the mouse it was filed under in the original analysis project
(samplesheets/cohort2_barbulescu2024.csv). Four traces were filed under two mice; they are listed twice and
analysed as filed, through temporary copies, because SHMtrace accepts a file in one sample only.
Writes docs/calibration/crm/out/sanger_cohort2.json: the project result, the run time, and per-mouse values
in the units of the published Fig. 2g-i, for crm/figs/figS9_sanger_cohort2.py and the manuscript numbers.
"""
import csv
import json
import os
import shutil
import tempfile
import time
from collections import Counter
from pathlib import Path

from shmtrace import project, samplesheet

CAL = Path(__file__).resolve().parents[1]
REF = CAL.parents[1] / "src/shmtrace/references"
SHEET = CAL / "samplesheets/cohort2_barbulescu2024.csv"
OUT = Path(__file__).parent / "out/sanger_cohort2.json"
BINS = [(0, 0), (1, 1), (2, 5), (6, 10), (11, 20), (21, 10 ** 6)]
TS, TV = {"C>T", "G>A"}, {"C>A", "C>G", "G>C", "G>T"}


def main():
    with tempfile.TemporaryDirectory() as tmp:
        rows = list(csv.DictReader(open(SHEET)))
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
        per[m] = {"group": s["group"], "clones": s["reads"], "freq_x1e4": 10 * (s["mutation_frequency_x1e3"] or 0),
                  "load_bins": [sum(v for k, v in load.items() if lo <= k <= hi) for lo, hi in BINS],
                  "spectrum": sp,
                  "gc_ts_pct": 100 * sum(sp.get(k, 0) for k in TS) / tot if tot else None,
                  "gc_tv_pct": 100 * sum(sp.get(k, 0) for k in TV) / tot if tot else None}
    qc = Counter(f["check"] for f in pub["qc"]["flags"])
    json.dump({"seconds": seconds, "traces": len(rows), "clones_analysed": pub["qc"]["clones_analysed"],
               "excluded": Counter(e["reason"] for e in pub["qc"]["excluded"]), "qc_flags": qc,
               "per_mouse": per, "project": pub}, OUT.open("w"), default=str)
    print(f"{len(rows)} traces, {pub['qc']['clones_analysed']} clones analysed in {seconds:.0f} s; QC flags {dict(qc)}")


if __name__ == "__main__":
    main()
