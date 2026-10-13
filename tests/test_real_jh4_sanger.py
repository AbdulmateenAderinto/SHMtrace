"""Regression: the lab's real JH4 Sanger traces give the manuscript's numbers (Ver17).

The traces are not in the repository; the test runs where they are (the sample
sheet's paths) and is skipped elsewhere. A change to calling, trimming,
alignment or the germline rule that moves a published number fails here.
"""
import os
from pathlib import Path

import pytest

pytest.importorskip("mappy", reason="alignment backend not installed")

from shmtrace import project, samplesheet  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SHEET = ROOT / "docs/calibration/samplesheets/jh4_manuscript.csv"
REF = ROOT / "src/shmtrace/references"
# mouse: (clones analysed, unique calls in 101-629) as in the manuscript
EXPECTED = {"WT2": (37, 135), "WT3": (40, 230), "WT4": (26, 106),
            "Aicda_A1": (11, 1), "Aicda_A2": (10, 1), "Aicda_A3": (10, 2)}


def _first_file():
    import csv
    row = next(csv.DictReader(open(SHEET)))
    return os.path.expanduser(row["file"])


pytestmark = pytest.mark.skipif(not SHEET.exists() or not os.path.exists(_first_file()),
                                reason="the lab's JH4 Sanger traces are not on this machine")


@pytest.fixture(scope="module")
def result():
    sh = samplesheet.load(str(SHEET))
    return project.run_project(sh.project(str(REF / "mouse_jh4_intron.fasta"), "jh4_intron",
                                          contaminants=[str(REF / "mouse_jh4_contaminant.fasta")]))


def test_manuscript_clone_and_call_counts(result):
    got = {n: (len(d["_aligned"]), sum(101 <= v.pos <= 629 for v in d["variants"]))
           for n, d in result["_details"].items()}
    assert got == EXPECTED


def test_qc_finds_the_known_bookkeeping_problems(result):
    flags = result["qc"]["flags"]
    dup = {f["clone"] for f in flags if f["check"] == "duplicate_clone"}
    # wt3-3, filed with WT2, identical to WT3 clone 434141-F9 (wt3-4 and wt3-6 were moved
    # to WT3 on their founder mutations, so they no longer cross animals)
    assert "434138-C11" in dup and not {"434138-D11", "434138-F11"} & dup
    mism = [f for f in flags if f["check"] == "label_mismatch"]
    assert sum(f["sample"] == "WT2" and f["other_sample"] == "WT3" for f in mism) >= 12
    assert sum(f["sample"] == "WT3" and f["other_sample"] == "WT2" for f in mism) >= 13
    # no Aicda-/- clone carries more mutations than its qualities explain
    assert not [f for f in flags if f["check"] == "control_mutations"]
    # the two never-filed contaminants are listed as excluded by the sheet
    assert sum(e["by"] == "sample sheet" for e in result["qc"]["excluded"]) == 2
