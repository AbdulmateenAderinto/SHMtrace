"""The bundled Sanger example (examples_sanger): Barbulescu et al. 2024's JH4 clones load and give the curated result."""
import math
from pathlib import Path
from statistics import mean

import pytest

from shmtrace import project, samplesheet

ROOT = Path(__file__).resolve().parents[1]
SHEET = ROOT / "examples_sanger" / "samples.csv"
REF = ROOT / "src" / "shmtrace" / "references"
pytestmark = pytest.mark.skipif(not SHEET.is_file(), reason="no Sanger example")


def test_the_sheet_names_twelve_mice_and_every_trace_is_there():
    sh = samplesheet.load(str(SHEET))
    names = {s["name"]: s for s in sh.samples()}
    assert len(names) == 12
    assert {n for n, s in names.items() if s["control"]} == {"Aicda_2a1", "Aicda_2a2", "Aicda_2a3", "Aicda_2a4"}
    assert all(Path(f).is_file() for s in sh.samples() for f, _clone in s["path"])
    assert sum(r.shared_trace for r in sh.used) == 8


def test_the_example_reproduces_the_curated_analysis():
    sh = samplesheet.load(str(SHEET))
    res = project.run_project(sh.project(str(REF / "mouse_jh4_intron.fasta"), "jh4_intron",
                                         contaminants=[str(REF / "mouse_jh4_contaminant.fasta")]))
    assert res["qc"]["clones_analysed"] == 290 and len(res["qc"]["excluded"]) == 3
    f = {g: [s["mutation_frequency_x1e3"] for s in res["samples"] if s["group"] == g] for g in ("WT", "Maea_KO")}
    gm = {g: 10 ** mean(math.log10(v) for v in f[g]) for g in f}
    assert gm["WT"] / gm["Maea_KO"] == pytest.approx(4.66, abs=0.01)
    assert max(s["mutation_frequency_x1e3"] for s in res["samples"] if s["control"]) < 0.5
