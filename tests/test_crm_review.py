"""Regression checks for the scientific review's simulation and window fixes."""
import importlib.util
import random
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

CRM = Path(__file__).resolve().parents[1] / "docs/calibration/crm"
sys.path.insert(0, str(CRM))


def load(name):
    spec = importlib.util.spec_from_file_location(name, CRM / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_shared_artifact_has_identical_alleles_across_libraries(monkeypatch):
    fv = load("fdr_validation")
    monkeypatch.setattr(fv.simulate, "_ref_to_query", lambda hit, n: dict(enumerate(range(n))))
    aligner = SimpleNamespace(map=lambda seq: iter([SimpleNamespace(is_primary=True, strand=1)]))
    reads = [("r", "A" * 12, "I" * 12)]
    alleles = {2: "G", 5: "C"}
    rng = random.Random(7)
    results = [fv.add_artifact(reads, "A" * 12, aligner, [2, 5], 1, rng, alleles)
               for _ in range(4)]
    assert all(result[0][1] == "AAGAACAAAAAA" for result in results)
    assert all(result[0][2] == "I" * 12 for result in results)


def test_shared_artifact_rejects_reference_alleles():
    fv = load("fdr_validation")
    with pytest.raises(ValueError, match="non-reference"):
        fv.add_artifact([], "AAAA", None, [2], 1, random.Random(1), {2: "A"})


def test_artifact_injection_preserves_reverse_read_orientation_and_quality(monkeypatch):
    fv = load("fdr_validation")
    monkeypatch.setattr(fv.simulate, "_ref_to_query", lambda hit, n: dict(enumerate(range(n))))
    aligner = SimpleNamespace(map=lambda seq: iter([SimpleNamespace(is_primary=True, strand=-1)]))
    ref = "AACCGGTTACGA"
    reads = [("reverse", "TCGTAACCGGTT", "ABCDEFGHIJKL")]
    unchanged = fv.add_artifact(reads, ref, aligner, [2], 0, random.Random(7), {2: "A"})
    assert unchanged == reads
    changed = fv.add_artifact(reads, ref, aligner, [2], 1, random.Random(7), {2: "A"})
    assert changed[0][2] == reads[0][2]
    forward = changed[0][1].translate(fv._COMP)[::-1]
    assert forward[2] == "A"
    assert sum(a != b for a, b in zip(forward, ref, strict=True)) == 1


def test_two_control_design_uses_both_source_libraries():
    fv = load("fdr_validation")
    controls = dict.fromkeys(["KO2a", "KO2b", "KO3a", "KO3b"])
    assert fv.select_controls(controls, 2) == ["KO2a", "KO3a"]
    assert set(fv.select_controls(controls, 3)) == {"KO2a", "KO2b", "KO3a"}
    with pytest.raises(ValueError):
        fv.select_controls(controls, 5)


def test_benchmark_window_is_zero_based_and_includes_both_boundaries():
    bench = load("benchmark_real")
    read = SimpleNamespace(pairs=[(99, "A", 30), (100, "C", 30), (628, "G", 30),
                                  (629, "T", 30), (200, "N", 30), (201, None, 30)])
    assert bench.window_bases([read], 100, 628) == 2


def test_benchmark_fitted_null_carries_decoys_into_control_threshold(tmp_path):
    from shmtrace.align import AlignedRead
    from shmtrace.config import Params
    from shmtrace.null_model import fit_aligned
    from shmtrace.variants import call_variants

    ref = "ACGTACGTACGT"
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">amp\n{ref}\n")
    reads = [[AlignedRead(f"k{k}_{i}", "amp", 0, len(ref),
                          [(j, b, 40) for j, b in enumerate(ref)], ref)
              for i in range(4)] for k in range(2)]
    p = Params(error_model="qcal", multiple_testing="control", analysis_window=(1, 12),
               end_trim_bp=0, recalibrate=False)
    null, _stats = fit_aligned(reads, ref, p, 0, 11, calibrate=True, ref_name="amp")
    assert len(null.decoys) == 2
    _calls, stats = call_variants(reads[0], str(fasta), "amp", p, null)
    assert stats["multiple_testing"] == "control"


def test_high_fraction_mask_is_per_allele_not_combined_nonreference(tmp_path):
    from shmtrace.align import AlignedRead
    from shmtrace.config import Params
    from shmtrace.variants import call_variants

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">amp\nA\n")
    reads = [AlignedRead(str(i), "amp", 0, 1, [(0, b, 40)], b)
             for i, b in enumerate("AAACCCGGGG")]
    p = Params(error_model="geneious", multiple_testing="none", max_pvalue=1,
               analysis_window=(1, 1), end_trim_bp=0, germline_freq_threshold=0.5)
    calls, stats = call_variants(reads, str(fasta), "amp", p)
    assert {v.alt for v in calls} == {"C", "G"}
    assert stats["germline_positions_removed"] == 0
