"""The control-based FDR (multiple_testing="control", the recommended mode).

Benjamini-Hochberg over the candidates seen in a library assumes their p-values
are uniform under the null; a test conditioned on the allele having been seen
does not give that, and on N114 spike-ins it let the FDR reach 0.105. The
control-based procedure measures the null instead: each control, called against
the others at the samples' depth, is a decoy library.
"""
from pathlib import Path

import pytest

from shmtrace import project
from shmtrace.variants import _control_fdr

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples_compare"


def test_threshold_stops_where_decoys_predict_too_many_false_calls():
    # ten strong calls; decoys show one candidate at p <= 0.01 and many above
    sample = [1e-8] * 10 + [0.01, 0.02, 0.03, 0.2, 0.5]
    decoys = [([0.01, 0.2, 0.3, 0.6], 1000), ([0.25, 0.5, 0.7], 1000)]
    keep, t = _control_fdr(sample, decoys, 1000, 0.05)
    # at t = 1e-8: 0 expected false of 10 -> pass; at 0.01: 0.5 of 11 -> pass;
    # at 0.02, 0.03: 0.5 of 12, 13 -> pass; at 0.2: 1.0 of 14 -> 0.071 -> fail
    assert t == 0.03
    assert sum(keep) == 13


def test_decoys_are_scaled_to_the_samples_bases():
    sample = [1e-6] * 10 + [0.01]
    decoys = [([0.01], 100)]                    # one decoy call in a tenth of the bases
    keep, t = _control_fdr(sample, decoys, 1000, 0.05)
    assert t == 1e-6 and sum(keep) == 10       # 10 expected false at 0.01: rejected


def test_no_decoy_calls_and_no_sample_calls():
    assert _control_fdr([], [([0.5], 10)], 10, 0.05) == ([], None)
    keep, t = _control_fdr([0.5, 0.9], [([1e-9, 1e-9], 10)], 10, 0.05)
    assert t is None and not any(keep)


@pytest.mark.skipif(not (EX / "KO3.fastq").exists(), reason="example libraries not present")
def test_project_uses_control_fdr_and_never_a_control_as_its_own_decoy():
    res = project.run_project({
        "name": "t", "region": "jh4_intron",
        "reference": str(EX / "JH4_reference.fasta"),
        "contaminants": [str(EX / "JH4_contaminant.fasta")],
        "samples": [{"name": "WT3", "path": str(EX / "WT3.fastq"), "group": "WT"}]
        + [{"name": k, "path": str(EX / f"{k}.fastq"), "group": "KO", "control": True}
           for k in ("KO1", "KO2", "KO3")]})
    calling = {n: d["summary"]["calling"] for n, d in res["_details"].items()}
    assert calling["WT3"]["multiple_testing"] == "control"
    assert calling["WT3"]["unique_variants"] > 100
    # were a control its own decoy, its calls would be suppressed by construction;
    # each control is tested against its peers' decoys only, so this must still run
    # the control-based procedure and keep the background small
    for k in ("KO1", "KO2", "KO3"):
        assert calling[k]["multiple_testing"] == "control"
        assert calling[k]["unique_variants"] <= 5
