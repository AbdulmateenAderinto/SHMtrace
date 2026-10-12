"""End-to-end behaviour on the bundled example.

These are the tests that matter scientifically: the tool must keep recovering
the AID-dependent difference between wild-type and Aicda-/- reads.
"""
import json

import pytest

from shmtrace import pipeline
from shmtrace.config import load_preset

pytest.importorskip("mappy", reason="alignment backend not installed")


@pytest.fixture(scope="module")
def wt_result(wt_fastq, reference, contaminant):
    p = load_preset("jh4_intron")
    return pipeline.run(str(wt_fastq), str(reference), [str(contaminant)], p)


@pytest.fixture(scope="module")
def ko_result(ko_fastq, reference, contaminant):
    p = load_preset("jh4_intron")
    return pipeline.run(str(ko_fastq), str(reference), [str(contaminant)], p)


def test_summary_has_the_documented_shape(wt_result):
    s = wt_result["summary"]
    for key in ("params", "reads", "calling", "unique_variant_density",
                "spectrum_unique", "spectrum_weighted", "gc_bias_unique"):
        assert key in s, f"missing summary key {key}"


def test_reads_are_accounted_for(wt_result):
    r = wt_result["summary"]["reads"]
    assert r["reads_in"] > 0
    assert r["passed_prefilter"] <= r["reads_in"]
    assert r["aligned_to_target"] <= r["passed_prefilter"]


def test_wild_type_carries_many_more_variants_than_aid_deficient(wt_result, ko_result):
    """The core biology. If this fails the tool is not measuring SHM."""
    wt = wt_result["summary"]["calling"]["unique_variants"]
    ko = ko_result["summary"]["calling"]["unique_variants"]
    assert wt > ko, f"WT {wt} should exceed Aicda-/- {ko}"
    assert wt > 3 * ko, f"expected a large AID-dependent gap, got WT {wt} vs KO {ko}"


def test_wild_type_density_exceeds_knockout(wt_result, ko_result):
    assert (wt_result["summary"]["unique_variant_density"]
            > ko_result["summary"]["unique_variant_density"])


def test_wild_type_hits_both_at_and_gc_substantially(wt_result):
    """SHM in the JH4 intron is not a G:C-only process.

    AID deaminates C, but error-prone repair by Pol-eta spreads mutations to
    A:T, and in the neutral JH4 intron A:T substitutions are a large fraction of
    the total. So the correct expectation is that both classes are well
    represented, not that G:C dominates.
    """
    g = wt_result["summary"]["gc_bias_unique"]
    total = g["at_GC"] + g["at_AT"]
    assert total > 0
    assert g["at_GC"] / total > 0.15, "G:C substitutions essentially absent"
    assert g["at_AT"] / total > 0.15, "A:T substitutions essentially absent"


def test_wild_type_base_signature_differs_from_the_aid_null_background(wt_result, ko_result):
    """The point of the Aicda-/- control.

    Real SHM and residual nanopore error do not target the same bases, so the
    G:C fraction of wild-type calls should differ from that of the background
    calls. If these ever converge, the background is no longer distinguishable
    from signal and the filtering needs re-deriving.
    """
    wt = wt_result["summary"]["gc_bias_unique"]["fraction_GC"]
    ko = ko_result["summary"]["gc_bias_unique"]["fraction_GC"]
    assert abs(wt - ko) > 0.10, f"WT {wt} and Aicda-/- {ko} signatures are too similar"


def test_bundled_example_regression(wt_result, ko_result):
    """Pins current output so caller changes show up as a deliberate diff.

    Values as of the phasing + Benjamini-Hochberg correction (31 July 2026).
    Before it: WT 322, KO 28. The KO drop is the point, see
    docs/calibration/FINDINGS.md. Update these numbers deliberately.

    October 2026: WT 281 -> 310 because phasing no longer drops a 4-site block
    confined to this mouse (likely an expanded clone; see phasing.py).
    """
    assert wt_result["summary"]["calling"]["unique_variants"] == pytest.approx(310, abs=20)
    assert ko_result["summary"]["calling"]["unique_variants"] == pytest.approx(9, abs=6)


def test_background_is_far_below_signal(wt_result, ko_result):
    """The headline property: the AID-null background must stay well separated."""
    wt = wt_result["summary"]["unique_variant_density"]
    ko = ko_result["summary"]["unique_variant_density"]
    assert wt / ko > 20, f"separation collapsed to {wt / ko:.1f}x"


def test_run_is_deterministic(wt_fastq, reference, contaminant):
    p = load_preset("jh4_intron")
    a = pipeline.run(str(wt_fastq), str(reference), [str(contaminant)], p)
    b = pipeline.run(str(wt_fastq), str(reference), [str(contaminant)], p)
    assert a["summary"]["calling"] == b["summary"]["calling"]
    assert [v.pos for v in a["variants"]] == [v.pos for v in b["variants"]]


def test_variants_respect_the_germline_threshold(wt_result):
    """Anything at or above the germline frequency is reference, not SHM."""
    threshold = wt_result["summary"]["params"]["germline_freq_threshold"]
    assert all(v.frequency < threshold for v in wt_result["variants"])


def test_variants_respect_the_p_value_cutoff(wt_result):
    cutoff = wt_result["summary"]["params"]["max_pvalue"]
    assert all(v.p_value <= cutoff for v in wt_result["variants"])


def test_no_indels_are_called(wt_result):
    assert all(v.ref in "ACGT" and v.alt in "ACGT" for v in wt_result["variants"])


def test_per_read_mutation_loads_are_reported(wt_result):
    loads = wt_result["per_read_mutations"]
    assert len(loads) > 0
    assert all(isinstance(n, int) and n >= 0 for n in loads)


def test_write_outputs_produces_readable_files(wt_result, tmp_path):
    pipeline.write_outputs(wt_result, str(tmp_path))
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["calling"]["unique_variants"] == wt_result["summary"]["calling"]["unique_variants"]
    rows = (tmp_path / "variants.csv").read_text().strip().splitlines()
    assert rows[0].startswith("position,ref,alt")
    assert len(rows) == len(wt_result["variants"]) + 1
