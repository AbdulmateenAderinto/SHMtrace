"""Density, spectrum and the G:C bias readout."""
import pytest

from shmtrace.config import Params
from shmtrace.metrics import gc_at_summary, spectrum_matrix, unique_variant_density
from shmtrace.variants import Variant


def v(ref, alt, reads=1, pos=1):
    return Variant(pos=pos, ref=ref, alt=alt, supporting_reads=reads,
                   depth=500, frequency=reads / 500, p_value=1e-6)


def test_density_formula():
    """variants / bases at callable positions x scale"""
    p = Params(density_scale=1e5)
    got = unique_variant_density([v("G", "A")] * 10, {"analysed_bases": 529 * 500}, p)
    assert got == pytest.approx(10 / (529 * 500) * 1e5)


def test_density_falls_back_to_positions_times_depth_for_older_results():
    p = Params(density_scale=1e5)
    got = unique_variant_density([v("G", "A")] * 10,
                                 {"analysed_positions": 529, "mean_depth": 500}, p)
    assert got == pytest.approx(10 / (529 * 500) * 1e5)


def test_density_is_zero_without_analysed_bases():
    assert unique_variant_density([v("G", "A")], {"analysed_bases": 0}, Params()) == 0.0


def test_density_scales_linearly_with_variant_count():
    p = Params()
    one = unique_variant_density([v("G", "A")], {"analysed_bases": 529 * 500}, p)
    ten = unique_variant_density([v("G", "A")] * 10, {"analysed_bases": 529 * 500}, p)
    assert ten == pytest.approx(10 * one)


def test_spectrum_counts_unique_variants_once():
    mat = spectrum_matrix([v("G", "A", reads=99)], weighted=False)
    assert mat["G"]["A"] == 1


def test_spectrum_weighted_uses_supporting_reads():
    mat = spectrum_matrix([v("G", "A", reads=99)], weighted=True)
    assert mat["G"]["A"] == 99


def test_spectrum_ignores_non_acgt():
    mat = spectrum_matrix([v("N", "A"), v("G", "-")], weighted=False)
    assert sum(sum(r.values()) for r in mat.values()) == 0


def test_gc_fraction_is_the_aid_readout():
    variants = [v("G", "A"), v("C", "T"), v("A", "G")]
    s = gc_at_summary(variants)
    assert s["at_GC"] == 2 and s["at_AT"] == 1
    assert s["fraction_GC"] == pytest.approx(2 / 3, abs=1e-4)


def test_gc_fraction_empty_is_zero_not_a_crash():
    assert gc_at_summary([])["fraction_GC"] == 0.0


def _v(pos, ref="G", alt="A"):
    from shmtrace.variants import Variant
    return Variant(pos, ref, alt, 3, 500, 0.006, 1e-9)


def test_scattered_calls_are_not_flagged_as_clustered():
    import random

    from shmtrace import metrics
    rng = random.Random(0)
    for _ in range(200):
        vs = [_v(p) for p in rng.sample(range(1, 530), 25)]
        assert metrics.positional_clusters(vs, 529) == []


def test_a_tight_run_of_calls_is_flagged():
    from shmtrace import metrics
    vs = [_v(p) for p in (164, 164, 165, 166, 166, 171, 171, 171, 172)]
    vs += [_v(p) for p in (30, 70, 110, 200, 215, 240, 260, 300)]
    hits = metrics.positional_clusters(vs, 529)
    assert len(hits) == 1
    assert (hits[0]["start"], hits[0]["end"], hits[0]["calls"]) == (164, 172, 9)
    assert hits[0]["p_adjusted"] < 1e-3


def test_two_alts_at_one_position_count_as_two_calls():
    from shmtrace import metrics
    one = [_v(100, alt="A"), _v(101), _v(102)]
    two = [_v(100, alt="A"), _v(100, alt="T"), _v(101), _v(102)]
    assert metrics.positional_clusters(one, 529) == []
    assert len(metrics.positional_clusters(two + [_v(100, alt="C"), _v(103)], 529)) == 1


def test_spectrum_note_counts_only_the_expected_classes():
    from shmtrace import metrics
    vs = [_v(10, "C", "T"), _v(20, "G", "A"), _v(30, "G", "C"), _v(40, "A", "T")]
    assert metrics.spectrum_note(vs, ("C>T", "G>A")) == {"calls": 4, "in_expected": 2,
                                                         "fraction": 0.5}
    assert metrics.spectrum_note([], ("C>T",))["fraction"] is None
