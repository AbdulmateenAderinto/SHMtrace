"""Density, spectrum and the G:C bias readout."""
import pytest

from nanoshm.config import Params
from nanoshm.metrics import gc_at_summary, spectrum_matrix, unique_variant_density
from nanoshm.variants import Variant


def v(ref, alt, reads=1, pos=1):
    return Variant(pos=pos, ref=ref, alt=alt, supporting_reads=reads,
                   depth=500, frequency=reads / 500, p_value=1e-6)


def test_density_formula():
    """variants / (bp x depth) x scale"""
    p = Params(density_scale=1e5)
    got = unique_variant_density([v("G", "A")] * 10, analysed_bp=529, mean_depth=500, p=p)
    assert got == pytest.approx(10 / (529 * 500) * 1e5)


def test_density_is_zero_without_analysed_bases():
    assert unique_variant_density([v("G", "A")], 0, 500, Params()) == 0.0


def test_density_scales_linearly_with_variant_count():
    p = Params()
    one = unique_variant_density([v("G", "A")], 529, 500, p)
    ten = unique_variant_density([v("G", "A")] * 10, 529, 500, p)
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
