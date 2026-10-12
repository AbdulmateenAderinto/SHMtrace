"""Platform concordance statistics (concordance.py)."""
import math

import numpy as np
import pytest

from nanoshm import concordance

A = {"m1": {"G>A": 30, "C>T": 25, "A>G": 20, "T>C": 10, "A>T": 8, "G>C": 7},
     "m2": {"G>A": 40, "C>T": 30, "A>G": 25, "T>C": 12, "A>T": 9, "G>C": 4},
     "m3": {"G>A": 20, "C>T": 18, "A>G": 15, "T>C": 9, "A>T": 6, "G>C": 3}}


def test_identical_spectra_agree_perfectly():
    s = concordance.spectrum(A, A)
    assert s["pearson_r"] == pytest.approx(1)
    assert s["slope"] == pytest.approx(1)
    assert s["max_abs_diff_pp"] == pytest.approx(0)
    assert all(t["p_sidak"] == 1.0 for t in s["paired"]["per_type"])


def test_wilson_interval_matches_the_textbook_value():
    lo, hi = concordance.wilson(10, 100)
    assert lo == pytest.approx(5.52, abs=0.01) and hi == pytest.approx(17.44, abs=0.01)


def test_a_platform_that_shifts_one_type_is_detected():
    B = {k: dict(v) for k, v in A.items()}
    for v in B.values():
        v["A>T"] = v["A>T"] + 25            # B finds many more A>T in every animal
    s = concordance.spectrum(A, B)
    biggest = max(s["paired"]["per_type"], key=lambda t: abs(t["mean_diff_pp"]))
    assert biggest["type"] == "A>T" and biggest["mean_diff_pp"] < 0
    assert s["paired"]["interaction_p"] < 0.05


def test_rm_interaction_is_zero_for_additive_platform_offsets():
    a = np.array([[10.0, 20, 30], [12, 22, 32], [9, 19, 29]])
    f, d1, d2, p = concordance.rm_interaction(a, a + 0.0)
    assert (d1, d2) == (2, 4)
    assert math.isnan(f) or f == pytest.approx(0)


def test_positional_profile_compares_the_same_bins():
    p = concordance.positional([110, 112, 310, 311, 312], [111, 113, 315, 316, 317], 101, 400, 25)
    assert p["pearson_r"] == pytest.approx(1)
    assert len(p["starts"]) == 12
