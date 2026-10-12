"""Benjamini-Hochberg and the empirical binomial null."""

from shmtrace.variants import _benjamini_hochberg, _binom_sf


def test_bh_rejects_nothing_when_all_null():
    assert not any(_benjamini_hochberg([0.2, 0.4, 0.6, 0.8, 0.99], 0.05))


def test_bh_keeps_a_clear_signal():
    keep = _benjamini_hochberg([1e-12, 0.4, 0.6, 0.8, 0.99], 0.05)
    assert keep[0] and not any(keep[1:])


def test_bh_rejects_essentially_nothing_under_the_null():
    """The whole point.

    Under a true null the p-values are uniform. Uncorrected p <= 0.05 across the
    1,587 tests the caller runs per sample rejects ~79 of them by chance, which
    is the same order as the entire observed Aicda-/- background. BH rejects
    almost none.
    """
    import random
    random.seed(0)
    pvals = [random.random() for _ in range(1587)]
    uncorrected = sum(p <= 0.05 for p in pvals)
    corrected = sum(_benjamini_hochberg(pvals, 0.05))
    assert 50 < uncorrected < 120, uncorrected
    assert corrected <= 1, corrected


def test_bh_step_up_keeps_the_whole_prefix():
    """If the k-th smallest passes, everything smaller is kept too."""
    pvals = [0.001, 0.002, 0.003, 0.9, 0.95]
    keep = _benjamini_hochberg(pvals, 0.05)
    assert keep[:3] == [True, True, True]
    assert keep[3:] == [False, False]


def test_bh_empty_input():
    assert _benjamini_hochberg([], 0.05) == []


def test_binom_sf_is_a_survival_function():
    assert _binom_sf(0, 500, 0.001) == 1.0
    assert _binom_sf(1, 500, 0.001) > _binom_sf(5, 500, 0.001)
    assert _binom_sf(50, 500, 0.001) < 1e-10


def test_binom_sf_bounded():
    for k in (1, 5, 20):
        assert 0.0 <= _binom_sf(k, 500, 0.002) <= 1.0


def test_empirical_null_is_stricter_where_error_is_known():
    """A position known to err at 2% should not call a 2% observation."""
    observed, depth = 10, 500
    assert _binom_sf(observed, depth, 0.02) > 0.05     # matches known error -> not called
    assert _binom_sf(observed, depth, 0.0005) < 0.05   # clean position -> called


class TestEmpiricalNullFloor:
    """An unusually clean control must not make the caller more permissive.

    Fitting the null on near-zero-mutation controls drives the estimated rate
    below the sequencer's own quality-implied error. Without a floor the
    binomial test then becomes *easier* to pass than the read-quality model,
    which is how a cross-run null added calls instead of removing them.
    """

    def test_floor_keeps_the_stricter_of_the_two(self):
        from shmtrace.variants import _binom_sf

        depth, observed = 500, 3
        quality_implied = 0.002          # what the base qualities say
        over_clean_null = 0.00008        # what a pristine control would fit
        p_optimistic = _binom_sf(observed, depth, over_clean_null)
        p_floored = _binom_sf(observed, depth, max(over_clean_null, quality_implied))
        assert p_optimistic < 0.05, "without a floor this call passes"
        assert p_floored > 0.05, "with the floor it does not"

    def test_floor_does_not_block_a_genuinely_noisy_position(self):
        """Where the control shows real error, the empirical rate must win."""
        from shmtrace.variants import _binom_sf

        depth, observed = 500, 3
        quality_implied = 0.0005
        measured_noisy = 0.02            # this position really does err
        assert _binom_sf(observed, depth, max(measured_noisy, quality_implied)) > 0.05
