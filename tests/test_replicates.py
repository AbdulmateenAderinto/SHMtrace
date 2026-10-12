"""Group tests with the animal as the unit (replicates.py), and power (power.py)."""
import numpy as np
import pytest

from nanoshm import power, replicates

N = 262913


def _animals(rates, rng, size=5.0):
    lam = np.asarray(rates) * N
    return [(int(k), N) for k in rng.negative_binomial(size, size / (size + lam))]


def test_negative_binomial_recovers_group_rates():
    rng = np.random.default_rng(1)
    a = _animals([3e-3] * 8, rng)
    b = _animals([6e-3] * 8, rng)
    (ra, rb), phi, _ll = replicates.fit_nb([tuple(zip(*a, strict=True)), tuple(zip(*b, strict=True))])
    assert ra == pytest.approx(3e-3, rel=0.25) and rb == pytest.approx(6e-3, rel=0.25)


def test_a_large_difference_is_found_and_pooling_reads_overstates_it():
    rng = np.random.default_rng(2)
    ctrl = _animals([5e-5, 8e-5, 3e-5], rng)
    case = _animals([3e-3, 9e-3, 1.5e-3], rng)
    c = replicates.compare(ctrl, case, boot=199)
    assert c["fold"] > 30 and c["fold_ci"][0] > 1
    assert c["p_negbinomial"] < 0.05 and c["p_welch_log"] < 0.05
    assert c["p_permutation"] == pytest.approx(0.1)          # the floor for 3 vs 3
    assert c["p_pooled_reads"] < 1e-10                        # pseudoreplication


def test_identical_groups_are_not_called_different():
    rng = np.random.default_rng(3)
    ps = []
    for _ in range(40):
        c = replicates.compare(_animals([3e-3] * 4, rng), _animals([3e-3] * 4, rng), boot=99)
        ps.append(c["p_negbinomial"])
    assert np.mean(np.array(ps) < 0.05) <= 0.15


def test_permutation_is_exact_for_small_groups():
    assert replicates.permutation_p([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) == pytest.approx(0.1)


def test_power_grows_with_animals_not_reads_when_animals_vary():
    lo = power.simulate(3e-3, 2.0, 3, 500, 529, 1.0, 5.0, sims=400)
    hi = power.simulate(3e-3, 2.0, 12, 500, 529, 1.0, 5.0, sims=400)
    deep = power.simulate(3e-3, 2.0, 3, 4000, 529, 1.0, 5.0, sims=400)
    assert hi > lo + 0.2
    assert abs(deep - lo) < 0.1


def test_pilot_reads_a_project_result():
    res = {"samples": [{"sample": f"W{i}", "group": "WT", "mutations_total": k,
                        "frequency_bases": N, "unique_mutations": 50, "density_bases": N,
                        "analysed_bp": 529, "reads": 497} for i, k in enumerate((500, 1000, 2000))],
           "mutation_load": {f"W{i}": {"0": 300, "1": 100, "5": 50, "20": 47} for i in range(3)}}
    pl = power.pilot(res, "WT")
    assert pl["rate"] == pytest.approx(np.mean([500, 1000, 2000]) / N)
    assert 0.5 < pl["cv_between"] < 1.0 and pl["dispersion"] > 1
