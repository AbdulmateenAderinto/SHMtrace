"""The expected-error / length gate and deterministic downsampling."""
from nanoshm.config import Params
from nanoshm.reads_filter import downsample, prefilter

Q30 = chr(30 + 33)
Q10 = chr(10 + 33)


def _read(name, n=100, q=Q30):
    return (name, "A" * n, q * n)


def test_clean_reads_pass():
    kept, stats = prefilter([_read("a"), _read("b")], Params())
    assert stats["passed_prefilter"] == 2
    assert stats["dropped_ee"] == 0 and stats["dropped_length"] == 0


def test_noisy_reads_are_dropped_on_expected_error():
    kept, stats = prefilter([_read("clean"), _read("noisy", q=Q10)], Params())
    assert [r[0] for r in kept] == ["clean"]
    assert stats["dropped_ee"] == 1


def test_overlong_reads_are_dropped_as_concatemers():
    p = Params(max_total_len=200)
    kept, stats = prefilter([_read("ok", 150), _read("chimera", 500)], p)
    assert [r[0] for r in kept] == ["ok"]
    assert stats["dropped_length"] == 1


def test_length_gate_is_checked_before_quality():
    """An over-length read counts as length-dropped, not double-counted."""
    p = Params(max_total_len=200)
    _kept, stats = prefilter([_read("both", 500, Q10)], p)
    assert stats["dropped_length"] == 1
    assert stats["dropped_ee"] == 0


def test_stats_account_for_every_input_read():
    p = Params(max_total_len=200)
    reads = [_read("a"), _read("b", q=Q10), _read("c", 500)]
    kept, s = prefilter(reads, p)
    assert s["reads_in"] == 3
    assert s["passed_prefilter"] + s["dropped_ee"] + s["dropped_length"] == s["reads_in"]
    assert len(kept) == s["passed_prefilter"]


def test_downsampling_is_deterministic_for_a_seed():
    reads = [_read(f"r{i}") for i in range(50)]
    p = Params(downsample_to=10, downsample_seed=7)
    a, _ = downsample(reads, p)
    b, _ = downsample(reads, p)
    assert [r[0] for r in a] == [r[0] for r in b]


def test_different_seeds_give_different_samples():
    reads = [_read(f"r{i}") for i in range(50)]
    a, _ = downsample(reads, Params(downsample_to=10, downsample_seed=1))
    b, _ = downsample(reads, Params(downsample_to=10, downsample_seed=2))
    assert [r[0] for r in a] != [r[0] for r in b]


def test_downsampling_is_a_noop_when_shallower_than_target():
    reads = [_read(f"r{i}") for i in range(5)]
    out, stats = downsample(reads, Params(downsample_to=500))
    assert len(out) == 5
    assert stats["downsample_applied"] is False


def test_downsampling_disabled_keeps_everything():
    reads = [_read(f"r{i}") for i in range(50)]
    out, stats = downsample(reads, Params(downsample_to=None))
    assert len(out) == 50
    assert stats["downsample_applied"] is False
