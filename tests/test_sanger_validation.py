"""End to end on simulated Sanger projects with known truth (sanger_sim.py).

Each project has six animals of ~30 clones, unusable traces, and one of every
bookkeeping error the QC is meant to catch. The thresholds are below what 40
simulated projects achieved (docs/calibration/SANGER_PIPELINE.md),
so a change that makes calling or QC worse fails here.
"""
import pytest

pytest.importorskip("mappy", reason="alignment backend not installed")

from shmtrace import sanger_sim  # noqa: E402


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    d = tmp_path_factory.mktemp("sanger_sim")
    return [sanger_sim.run_once(str(d / f"rep{s}"), seed=s) for s in (3, 9)]


def test_calls_are_found_and_germline_is_not_called(runs):
    s = sanger_sim.summarise(runs)
    assert s["recall_mean"] >= 0.97
    assert s["precision_mean_mutated"] >= 0.88
    assert s["germline_called"] == 0
    # the single-pass Sanger error floor at Q20; real Aicda-/- clones: ~30 per 1e5
    assert s["false_calls_per_1e5_bases"] < 45


def test_failed_and_paralogue_traces_are_not_analysed(runs):
    assert all(r["excluded_ok"] for r in runs)


def test_every_planted_bookkeeping_error_is_flagged(runs):
    for r in runs:
        for p in r["planted"]:
            if p["error"] == "misfiled_relabelled" and p["shared_with_mates"] < 6:
                continue          # below the measured detection limit
            assert p["detected"], p


def test_correctly_filed_clones_are_not_flagged(runs):
    assert all(not r["false_flags"] for r in runs)
    assert all(not r["germline_flags_without_lineage"] for r in runs)


def test_mutations_lost_to_the_germline_rule_are_flagged(runs):
    assert all(x["flagged"] for r in runs for x in r["lineage_germline"])


def test_nothing_is_flagged_when_nothing_is_wrong(tmp_path):
    r = sanger_sim.run_once(str(tmp_path / "clean"), seed=5, plant=False)
    assert not r["false_flags"]
