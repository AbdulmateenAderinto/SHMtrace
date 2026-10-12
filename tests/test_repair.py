"""Repair-pathway indices and spectrum decomposition (repair.py)."""
import numpy as np
import pytest

from shmtrace import repair

REF = "AAGCTAGCTGACATTGCAGCTAAGCTTGCAAGCTGCA" * 10


def test_indices_read_the_routes_off_the_calls():
    calls = [(3, "G", "A"), (4, "C", "T"), (8, "C", "G"), (10, "G", "C"), (1, "A", "G"), (2, "A", "T")]
    i = repair.indices(calls, REF)
    assert i["mutations"] == 6
    assert i["pct_AT"] == pytest.approx(100 * 2 / 6, abs=0.01)
    assert i["pct_GC_Ts_of_GC"] == 50.0 and i["pct_GC_Tv_of_GC"] == 50.0
    assert i["pct_CG_of_GC_Tv"] == 100.0


def _spec(props, n=4000, rng=None):
    rng = rng or np.random.default_rng(0)
    p = np.array([props.get(k, 0) for k in repair.SUBS12], float)
    p /= p.sum()
    return dict(zip(repair.SUBS12, rng.multinomial(n, p).tolist(), strict=True))


def test_a_mixture_of_reference_spectra_is_recovered():
    ung_msh2 = {"C>T": 1, "G>A": 1}                                     # replication over U only
    msh2 = {"C>T": 2, "G>A": 2, "C>G": 1.5, "G>C": 1.5, "C>A": 0.5, "G>T": 0.5}
    full = {**{k: 0 for k in repair.SUBS12}, "A>G": 3, "T>C": 2, "A>T": 2, "T>A": 1, "A>C": 1,
            "T>G": 0.5}
    mix = {k: 0.2 * ung_msh2.get(k, 0) / 2 + 0.3 * msh2.get(k, 0) / 7.5 + 0.5 * full.get(k, 0) / 9.5
           for k in repair.SUBS12}
    d = repair.decompose(_spec(mix, 20000), {"ung_msh2": _spec(ung_msh2, 20000),
                                             "msh2": _spec(msh2, 20000), "at": _spec(full, 20000)})
    assert d["weights"]["ung_msh2"] == pytest.approx(0.2, abs=0.05)
    assert d["weights"]["msh2"] == pytest.approx(0.3, abs=0.05)
    assert d["weights"]["at"] == pytest.approx(0.5, abs=0.05)
    assert d["r2"] > 0.95


def test_decomposition_intervals_cover_the_truth():
    """95% intervals should cover the true weight in most draws (4 animals, 150 mutations each)."""
    a = {"C>T": 1, "G>A": 1}
    b = {"A>G": 1, "T>C": 1, "A>T": 1}
    mix = {k: 0.4 * a.get(k, 0) / 2 + 0.6 * b.get(k, 0) / 3 for k in repair.SUBS12}
    covered = 0
    for seed in range(8):
        rng = np.random.default_rng(seed)
        animals = [_spec(mix, 150, rng) for _ in range(4)]
        refs = {"a": [_spec(a, 300, rng) for _ in range(3)],
                "b": [_spec(b, 300, rng) for _ in range(3)]}
        d = repair.decompose_animals(animals, refs, boot=200, seed=seed)
        covered += d["ci"]["a"][0] <= 0.4 <= d["ci"]["a"][1]
    assert covered >= 6


def test_the_app_serves_repair_routes_and_power(tmp_path):
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient

    from shmtrace.webapp import main
    calls = {"A1": [(3, "G", "A"), (4, "C", "T"), (8, "C", "G"), (1, "A", "G"), (2, "A", "T")] * 6,
             "A2": [(3, "G", "A"), (8, "C", "G"), (1, "A", "G"), (2, "A", "T")] * 6,
             "B1": [(3, "G", "A"), (4, "C", "T"), (10, "G", "C")] * 8,
             "B2": [(4, "C", "T"), (8, "C", "G"), (10, "G", "C")] * 8}
    group = {"A1": "A", "A2": "A", "B1": "B", "B2": "B"}
    res = {"reference": REF, "group_order": ["A", "B"],
           "samples": [{"sample": n, "group": group[n], "control": False, "mutations_total": 40 + 10 * i,
                        "frequency_bases": 100000, "unique_mutations": 30, "density_bases": 100000,
                        "analysed_bp": 400, "reads": 250} for i, n in enumerate(calls)],
           "variants": [{"sample": n, "pos": p, "ref": r, "alt": a} for n, cs in calls.items()
                        for p, r, a in cs],
           "mutation_load": {n: {"0": 200, "1": 30, "3": 20} for n in calls}}
    main.JOBS["t1"] = {"state": "done", "progress": 1, "message": "", "error": None, "result": res}
    c = TestClient(main.app)
    r = c.post("/api/repair/t1", json={"case": "B", "control": "A", "references": {}}).json()
    assert {x["sample"] for x in r["per_animal"]} == set(calls)
    assert any(x["index"] == "pct_GC_Tv_of_GC" for x in r["comparison"])
    p = c.post("/api/power", json={"job": "t1", "group": "A", "metric": "frequency", "fold": 2,
                                   "animals": [3, 6], "reads": [250], "sims": 200}).json()
    assert len(p["rows"]) == 2 and p["pilot"]["group"] == "A"
    m = c.post("/api/power", json={"rate": 0.003, "cv": 1, "bp": 500, "fold": 4,
                                   "animals": [3, 10], "reads": [500], "sims": 200}).json()
    assert m["rows"][1]["power"] > m["rows"][0]["power"]
