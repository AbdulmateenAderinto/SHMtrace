"""Inline-barcode demultiplexing and screen ranking (demux.py)."""
import csv
import random

import pytest

from shmtrace import demux, samplesheet

rng = random.Random(5)
REF = "".join(rng.choice("ACGT") for _ in range(600))
BCS = {"S1": "ACGTACGTAC", "S2": "TTGCAAGCTG", "S3": "GGATCCTAGA", "S4": "CATGGTCAAC"}
BCS2 = {"S1": "AGCTTCGACA", "S2": "CTAGGATCGT", "S3": "TCCAGGTTAC", "S4": "GTGACCATGA"}


def _noisy(s, rate=0.03):
    out = []
    for b in s:
        u = rng.random()
        if u < rate / 3:
            continue
        if u < 2 * rate / 3:
            out.append(rng.choice("ACGT"))
        out.append(rng.choice([x for x in "ACGT" if x != b]) if rng.random() < rate / 3 else b)
    return "".join(out)


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _pool(path, n=60, dual=False):
    truth = {}
    with open(path, "w") as fh:
        for name in BCS:
            for i in range(n):
                s = _noisy("TTAG" + BCS[name] + REF + (BCS2[name] if dual else "") + "GATC")
                if rng.random() < 0.5:
                    s = _rc(s)
                rid = f"{name}_{i}"
                truth[rid] = name
                fh.write(f"@{rid}\n{s}\n+\n{'5' * len(s)}\n")
    return truth


@pytest.mark.parametrize("dual", [False, True])
def test_reads_go_to_their_sample_and_doubtful_ones_to_none(tmp_path, dual):
    truth = _pool(tmp_path / "pool.fq", dual=dual)
    samples = {n: (BCS[n], BCS2[n] if dual else None) for n in BCS}
    d = demux.demultiplex(str(tmp_path / "pool.fq"), samples, REF[:20], REF[-20:],
                          str(tmp_path / "out"))
    wrong = right = 0
    for n, path in d["paths"].items():
        for line in open(path).read().split("\n")[::4]:
            if line.startswith("@"):
                if truth[line[1:]] == n:
                    right += 1
                else:
                    wrong += 1
    assert right >= 0.9 * len(truth)
    assert wrong <= 0.01 * len(truth)


def test_a_sample_sheet_with_barcodes_is_demultiplexed(tmp_path):
    _pool(tmp_path / "pool.fq")
    (tmp_path / "ref.fa").write_text(f">ref\n{REF}\n")
    with open(tmp_path / "sheet.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "sample", "group", "control", "barcode"])
        for n, bc in BCS.items():
            w.writerow(["pool.fq", n, "C" if n == "S1" else "T", "yes" if n == "S1" else "no", bc])
    sh = samplesheet.load(str(tmp_path / "sheet.csv"))
    proj = sh.project(str(tmp_path / "ref.fa"), "generic_amplicon", platform="ont")
    assert all(s["path"].endswith(f"{s['name']}.fastq") for s in proj["samples"])
    assert sum(proj["demux"]["pool.fq"]["per_sample"].values()) >= 0.9 * 4 * 60


def test_screen_ranks_conditions_with_bh_q_values():
    res = {"animal_tests": [
        {"metric": "unique_density", "group": g, "control": "EV", "animals_case": 3,
         "fold": f, "fold_ci": None, "p_negbinomial": p}
        for g, f, p in (("A", 4.0, 0.001), ("B", 1.1, 0.8), ("C", 0.3, 0.004), ("D", 1.5, 0.2))]}
    rows = demux.screen_table(res)
    assert [r["condition"] for r in rows[:2]] == ["A", "C"]
    assert rows[0]["q"] == pytest.approx(0.004) and rows[1]["q"] == pytest.approx(0.008)


def test_bh_matches_the_textbook():
    assert demux.bh([0.01, 0.04, 0.03, 0.005]) == pytest.approx([0.02, 0.04, 0.04, 0.02])


def test_the_app_splits_a_pooled_run_and_ranks_a_screen(tmp_path):
    """Sheet with barcodes -> /api/run splits the pool -> /api/screen ranks the conditions."""
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    pytest.importorskip("mappy")
    import time

    from fastapi.testclient import TestClient

    from shmtrace.webapp.main import app
    r = random.Random(8)
    ref = "".join(r.choice("ACGT") for _ in range(500))
    (tmp_path / "ref.fa").write_text(f">amp\n{ref}\n")
    loads = {"EV": 0.0, "up": 3.0}
    rows, codes = [], {}
    with open(tmp_path / "pool.fq", "w") as fh:
        for cond, load in loads.items():
            for rep in (1, 2, 3):
                name = f"{cond}_{rep}"
                codes[name] = "".join(r.choice("ACGT") for _ in range(10))
                rows.append(["pool.fq", name, cond, "yes" if cond == "EV" else "no", codes[name]])
                clones = []
                for _ in range(30):
                    m = {}
                    for _ in range(int(r.expovariate(1 / load)) if load else 0):
                        i = r.randrange(40, 460)
                        m[i] = r.choice([x for x in "ACGT" if x != ref[i]])
                    clones.append(m)
                for k in range(150):
                    s = list(ref)
                    for i, a in r.choice(clones).items():
                        s[i] = a
                    seq = _noisy(codes[name] + "".join(s), 0.01)
                    if r.random() < 0.5:
                        seq = _rc(seq)
                    fh.write(f"@{name}_{k}\n{seq}\n+\n{'?' * len(seq)}\n")
    with open(tmp_path / "sheet.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "sample", "group", "control", "barcode"])
        w.writerows(rows)
    c = TestClient(app)
    sh = c.post("/api/samplesheet", json={"path": str(tmp_path / "sheet.csv")}).json()
    assert len(next(iter(sh["pooled"].values()))) == 6
    job = c.post("/api/run", json={"project": {
        "name": "screen", "region": "generic_amplicon", "reference": str(tmp_path / "ref.fa"),
        "samples": [{k: s[k] for k in ("name", "path", "group", "control")} for s in sh["samples"]],
        "sample_sheet": sh["sheet"]}}).json()["job"]
    for _ in range(600):
        st = c.get(f"/api/job/{job}").json()
        if st["state"] in ("done", "error"):
            break
        time.sleep(0.2)
    assert st["state"] == "done", st
    res = c.get(f"/api/result/{job}").json()
    assert "split by barcode" in res["notes"][0]
    assert all(s["reads"] >= 100 for s in res["samples"])
    sc = c.post(f"/api/screen/{job}", json={"metric": "unique_density"}).json()
    up = sc["rows"][0]
    assert up["condition"] == "up" and up["fold"] > 5
    assert c.get(f"/api/screen/{job}/screen.csv").text.startswith("condition,")
