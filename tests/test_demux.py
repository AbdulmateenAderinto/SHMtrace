"""Inline-barcode demultiplexing and screen ranking (demux.py)."""
import csv
import random

import pytest

from nanoshm import demux, samplesheet

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
