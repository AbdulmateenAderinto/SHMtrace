"""Checking calls by hand: the Geneious-style tables, the reads behind a call, and review marks."""
import io
import zipfile

import pytest

from shmtrace import browse, geneious, report
from shmtrace.benchmark import read_geneious_export


def _res():
    variants = [
        {"sample": "WT1", "group": "Wt", "pos": 12, "ref": "C", "alt": "T", "reads": 3, "depth": 300,
         "frequency": 0.01, "p_value": 1.2e-6, "plus_reads": 2, "minus_reads": 1, "context": "AGCTA",
         "homopolymer": 1, "WRC_GYW": True},
        {"sample": "WT1", "group": "Wt", "pos": 5, "ref": "A", "alt": "C", "reads": 1, "depth": 290,
         "frequency": 0.00345, "p_value": 0.03, "plus_reads": 1, "minus_reads": 0, "context": "TTACG",
         "homopolymer": 2, "WRC_GYW": False},
        {"sample": "KO1", "group": "KO", "pos": 12, "ref": "C", "alt": "T", "reads": 1, "depth": 310,
         "frequency": 0.0032, "p_value": 0.04, "plus_reads": 1, "minus_reads": 0, "context": "AGCTA",
         "homopolymer": 1, "WRC_GYW": True},
    ]
    return {"samples": [{"sample": "WT1"}, {"sample": "KO1"}], "variants": variants,
            "spectrum_group": {}, "groups": [], "positional": [], "mutation_load": {}}


def test_annotation_rows_are_geneious_shaped_and_ordered_by_position():
    rows = geneious.rows(_res(), "WT1")
    assert [r["Minimum"] for r in rows] == [5, 12]
    r = rows[1]
    assert r["Name"] == "C->T" and r["Change"] == "C -> T" and r["Polymorphism Type"] == "SNP (transition)"
    assert r["Variant Frequency"] == "1.00%" and r["Variant Raw Frequency"] == 3 and r["Coverage"] == 300
    assert rows[0]["Polymorphism Type"] == "SNP (transversion)"


def test_the_export_reads_back_as_the_same_calls(tmp_path):
    path = tmp_path / geneious.file_name("WT1")
    path.write_text(geneious.csv_text(geneious.rows(_res(), "WT1")))
    assert read_geneious_export(str(path), max_p=None, max_freq=None) == {(5, "A", "C"), (12, "C", "T")}


def test_review_marks_reach_the_tables_and_can_drop_rejected_calls():
    res = _res()
    res["review"] = {"WT1|5|A>C": {"mark": "rejected", "note": "one read, end of read"},
                     "WT1|12|C>T": {"mark": "accepted", "note": ""}}
    assert [r["Review"] for r in geneious.rows(res, "WT1", res["review"])] == ["rejected", "accepted"]
    assert [r["Minimum"] for r in geneious.rows(res, "WT1", res["review"], drop_rejected=True)] == [12]
    text = report.tables(res)["variants.csv"]
    assert text.splitlines()[0].endswith("review,review_note")
    assert "one read, end of read" in text
    z = zipfile.ZipFile(io.BytesIO(geneious.zip_bytes(res, drop_rejected=True)))
    assert sorted(z.namelist()) == ["KO1 SHMtrace Annotations.csv", "WT1 SHMtrace Annotations.csv"]
    assert "A -> C" not in z.read("WT1 SHMtrace Annotations.csv").decode()


def test_file_names_survive_awkward_sample_names():
    assert geneious.file_name("Aicda-/- 1") == "Aicda-_- 1 SHMtrace Annotations.csv"


def _alignment():
    reads = [{"name": f"r{i}", "fate": "kept", "start": 1, "end": 30, "strand": 1,
              "diffs": [[12, "T", 35]] if i < 3 else ([[12, "G", 9]] if i == 3 else []), "indels": []}
             for i in range(50)]
    reads.append({"name": "short", "fate": "kept", "start": 15, "end": 30, "strand": 1, "diffs": [], "indels": []})
    reads.append({"name": "gone", "fate": "low quality", "start": 1, "end": 30, "strand": 1,
                  "diffs": [[12, "T", 30]], "indels": []})
    return {"reference_name": "ref", "reference": "ACGT" * 8, "reads": reads,
            "counts": {}, "truncated": 0, "called": [{"pos": 12, "ref": "T", "alt": "T", "reads": 3}]}


def test_call_reads_puts_the_supporting_reads_first_and_counts_only_analysed_reads():
    d = browse.call_reads(_alignment(), pos=12, alt="T", flank=5, others=10)
    assert len(d["carrying"]) == 3 and d["covering"] == 50       # "short" starts after 12, "gone" was dropped
    assert d["others_total"] == 47 and len(d["others"]) == 10
    assert d["other_alt_reads"] == 1                             # the read with G
    assert d["start"] == 7 and d["end"] == 17 and len(d["sequence"]) == 11


fastapi = pytest.importorskip("fastapi")


def test_review_endpoint_records_clears_and_refuses_unknown_calls():
    from fastapi.testclient import TestClient

    from shmtrace.webapp import main
    main.JOBS["t-review"] = {"state": "done", "result": _res(), "progress": 1, "message": "", "error": None}
    c = TestClient(main.app)
    r = c.post("/api/review/t-review", json={"marks": {"WT1|12|C>T": {"mark": "accepted"}}}).json()
    assert r["reviewed"] == 1 and r["accepted"] == 1 and r["calls"] == 3
    assert c.get("/api/review/t-review").json() == {"WT1|12|C>T": {"mark": "accepted", "note": ""}}
    assert c.post("/api/review/t-review", json={"marks": {"WT1|12|C>T": {"mark": ""}}}).json()["reviewed"] == 0
    assert c.post("/api/review/t-review", json={"marks": {"WT1|99|C>T": {"mark": "accepted"}}}).status_code == 400
    assert c.post("/api/review/t-review", json={"marks": {"WT1|12|C>T": {"mark": "maybe"}}}).status_code == 400
    z = c.get("/api/export/t-review/geneious.zip")
    assert z.status_code == 200 and z.headers["content-type"] == "application/zip"
    main.JOBS.pop("t-review")
