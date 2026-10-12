"""The sample sheet: explicit, checked assignment of files to samples."""
import csv

import pytest

from shmtrace import samplesheet
from shmtrace.io_utils import read_fastq, write_abif


def _sheet(path, rows, header=("file", "sample", "group", "control", "order", "well", "label")):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return path


@pytest.fixture
def traces(tmp_path):
    d = tmp_path / "traces"
    d.mkdir()
    for name, seq in (("a", "ACGTACGTAA"), ("b", "ACGTACGTCC"), ("c", "ACGTACGTGG")):
        write_abif(str(d / f"{name}.ab1"), seq, [40] * len(seq))
    return d


def test_samples_are_built_from_rows_and_named_by_order_and_well(tmp_path, traces):
    sh = samplesheet.load(_sheet(tmp_path / "s.csv", [
        ["traces/a.ab1", "WT1", "WT", "no", "100", "A1", "w1-1"],
        ["traces/b.ab1", "WT1", "WT", "no", "100", "B1", "w1-2"],
        ["traces/c.ab1", "KO1", "KO", "yes", "100", "C1", "k1-1"]]))
    s = {x["name"]: x for x in sh.samples()}
    assert s["KO1"]["control"] and not s["WT1"]["control"]
    assert [cid for _f, cid in s["WT1"]["path"]] == ["100-A1", "100-B1"]
    # the reader names each clone by its id, so tables carry order-well
    assert [n for n, _s, _q in read_fastq(s["WT1"]["path"])] == ["100-A1", "100-B1"]
    assert sh.provenance()["100-C1"]["label"] == "k1-1"


def test_every_problem_is_reported_with_its_line(tmp_path, traces):
    with pytest.raises(samplesheet.SheetError) as e:
        samplesheet.load(_sheet(tmp_path / "s.csv", [
            ["traces/a.ab1", "WT1", "WT", "no", "", "", ""],
            ["traces/a.ab1", "WT2", "WT", "no", "", "", ""],        # same file twice
            ["traces/missing.ab1", "WT1", "WT", "no", "", "", ""],  # no such file
            ["traces/b.ab1", "WT1", "KO", "maybe", "", "", ""],     # bad control value
            ["traces/c.ab1", "", "WT", "no", "", "", ""]]))         # no sample
    msg = "\n".join(e.value.problems)
    assert "line 3" in msg and "already listed on line 2" in msg
    assert "line 4: file not found" in msg
    assert "line 5: control is 'maybe'" in msg
    assert "line 6: no sample" in msg


def test_a_sample_must_have_one_group_and_control_status(tmp_path, traces):
    with pytest.raises(samplesheet.SheetError, match="sample WT1 is group 'KO'"):
        samplesheet.load(_sheet(tmp_path / "s.csv", [
            ["traces/a.ab1", "WT1", "WT", "no", "", "", ""],
            ["traces/b.ab1", "WT1", "KO", "no", "", "", ""]]))


def test_a_copied_trace_and_a_rerun_well_are_warned_about(tmp_path, traces):
    (traces / "a_copy.ab1").write_bytes((traces / "a.ab1").read_bytes())
    sh = samplesheet.load(_sheet(tmp_path / "s.csv", [
        ["traces/a.ab1", "WT1", "WT", "no", "100", "A1", ""],
        ["traces/a_copy.ab1", "WT2", "WT", "no", "100", "A2", ""],
        ["traces/b.ab1", "WT1", "WT", "no", "100", "B1", ""],
        ["traces/c.ab1", "WT1", "WT", "no", "100", "B1", ""]]))
    text = "\n".join(sh.warnings)
    assert "byte-for-byte the same file as line 2" in text
    assert "well B1 is also line 4" in text
    # the two traces of the re-run well keep distinct names
    ids = [cid for s in sh.samples() for _f, cid in s["path"]]
    assert len(ids) == len(set(ids))


def test_excluded_rows_are_kept_out_with_their_reason(tmp_path, traces):
    sh = samplesheet.load(_sheet(tmp_path / "s.csv", [
        ["traces/a.ab1", "WT1", "WT", "no", "", "", "", ""],
        ["traces/b.ab1", "WT1", "WT", "no", "", "", "", "contaminant"]],
        header=("file", "sample", "group", "control", "order", "well", "label", "exclude")))
    assert [cid for _f, cid in sh.samples()[0]["path"]] == ["a"]
    assert sh.excluded()[0]["reason"] == "contaminant"


def test_draft_lists_every_read_file_by_subfolder(tmp_path, traces):
    out = samplesheet.draft(str(tmp_path), str(tmp_path / "draft.csv"))
    rows = list(csv.DictReader(open(out)))
    assert {r["sample"] for r in rows} == {"traces"}
    assert len(rows) == 3


def test_tab_separated_sheets_are_read(tmp_path, traces):
    p = tmp_path / "s.tsv"
    p.write_text("file\tsample\ntraces/a.ab1\tWT1\n")
    assert samplesheet.load(str(p)).samples()[0]["name"] == "WT1"


def test_the_app_loads_a_sheet_and_reports_its_problems(tmp_path, traces):
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient

    from shmtrace.webapp.main import app
    c = TestClient(app)
    good = _sheet(tmp_path / "s.csv", [["traces/a.ab1", "WT1", "WT", "no", "100", "A1", "w1-1"]])
    r = c.post("/api/samplesheet", json={"path": str(good)})
    assert r.status_code == 200
    assert r.json()["samples"][0]["path"][0][1] == "100-A1"
    assert r.json()["sheet"]["provenance"]["100-A1"]["label"] == "w1-1"
    bad = _sheet(tmp_path / "b.csv", [["traces/none.ab1", "WT1", "WT", "no", "", "", ""]])
    r = c.post("/api/samplesheet", json={"path": str(bad)})
    assert r.status_code == 400 and "file not found" in r.json()["detail"]["problems"][0]
