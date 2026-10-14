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


def test_the_platform_is_inferred_from_the_files(tmp_path, traces):
    # a sheet of chromatograms is Sanger, a sheet of reads is nanopore; `run` once
    # defaulted to Sanger, which tested read libraries with the wrong error model
    sanger = samplesheet.load(_sheet(tmp_path / "s.csv", [["traces/a.ab1", "WT1", "WT", "no", "", "", ""]]))
    assert sanger.platform() == "sanger"
    assert sanger.project("ref.fasta", "jh4_intron")["platform"] == "sanger"
    reads = tmp_path / "reads"
    reads.mkdir()
    for n in ("WT1.fastq", "KO1.fastq.gz"):
        (reads / n).write_text("@r\nACGT\n+\nIIII\n")
    ont = samplesheet.load(_sheet(tmp_path / "o.csv", [["reads/WT1.fastq", "WT1", "WT", "no", "", "", ""],
                                                       ["reads/KO1.fastq.gz", "KO1", "KO", "yes", "", "", ""]]))
    assert ont.platform() == "ont"
    assert ont.project("ref.fasta", "jh4_intron")["platform"] == "ont"
    # FASTA clone sequences follow the other files; alone they are Sanger, as before
    (reads / "c1.fasta").write_text(">c1\nACGT\n")
    (traces / "c2.fasta").write_text(">c2\nACGT\n")
    fa = samplesheet.load(_sheet(tmp_path / "f.csv", [["traces/c2.fasta", "WT1", "WT", "no", "", "", ""]]))
    assert fa.platform() == "sanger"
    with_traces = samplesheet.load(_sheet(tmp_path / "t.csv", [["traces/a.ab1", "WT1", "WT", "no", "", "", ""],
                                                               ["traces/c2.fasta", "WT1", "WT", "no", "", "", ""]]))
    assert with_traces.platform() == "sanger"
    # an explicit platform wins
    assert ont.project("ref.fasta", "jh4_intron", "illumina")["platform"] == "illumina"


def test_a_sheet_mixing_traces_and_reads_is_refused(tmp_path, traces):
    (tmp_path / "WT2.fastq").write_text("@r\nACGT\n+\nIIII\n")
    mixed = samplesheet.load(_sheet(tmp_path / "m.csv", [["traces/a.ab1", "WT1", "WT", "no", "", "", ""],
                                                         ["WT2.fastq", "WT2", "WT", "no", "", "", ""]]))
    with pytest.raises(samplesheet.SheetError, match="mixes"):
        mixed.platform()
    with pytest.raises(SystemExit, match="mixes"):
        samplesheet.main(["run", str(tmp_path / "m.csv"), "--reference", "r.fasta", "--region", "jh4_intron",
                          "-o", str(tmp_path / "out")])


def test_the_control_column_names_the_kind_of_control(tmp_path):
    from shmtrace import samplesheet
    for n in ("a", "b", "c", "d"):
        (tmp_path / f"{n}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "s.csv").write_text("file,sample,group,control\n"
                                    "a.fastq,A,wt,no\nb.fastq,B,ko,Empty vector\n"
                                    "c.fastq,C,d0,baseline\nd.fastq,D,ko,yes\n")
    ss = {s["name"]: s for s in samplesheet.load(str(tmp_path / "s.csv")).samples()}
    assert not ss["A"]["control"] and "control_kind" not in ss["A"]
    assert ss["B"]["control"] and ss["B"]["control_kind"] == "empty_vector"
    assert ss["C"]["control"] and ss["C"]["control_kind"] == "baseline"
    assert ss["D"]["control"] and "control_kind" not in ss["D"]
    (tmp_path / "bad.csv").write_text("file,sample,group,control\na.fastq,A,wt,maybe\n")
    import pytest
    with pytest.raises(samplesheet.SheetError, match="kind of control"):
        samplesheet.load(str(tmp_path / "bad.csv"))
