"""Hand edits to Sanger base calls (edits.py), on a clone from the bundled Sanger example."""
from pathlib import Path

import pytest

from shmtrace import edits
from shmtrace.io_utils import read_abif, read_fastq

ROOT = Path(__file__).resolve().parents[1]
CLONE = ROOT / "examples_sanger/traces/berrumar_391861_chromats/7-2.1-7-user_added.ab1"
REF = str(ROOT / "src/shmtrace/references/mouse_jh4_intron.fasta")
REF_SEQ = "".join(open(REF).read().split("\n")[1:])
pytestmark = pytest.mark.skipif(not CLONE.is_file(), reason="no Sanger example")


def test_an_edit_records_what_the_trace_called_and_refuses_no_ops():
    e = edits.record([], "WT_2w1", "391861-7", str(CLONE), REF, REF_SEQ, 104, "g", "shoulder")
    assert e == [{**e[0], "sample": "WT_2w1", "clone": "391861-7", "pos": 104, "reference": "G",
                  "called": "A", "to": "G", "note": "shoulder"}]
    with pytest.raises(ValueError, match="already calls A"):
        edits.record([], "WT_2w1", "391861-7", str(CLONE), REF, REF_SEQ, 104, "A")
    with pytest.raises(ValueError, match="one of A, C, G, T or N"):
        edits.record([], "WT_2w1", "391861-7", str(CLONE), REF, REF_SEQ, 104, "X")
    assert edits.record(e, "WT_2w1", "391861-7", "", REF, REF_SEQ, 104, "") == []


def test_the_edited_copy_changes_one_base_on_the_reads_own_strand(tmp_path):
    e = edits.record([], "WT_2w1", "391861-7", str(CLONE), REF, REF_SEQ, 104, "G")
    proj = {"reference": REF, "samples": [{"name": "WT_2w1", "path": [[str(CLONE), "391861-7"]]},
                                          {"name": "KO", "path": "/elsewhere.fastq"}]}
    new = edits.edited_project(proj, e, tmp_path)
    fq, rid = new["samples"][0]["path"][0]
    assert rid == "391861-7" and fq.endswith("391861-7.edited.fastq")
    assert new["samples"][1] == proj["samples"][1] and new["original_samples"] == proj["samples"]
    (name, seq, qual), = list(read_fastq([[fq, rid]]))
    _n, orig, oq = read_abif(str(CLONE))
    assert name == "391861-7" and qual == oq
    changed = [(a, b) for a, b in zip(orig, seq) if a != b]
    assert changed == [("T", "C")]                 # A->G along the reference, on a reverse clone
    again = edits.edited_project(new, e, tmp_path / "2")      # a re-run of a re-run edits the trace
    assert again["samples"][0]["path"][0][0] != fq and again["original_samples"] == proj["samples"]
    assert "391861-7 104 A→G" in edits.notes(e)[0]
