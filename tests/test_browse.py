"""Read alignments kept for inspection (browse.py) and the fates that label them."""
import pytest

from shmtrace import browse, project


@pytest.fixture(scope="module")
def run(examples):
    return project.run_project({
        "name": "t", "region": "jh4_intron",
        "reference": str(examples / "JH4_reference.fasta"),
        "contaminants": [str(examples / "JH4_contaminant.fasta")],
        "browse": True,
        "samples": [{"name": "WT", "path": str(examples / "WT3_JH4_sub500.fastq"),
                     "group": "WT"},
                    {"name": "KO", "path": str(examples / "KO3_JH4_sub500.fastq"),
                     "group": "KO", "control": True}]})


def test_every_input_read_is_accounted_for(run):
    a = run["_alignments"]["WT"]
    analysed = next(s["reads"] for s in run["samples"] if s["sample"] == "WT")
    assert a["counts"]["kept"] == analysed
    # nothing is recorded twice, and every row carries a known fate
    names = [r["name"] for r in a["reads"]]
    assert len(names) == len(set(names))
    assert {r["fate"] for r in a["reads"]} <= set(browse.FATE_ORDER)


def test_a_kept_read_carries_its_differences_not_its_sequence(run):
    kept = [r for r in run["_alignments"]["WT"]["reads"] if r["fate"] == "kept"]
    r = next(x for x in kept if x["diffs"])
    pos, base, qual = r["diffs"][0]
    assert r["start"] <= pos <= r["end"]
    assert base in "ACGT" and 0 <= qual <= 93
    assert "sequence" not in r            # only the differences travel


def test_window_clips_to_the_region_and_pages(run):
    a = run["_alignments"]["WT"]
    w = browse.window(a, 200, 300, limit=10)
    assert w["start"] == 200 and len(w["sequence"]) == 101
    assert len(w["reads"]) == 10 and w["total"] >= 10
    assert all(p >= 200 and p <= 300 for c in [w["called"]] for p in [x["pos"] for x in c])
    second = browse.window(a, 200, 300, limit=10, offset=10)
    assert second["reads"] != w["reads"]


def test_window_filters_by_fate(run):
    a = run["_alignments"]["WT"]
    only = browse.window(a, 1, 1230, fates={"kept"}, limit=5000)
    assert {r["fate"] for r in only["reads"]} == {"kept"}
    assert only["total"] == a["counts"]["kept"]


def test_each_stage_records_why_it_dropped_a_read(tmp_path, examples):
    """A read rejected at any stage is labelled with the stage that rejected it."""
    from shmtrace import pipeline
    from shmtrace.io_utils import read_fasta
    from shmtrace.project import build_params
    ref = str(examples / "JH4_reference.fasta")
    seq = read_fasta(ref)[0][1]
    good, bad = seq[40:900], seq[40:900]
    fq = tmp_path / "mixed.fastq"
    fq.write_text(
        f"@clean\n{good}\n+\n{chr(33 + 40) * len(good)}\n"            # passes
        f"@noisy\n{bad}\n+\n{chr(33 + 2) * len(bad)}\n"               # fails the EE gate
        f"@alien\n{'ACGT' * 220}\n+\n{chr(33 + 40) * 880}\n"          # does not align
        f"@stub\n{seq[300:460]}\n+\n{chr(33 + 40) * 160}\n")          # aligns too short
    p = build_params({"region": "jh4_intron"})
    fates: dict = {}
    kept, _stats = pipeline.prepare_reads(str(fq), ref, [], p, fates=fates)
    assert [a.qname for a in kept] == ["clean"]
    assert fates["noisy"]["fate"] == "low quality"
    assert fates["alien"]["fate"] == "unmapped"
    assert fates["stub"]["fate"] == "alignment too short"
    # a read rejected after alignment still knows where it aligned
    assert fates["stub"]["r_st"] < fates["stub"]["r_en"]


def test_dropped_reads_keep_their_place_against_the_reference(run):
    """Depth normalisation leaves the reads it did not sample in the record."""
    rows = [r for r in run["_alignments"]["WT"]["reads"] if r["fate"] != "kept"]
    assert all(r["end"] >= r["start"] or r["end"] == 0 for r in rows)


def test_the_cap_keeps_the_analysed_reads_first():
    class R:
        def __init__(self, n):
            self.qname, self.ctg, self.r_st, self.r_en = n, "r", 0, 10
            self.pairs, self.strand, self.indels = [], 1, []
    import tempfile, pathlib
    fa = pathlib.Path(tempfile.mkdtemp()) / "r.fasta"
    fa.write_text(">r\n" + "A" * 40 + "\n")
    kept = [R(f"k{i}") for i in range(5)]
    fates = {f"d{i}": {"fate": "unmapped"} for i in range(50)}
    out = browse.alignment(kept, fates, str(fa), [], max_reads=5)
    assert {r["fate"] for r in out["reads"]} == {"kept"}
    assert out["truncated"] == 50


def test_public_says_whether_alignments_exist(run):
    assert project.public(run)["has_alignments"] is True
    assert "_alignments" not in project.public(run)
