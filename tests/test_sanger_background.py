"""The AID-null clones as a Sanger project's null: a measured error rate, not an assumed one.

A Sanger clone is one molecule read once, so a call carries no test. What the controls can
say is how often a base of a given quality is wrong, how many of a library's calls that
accounts for, and which positions the chemistry misreads.
"""
import pytest

pytest.importorskip("mappy", reason="alignment backend not installed")

from shmtrace import project, samplesheet  # noqa: E402
from shmtrace.io_utils import read_fasta  # noqa: E402

ERR_POS = 300          # 1-based: where two control clones mismatch
MUT_POS = [250, 400]   # real mutations, in the samples only
SWAP = {"A": "G", "G": "A", "C": "T", "T": "C"}


def _clone(path, ref, muts):
    seq = list(ref)
    for pos in muts:
        seq[pos - 1] = SWAP[seq[pos - 1]]
    s = "".join(seq)
    path.write_text(f"@{path.stem}\n{s}\n+\n{'I' * len(s)}\n")


@pytest.fixture(scope="module")
def res(tmp_path_factory, reference):
    d = tmp_path_factory.mktemp("sanger_bg")
    ref = read_fasta(str(reference))[0][1].upper()
    rows = ["file,sample,group,control,clone_id"]
    for i in range(6):                                   # controls: two share one bad position
        f = d / f"ko{i}.fastq"
        _clone(f, ref, [ERR_POS] if i < 2 else [])
        rows.append(f"{f},KO,ko,yes,ko{i}")
    for i in range(4):                                   # samples: real mutations, one at the bad one
        f = d / f"wt{i}.fastq"
        _clone(f, ref, ([MUT_POS[i % 2]] if i else [ERR_POS]))
        rows.append(f"{f},WT,wt,,wt{i}")
    sheet = d / "sheet.csv"
    sheet.write_text("\n".join(rows) + "\n")
    sh = samplesheet.load(str(sheet))
    return project.run_project(sh.project(str(reference), "jh4_intron", platform="sanger"))


def test_the_background_is_measured_in_the_controls(res):
    bg = res["sanger_background"]
    assert bg["clones"] == 6
    assert bg["mismatches"] == 2 and bg["bases"] > 2000
    assert bg["rate_per_1e5"] == pytest.approx(1e5 * 2 / bg["bases"], rel=1e-6)
    assert bg["min_quality"] == 20


def test_a_position_the_controls_misread_is_reported(res):
    assert str(ERR_POS) in res["sanger_background"]["recurrent_positions"]
    n_ctl, depth = res["sanger_background"]["recurrent_positions"][str(ERR_POS)]
    assert n_ctl == 2 and depth == 6


def test_every_library_is_told_how_many_calls_the_background_explains(res):
    for s in res["samples"]:
        assert s["expected_false_calls"] >= 0
        if s["unique_mutations"]:
            assert s["expected_false_pct"] == pytest.approx(          # both are rounded for reporting
                100 * s["expected_false_calls"] / s["unique_mutations"], abs=0.3)
    wt = next(s for s in res["samples"] if s["sample"] == "WT")
    assert 0 < wt["expected_false_calls"] < wt["unique_mutations"], "real mutations should dominate"


def test_a_call_at_a_misread_position_is_flagged_but_kept(res):
    flags = [f for f in res["qc"]["flags"] if f["check"] == "control_error_site"]
    assert flags, "the sample call at the controls' bad position should be flagged"
    assert all(f["severity"] == "info" for f in flags), "a note about a position, not a defect"
    assert any(v["pos"] == ERR_POS and v["sample"] == "WT" for v in res["variants"]), \
        "flagging must not remove the call"
