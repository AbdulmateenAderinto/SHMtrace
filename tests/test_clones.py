"""Clones from the V(D)J junction (clones.py)."""
from types import SimpleNamespace

import pytest

from nanoshm import clones, project
from nanoshm.config import Params

J = 52          # 0-based start of the J segment in these synthetic reads


def read(name, junction, muts=(), r_st=J, length=200):
    """A read whose alignment starts at the J segment, preceded by `junction`."""
    pairs = [(J + i, "A", 40) for i in range(length)]
    for pos in muts:
        pairs[pos - J] = (pos, "G", 40)
    return SimpleNamespace(qname=name, r_st=r_st, r_en=J + length, seq=junction + "A" * length,
                           q_st=len(junction), pairs=pairs, strand=1)


def test_edit_distance_with_a_cap():
    assert clones._edit("ACGTACGT", "ACGTACGT", 3) == 0
    assert clones._edit("ACGTACGT", "ACGAACGT", 3) == 1
    assert clones._edit("ACGTACGT", "ACGTCGT", 3) == 1          # an indel, not a shift
    assert clones._edit("AAAAAAAA", "CCCCCCCC", 3) == 4          # stops past the cap


def test_junction_needs_an_alignment_starting_at_j():
    j = "TTACTACGGTAGTAGCTAGG"
    assert clones.junction(read("r", j), J, 20) == j
    assert clones.junction(read("r", j, r_st=0), J, 20) is None         # unrearranged
    assert clones.junction(read("r", "ACGT"), J, 20) is None             # too little 5'


def test_clustering_tolerates_sequencing_error_but_not_a_different_clone():
    a = "TTACTACGGTAGTAGCTAGG"
    a_err = "TTACTACGGTAGTAGCTAGC"                 # one nanopore error
    b = "GGGAATACTACGGTAGTAGA"                     # a different clone
    out = clones.cluster({"r1": a, "r2": a, "r3": a_err, "r4": b}, 3)
    assert out["r1"] == out["r2"] == out["r3"] != out["r4"]


def test_a_clone_counts_its_shared_mutations_once():
    p = Params(junction_ref_start=J + 1)
    a, b = "TTACTACGGTAGTAGCTAGG", "GGGAATACTACGGTAGTAGA"
    reads = ([read(f"a{i}", a, muts=(100, 150)) for i in range(10)]     # one clone, 10 reads
             + [read("b0", b, muts=(120,))])
    variants = [SimpleNamespace(pos=p0 + 1, alt="G") for p0 in (100, 150, 120)]
    out = clones.call_clones(reads, variants, p, J, J + 199)
    assert out["clones"] == 2
    big = next(r for r in out["rows"] if r["reads"] == 10)
    assert big["mutations"] == 2                    # not 20
    assert out["clonal_mutation_events"] == 3
    # read-weighted would be dominated by the big clone; clone-weighted is not
    assert out["largest_clone_pct"] == pytest.approx(100 * 10 / 11, abs=0.1)
    assert out["sharing"]["within"] > out["sharing"]["between"]


def test_no_junction_setting_means_no_clones():
    assert clones.call_clones([read("r", "A" * 20)], [], Params(), J, J + 100) is None


def test_project_reports_clones(examples):
    res = project.run_project({
        "name": "t", "region": "jh4_intron",
        "reference": str(examples / "JH4_reference.fasta"),
        "contaminants": [str(examples / "JH4_contaminant.fasta")],
        "samples": [{"name": "WT", "path": str(examples / "WT3_JH4_sub500.fastq"), "group": "WT"},
                    {"name": "KO", "path": str(examples / "KO3_JH4_sub500.fastq"),
                     "group": "KO", "control": True}]})
    wt = next(s for s in res["samples"] if s["sample"] == "WT")
    assert wt["clones"] > 5 and wt["clone_weighted_frequency_x1e3"] > 0
    sh = res["clones"]["WT"]["sharing"]
    assert sh["within"] > 5 * sh["between"]        # the junction is a barcode


def test_sample_warnings_can_be_added_after_the_row_exists(examples, monkeypatch):
    """Regression: the coverage warning once referenced the row before it was built."""
    monkeypatch.setattr(project, "coverage_dips", lambda dep, start: "coverage dips (test)")
    res = project.run_project({
        "name": "t", "region": "jh4_intron",
        "reference": str(examples / "JH4_reference.fasta"),
        "samples": [{"name": "WT", "path": str(examples / "WT3_JH4_sub500.fastq"), "group": "WT"}]})
    assert "coverage dips (test)" in res["samples"][0]["warnings"]
