"""Indel calling against the controls' own indel rate (indels.py)."""
import random
from dataclasses import replace

from nanoshm import indels, project
from nanoshm.align import AlignedRead
from nanoshm.config import Params

REF = "GATCCATGCAAATTTGGGCCCATATCGCGATCGATCGTAGCTAGCATGCAGTCAGTCAGGATCCTAGC" * 4
P = replace(Params(), end_trim_bp=0, germline_freq_threshold=0.5, indel_min_reads=3)


def test_deletion_in_a_homopolymer_is_left_aligned():
    # deleting any one of the three As in GCAAAT is the same event
    pos = REF.index("GCAAAT") + 2
    assert {indels.normalise(REF, pos + i, "D", 1) for i in range(3)} == {(pos, "del", "A")}


def test_insertion_is_rotated_to_its_leftmost_form():
    pos = REF.index("GCAAAT") + 5                 # before the T
    assert indels.normalise(REF, pos, "I", "A") == (pos - 3, "ins", "A")


def _reads(n, events_per_read, rng, q=40):
    out = []
    for i in range(n):
        ev = events_per_read(i, rng)
        out.append(AlignedRead(f"r{i}", "ref", 0, len(REF), [], "", 1 if i % 2 else -1,
                               [(pos, op, x, q) for pos, op, x in ev]))
    return out


def test_depth_counts_reads_spanning_the_event():
    reads = [AlignedRead("a", "ref", 0, 50, [], ""), AlignedRead("b", "ref", 30, 100, [], ""),
             AlignedRead("c", "ref", 0, 100, [], "")]
    pile = indels.Pileup(reads, REF, 0, len(REF) - 1)
    assert pile.depth(40, "del", "AAAA") == 3       # needs bases 39 and 44
    assert pile.depth(49, "del", "A") == 2          # read a ends at 50
    assert pile.depth(30, "ins", "A") == 2          # read b starts at 30


def _noise(i, rng):
    # a recurrent 1-bp deletion at one site in 5% of reads (the systematic
    # error a control measures), plus sparse scattered noise
    out = []
    if rng.random() < 0.05:
        out.append((60, "D", 1))
    if rng.random() < 0.05:
        out.append((rng.randint(20, 200), "D", 1))
    return out


def test_controls_call_nothing_and_a_clonal_deletion_is_called():
    rng = random.Random(7)
    lo, hi = 0, len(REF) - 1
    ctrls = [indels.Pileup(_reads(400, _noise, rng), REF, lo, hi) for _ in range(3)]
    null = indels.IndelNull(ctrls[1:], REF)
    null.qscale, _ = indels.calibrate(ctrls, REF, P)
    assert indels.call(ctrls[0], null, P)[0] == []

    def sample(i, r):
        return _noise(i, r) + ([(120, "D", 12)] if i < 10 else [])
    called, stats = indels.call(indels.Pileup(_reads(400, sample, rng), REF, lo, hi), null, P)
    assert [(v.pos, v.kind, v.length, v.supporting_reads) for v in called] == [(121, "del", 12, 10)]
    assert stats["indel_model"] == "control_rate"


def test_too_few_reads_are_not_called_and_low_quality_is_ignored():
    rng = random.Random(8)
    lo, hi = 0, len(REF) - 1
    null = indels.IndelNull([indels.Pileup(_reads(300, _noise, rng), REF, lo, hi)], REF)
    two = indels.Pileup(_reads(300, lambda i, r: [(150, "D", 20)] if i < 2 else [], rng),
                        REF, lo, hi)
    assert indels.call(two, null, P)[0] == []
    lowq = indels.Pileup(_reads(300, lambda i, r: [(150, "D", 20)] if i < 10 else [], rng, q=8),
                         REF, lo, hi, min_quality=30)
    assert lowq.events == {}


def test_sanger_calls_every_indel_without_a_test():
    rng = random.Random(9)
    clones = indels.Pileup(_reads(6, lambda i, r: [(90, "I", "G")] if i == 0 else [], rng),
                           REF, 0, len(REF) - 1)
    called, stats = indels.call(clones, None, P, sanger=True)
    assert len(called) == 1 and stats["indel_model"] == "sanger_threshold"


def test_project_reports_indels(wt_fastq, ko_fastq, reference, contaminant):
    proj = {"name": "t", "region": "jh4_intron", "reference": str(reference),
            "contaminants": [str(contaminant)], "overrides": {"exclude_indels": False},
            "samples": [{"name": "WT", "path": str(wt_fastq), "group": "WT"},
                        {"name": "KO", "path": str(ko_fastq), "group": "KO", "control": True}]}
    res = project.run_project(proj)
    wt = next(r for r in res["samples"] if r["sample"] == "WT")
    assert "indels_called" in wt and isinstance(res["indels"], list)
    assert any("Indels need a control" in n for n in res["notes"])
