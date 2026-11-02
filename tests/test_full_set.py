"""The reads downsampling leaves out: counted per call, and kept for Check calls to show."""
from shmtrace import browse
from shmtrace.align import AlignedRead
from shmtrace.config import Params
from shmtrace.project import full_counts, orientation_notes, phasing_notes

REF = "ACGT" * 25


def _read(name, alts=None, strand=1):
    alts = alts or {}
    pairs = [(i, alts.get(i, REF[i]), 30) for i in range(len(REF))]
    return AlignedRead(name, "ref", 0, len(REF), pairs, REF, strand, [], 0, len(REF))


class _V:
    def __init__(self, pos, alt):
        self.pos, self.alt = pos, alt


def test_full_counts_counts_carriers_and_depth_over_the_given_reads():
    reads = [_read(f"r{i}", {9: "T"} if i < 3 else {}) for i in range(10)]
    out = full_counts(reads, [_V(10, "T"), _V(20, "A")], Params(min_base_quality=0))
    assert out == {"10T": [3, 10], "20A": [0, 10]}


def test_quality_gate_applies_as_it_does_to_the_pileup():
    r = _read("low", {9: "T"})
    r.pairs[9] = (9, "T", 5)
    assert full_counts([r], [_V(10, "T")], Params(min_base_quality=10)) == {"10T": [0, 0]}


def test_undrawn_reads_join_the_alignment_after_the_analysed_ones():
    data = {"reads": [{"name": "a", "fate": "kept", "start": 1, "end": 100, "strand": 1, "diffs": [], "indels": []},
                      {"name": "b", "fate": "not sampled", "start": 1, "end": 0, "strand": 1, "diffs": [], "indels": []}],
            "counts": {}}
    browse.add_undrawn(data, [_read("b", {9: "T"}), _read("c")], REF)
    assert [(r["name"], r["fate"]) for r in data["reads"]] == [("a", "kept"), ("b", "not sampled"), ("c", "not sampled")]
    assert data["reads"][1]["diffs"] == [[10, "T", 30]]
    assert data["counts"] == {"kept": 1, "not sampled": 2}


def test_a_library_in_one_orientation_is_noted():
    one = {"X": [_read(f"r{i}") for i in range(60)],
           "Y": [_read(f"s{i}", strand=1 if i % 2 else -1) for i in range(60)]}
    notes = orientation_notes(one)
    assert len(notes) == 1 and "X" in notes[0] and "Y" not in notes[0]
    assert orientation_notes({"Z": [_read("t")]}) == []          # too few reads to say


def test_phasing_notes_name_a_capped_block_and_a_polyclonal_drop():
    stats = {"W1": {"aligned_before_phasing": 100, "blocks": [
        {"n_sites": 4, "reads": 41, "action": "kept (cap)", "why": "", "junctions": None},
        {"n_sites": 6, "reads": 12, "action": "drop", "why": "carried by reads of 9 different V(D)J junctions",
         "junctions": {"junctions": 9}}]}}
    notes = phasing_notes(stats, {"W1": []})
    assert len(notes) == 2
    assert "41 of 100" in notes[1] or "41 of 100" in notes[0]
    assert any("9 different V(D)J junctions" in n for n in notes)
