"""Reads hold their aligned columns packed, and everything that reads them sees the same columns."""
import pickle
import random

from shmtrace import browse
from shmtrace.align import AlignedRead, align_reads, build_index
from shmtrace.config import Params
from shmtrace.project import full_counts

REF = "ACGT" * 25


def _read(name, alts=None, strand=1):
    alts = alts or {}
    pairs = [(i, alts.get(i, REF[i]), 30 + i % 7) for i in range(len(REF))]
    return AlignedRead(name, "ref", 0, len(REF), pairs, REF, strand, [], 0, len(REF))


def test_a_compacted_read_is_the_same_read():
    plain, packed = _read("r", {9: "T"}), _read("r", {9: "T"}).compact()
    assert "_packed" in packed.__dict__ and "_pairs" not in packed.__dict__
    assert packed == plain and packed.pairs == plain.pairs and repr(packed) == repr(plain)


def test_a_read_built_from_a_list_keeps_that_list_so_it_can_be_edited():
    r = _read("r")
    r.pairs[9] = (9, "T", 5)
    assert r.pairs[9] == (9, "T", 5)
    r.pairs = [(0, "A", 40)]
    assert r.compact().pairs == [(0, "A", 40)]


def test_a_read_crosses_processes_packed_and_equal_with_what_was_set_on_it():
    for r in (_read("a"), _read("b").compact()):
        r.members = ["x", "y"]
        back = pickle.loads(pickle.dumps(r))
        assert back == r and back.members == ["x", "y"] and "_packed" in back.__dict__


def test_columns_that_cannot_be_packed_stay_as_they_are():
    odd = AlignedRead("r", "ref", 0, 2, [(0, "", 30), (1, "A", 30)], "A").compact()
    assert "_pairs" in odd.__dict__ and odd.pairs == [(0, "", 30), (1, "A", 30)]
    assert pickle.loads(pickle.dumps(odd)) == odd


def test_counts_and_the_read_view_do_not_depend_on_the_storage():
    class V:
        def __init__(self, pos, alt):
            self.pos, self.alt = pos, alt
    plain = [_read(f"r{i}", {9: "T"} if i < 3 else {}) for i in range(10)]
    packed = [_read(f"r{i}", {9: "T"} if i < 3 else {}).compact() for i in range(10)]
    calls = [V(10, "T"), V(20, "A")]
    assert full_counts(plain, calls, Params()) == full_counts(packed, calls, Params())
    assert browse.read_rows(plain, REF) == browse.read_rows(packed, REF)


def test_reads_come_out_of_the_aligner_packed(tmp_path):
    rng = random.Random(0)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    (tmp_path / "ref.fa").write_text(f">ref\n{ref}\n")
    aligner, targets, contams = build_index([str(tmp_path / "ref.fa")], [], preset="map-ont")
    read = ref[50:550]
    out = list(align_reads([("r1", read, "I" * len(read))], aligner, targets, contams,
                           Params(min_align_len=100, max_total_len=2000)))
    assert len(out) == 1 and "_packed" in out[0].__dict__
    assert [b for _r, b, _q in out[0].pairs] == list(read)
