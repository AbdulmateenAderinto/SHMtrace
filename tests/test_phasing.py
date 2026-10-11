"""Minor-haplotype detection: it must fire on co-occurrence, not on frequency."""
from dataclasses import dataclass

from nanoshm.config import Params
from nanoshm.phasing import candidate_sites, find_blocks, partition

REF = "ACGT" * 50          # 200 bp
DROP = Params(phase_mode="drop")


@dataclass
class FakeRead:
    pairs: list


def read(alts: dict, length=200, q=30):
    """A read matching REF except at the given {pos: base} overrides."""
    return FakeRead([(i, alts.get(i, REF[i]), q) for i in range(length)])


def test_haplotype_block_is_removed():
    block = {10: "T", 20: "T", 30: "T", 40: "T", 50: "T"}
    reads = [read({}) for _ in range(90)] + [read(block) for _ in range(10)]
    kept, stats = partition(reads, REF, DROP, 0, 199)
    assert stats["phasing_applied"]
    assert stats["reads_dropped_minor"] == 10
    assert len(kept) == 90


def test_low_abundance_block_is_still_removed():
    """The JH4 block sat at ~2% of reads in one mouse and ~17% in another, so
    detection must not depend on the block being frequent."""
    block = {10: "T", 20: "T", 30: "T", 40: "T", 50: "T"}
    reads = [read({}) for _ in range(980)] + [read(block) for _ in range(20)]
    _kept, stats = partition(reads, REF, DROP, 0, 199)
    assert stats["phasing_applied"], stats["reason"]
    assert stats["reads_dropped_minor"] == 20


def test_a_handful_of_reads_is_not_enough_to_seed():
    """Below a few reads a 'block' is indistinguishable from noise; do nothing."""
    block = {10: "T", 20: "T", 30: "T", 40: "T", 50: "T"}
    reads = [read({}) for _ in range(196)] + [read(block) for _ in range(4)]
    kept, stats = partition(reads, REF, DROP, 0, 199)
    assert not stats["phasing_applied"]
    assert len(kept) == 200


def test_scattered_mutations_are_not_phased_away():
    """Real SHM does not co-occur. Each read gets its own mutation."""
    reads = [read({10 + i: "T"}) for i in range(60)] + [read({}) for _ in range(40)]
    kept, stats = partition(reads, REF, DROP, 0, 199)
    assert not stats["phasing_applied"], stats["reason"]
    assert len(kept) == 100


def test_clean_reads_are_untouched():
    reads = [read({}) for _ in range(50)]
    kept, stats = partition(reads, REF, DROP, 0, 199)
    assert not stats["phasing_applied"]
    assert len(kept) == 50


def test_refuses_to_drop_more_than_the_cap():
    """A 'minor' population that is most of the data is a reference problem."""
    block = {10: "T", 20: "T", 30: "T", 40: "T", 50: "T"}
    # 48%: under the germline rule (so still a "variant") but over phase_max_drop
    reads = [read({}) for _ in range(520)] + [read(block) for _ in range(480)]
    kept, stats = partition(reads, REF, DROP, 0, 199)
    assert not stats["phasing_applied"]
    assert "refusing to phase" in stats["reason"], stats["reason"]
    assert len(kept) == 1000


def test_phasing_can_be_disabled():
    block = {10: "T", 20: "T", 30: "T", 40: "T", 50: "T"}
    reads = [read({}) for _ in range(90)] + [read(block) for _ in range(10)]
    kept, stats = partition(reads, REF, Params(phase_reads=False), 0, 199)
    assert stats["reason"] == "disabled"
    assert len(kept) == 100


def test_germline_positions_are_not_candidates():
    """An alt at >=50% is the reference being wrong, handled elsewhere."""
    reads = [read({10: "T"}) for _ in range(80)] + [read({}) for _ in range(20)]
    sites = candidate_sites(reads, REF, DROP, 0, 199)
    assert (10, "T") not in sites


BLOCK = {10: "T", 20: "T", 30: "T", 40: "T", 50: "T"}
CLONE = {110: "T", 120: "T", 130: "T", 140: "T"}


def test_auto_mode_flags_an_unknown_block_instead_of_dropping_it():
    """An expanded clone shares its mutations exactly like a haplotype does.
    Without outside evidence the block must be reported, not removed."""
    reads = [read({}) for _ in range(90)] + [read(BLOCK) for _ in range(10)]
    kept, stats = partition(reads, REF, Params(), 0, 199)
    assert not stats["phasing_applied"]
    assert len(kept) == 100
    assert stats["blocks"][0]["action"] == "flag"
    assert stats["blocks"][0]["n_sites"] == 5


def test_auto_mode_drops_a_block_seen_in_controls():
    reads = [read({}) for _ in range(90)] + [read(BLOCK) for _ in range(10)]
    known = [[(pos, alt) for pos, alt in BLOCK.items()]]
    kept, stats = partition(reads, REF, Params(), 0, 199, known_blocks=known)
    assert stats["phasing_applied"]
    assert stats["reads_dropped_minor"] == 10


def test_auto_mode_keeps_a_clone_while_dropping_a_known_haplotype():
    reads = ([read({}) for _ in range(80)] + [read(BLOCK) for _ in range(10)]
             + [read(CLONE) for _ in range(10)])
    known = [[(pos, alt) for pos, alt in BLOCK.items()]]
    kept, stats = partition(reads, REF, Params(), 0, 199, known_blocks=known)
    assert stats["reads_dropped_minor"] == 10
    actions = sorted(b["action"] for b in stats["blocks"])
    assert actions == ["drop", "flag"]
    assert len(kept) == 90


def test_two_independent_blocks_are_both_found():
    reads = ([read({}) for _ in range(80)] + [read(BLOCK) for _ in range(10)]
             + [read(CLONE) for _ in range(10)])
    blocks = find_blocks(reads, REF, Params(), 0, 199)
    assert len(blocks) == 2


def test_flag_mode_never_drops():
    reads = [read({}) for _ in range(90)] + [read(BLOCK) for _ in range(10)]
    known = [[(pos, alt) for pos, alt in BLOCK.items()]]
    kept, stats = partition(reads, REF, Params(phase_mode="flag"), 0, 199, known_blocks=known)
    assert len(kept) == 100
