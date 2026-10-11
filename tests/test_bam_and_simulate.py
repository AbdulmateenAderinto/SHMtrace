"""BAM input, spike-in validation helpers and SHMTool-style per-file counting."""
import random

import pytest

from nanoshm.compare import Unique
from nanoshm.config import Params
from nanoshm.io_utils import read_fastq
from nanoshm.simulate import choose_sites, evaluate, spike

pysam = pytest.importorskip("pysam")


def test_bam_reads_come_back_in_sequenced_orientation(tmp_path):
    path = str(tmp_path / "x.bam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "ref", "LN": 100}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for name, seq, rev in [("fwd", "ACGTTT", False), ("rev", "AAACGT", True),
                               ("rev", "AAAC", True)]:            # split record, shorter
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence = name, seq
            a.flag = 16 if rev else 0
            a.reference_id, a.reference_start = 0, 0
            a.cigar = [(0, len(seq))]
            a.query_qualities = pysam.qualitystring_to_array("I" * (len(seq) - 1) + "#")
            out.write(a)
    reads = {n: (s, q) for n, s, q in read_fastq(path)}
    assert reads["fwd"] == ("ACGTTT", "IIIII#")
    # stored reverse-complemented with qualities reversed: undone on read
    assert reads["rev"] == ("ACGTTT", "#IIIII")
    assert len(reads) == 2


def test_choose_sites_prefers_aid_hotspots():
    ref = "AGCTAGCTTTTTTTTTTTTTTTTTTAGCTAGCTTTTTTTTTTTTTTTTTT" * 4
    rng = random.Random(0)
    sites = choose_sites(ref, 0, len(ref) - 1, 20, rng)
    assert len({s[0] for s in sites}) == 20
    assert {s[2] for s in sites} <= {1, 2, 3, 5, 10, 20}


def test_spike_writes_planted_bases_on_both_strands(tmp_path):
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    (tmp_path / "r.fasta").write_text(f">ref\n{ref}\n")
    comp = str.maketrans("ACGT", "TGCA")
    reads = []
    for i in range(40):
        seq = ref if i % 2 else ref.translate(comp)[::-1]
        reads.append((f"r{i}", seq, "5" * len(seq)))
    p = Params(minimap2_preset="map-ont")
    alt = "A" if ref[300] != "A" else "C"
    spiked, planted = spike(reads, str(tmp_path / "r.fasta"), [(300, alt, 10)], p, rng)
    assert planted[(300, alt)] == 10
    changed = sum(1 for (_n, a, _q), (_m, b, _r) in zip(reads, spiked, strict=False) if a != b)
    assert changed == 10


def test_evaluate_scores_recall_and_false_calls():
    planted = {(10, "T"): 1, (20, "A"): 5, (30, "G"): 0}
    ev = evaluate({(20, "A"), (40, "C")}, planted)
    assert ev["true_calls"] == 1 and ev["false_calls"] == 1
    assert ev["recall_by_reads"] == {1: (0, 1), 5: (1, 1)}


def test_unique_mutations_are_counted_once_per_file_not_pooled():
    u = Unique()
    list.append(u, {(1, "C", "T"), (2, "G", "A")})
    list.append(u, {(1, "C", "T")})
    assert len(u) == 3                      # SHMTool: once per file
    assert sorted(u) == [(1, "C", "T"), (1, "C", "T"), (2, "G", "A")]
