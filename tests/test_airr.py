from shmtrace import airr
from shmtrace.align import AlignedRead
from shmtrace.variants import Variant


def test_denoised_alignment_keeps_only_called_mutations():
    ref = "ACGTACGTAC"
    # read: a called mutation at position 3 (G>T), an uncalled mismatch at 6 (C>G), a deletion at 8
    pairs = [(0, "A", 30), (1, "C", 30), (2, "T", 30), (3, "T", 30), (4, "A", 30), (5, "G", 30),
             (6, "G", 30), (8, "A", 30), (9, "C", 30)]
    r = AlignedRead("r1", "ref", 0, 10, pairs, "ACTTAGGAC", -1)
    v = Variant.__new__(Variant)
    v.pos, v.alt = 3, "T"
    row = airr.read_rows("S1", "WT", [r], [v], ref, {"r1": 0})[0]
    assert row["germline_alignment"] == "ACGTACGTAC"
    assert row["sequence_alignment_raw"] == "ACTTAGG-AC"      # every read base, deletion as '-'
    assert row["sequence_alignment"] == "ACTTACG-AC"          # only the called G>T kept; deletion as '-'
    assert len(row["sequence_alignment"]) == len(row["germline_alignment"]) == len(row["sequence_alignment_raw"])
    assert row["shmtrace_mutations"] == 1 and row["shmtrace_positions"] == "3G>T"
    assert row["rev_comp"] == "T" and row["clone_id"] == "S1_0" and row["repertoire_id"] == "S1"
    tsv = airr.write_tsv([row]).splitlines()
    assert tsv[0].split("\t")[:14] == airr.REQUIRED
