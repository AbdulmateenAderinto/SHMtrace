"""The analysed reads as SAM/BAM (samout.py)."""
import pytest

from shmtrace import samout
from shmtrace.io_utils import read_fastq


def test_reads_are_written_as_sorted_sam_with_soft_clips_and_strands(reference, wt_fastq):
    names = {n for n, _s, _q in list(read_fastq(str(wt_fastq)))[:40]}
    text, n = samout.sam(str(reference), str(wt_fastq), names, "ont", "WT3")
    lines = [l for l in text.splitlines() if not l.startswith("@")]
    assert n == len(lines) and 30 <= n <= 40
    assert text.startswith("@HD\tVN:1.6\tSO:coordinate") and "@RG\tID:WT3\tSM:WT3" in text
    pos = [int(l.split("\t")[3]) for l in lines]
    assert pos == sorted(pos)
    by_name = {n: (s, q) for n, s, q in read_fastq(str(wt_fastq))}
    for l in lines[:10]:
        f = l.split("\t")
        assert f[0] in names and f[1] in ("0", "16")
        # the CIGAR covers the whole read, soft clips included
        import re
        length = sum(int(k) for k, op in re.findall(r"(\d+)([MISX=])", f[5]))
        assert length == len(f[9]) == len(by_name[f[0]][0])


def test_bam_and_index_are_made_where_pysam_is_installed(reference, wt_fastq):
    pytest.importorskip("pysam")
    names = {n for n, _s, _q in list(read_fastq(str(wt_fastq)))[:10]}
    text, _n = samout.sam(str(reference), str(wt_fastq), names, "ont", "WT3")
    bam, bai = samout.bam_bytes(text)
    assert bam[:2] == b"\x1f\x8b" and bai[:3] == b"BAI"
