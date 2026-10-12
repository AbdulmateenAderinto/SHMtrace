"""Phred and expected-error maths, and the minimal FASTA/FASTQ readers."""
import math

import pytest

from shmtrace.io_utils import expected_error_rate, phred, read_fasta, read_fastq


def test_phred_offset_is_33():
    assert phred("!") == [0]
    assert phred("I") == [40]
    assert phred("+5?") == [10, 20, 30]


@pytest.mark.parametrize("q,expected", [(10, 0.1), (20, 0.01), (30, 0.001), (40, 0.0001)])
def test_expected_error_rate_matches_phred_definition(q, expected):
    ch = chr(q + 33)
    assert expected_error_rate(ch * 50) == pytest.approx(expected)


def test_expected_error_rate_is_the_mean_not_the_product():
    """Q20 throughout and Q10 throughout must average, not multiply."""
    qual = chr(20 + 33) * 10 + chr(10 + 33) * 10
    assert expected_error_rate(qual) == pytest.approx((0.01 + 0.1) / 2)


def test_expected_error_rate_empty_is_worst_case():
    assert expected_error_rate("") == 1.0


def test_q20_is_the_default_threshold_boundary():
    """max_ee_rate=0.01 is documented as Q20-equivalent; hold that link."""
    from shmtrace.config import Params

    assert expected_error_rate(chr(20 + 33) * 100) == pytest.approx(Params().max_ee_rate)


def test_read_fasta_handles_wrapped_lines(tmp_path):
    f = tmp_path / "r.fasta"
    f.write_text(">chr1 some description\nacgt\nACGT\n>chr2\nTTTT\n")
    recs = read_fasta(str(f))
    assert recs == [("chr1", "ACGTACGT"), ("chr2", "TTTT")]


def test_read_fastq_yields_name_seq_qual(tmp_path):
    f = tmp_path / "r.fastq"
    f.write_text("@read1 extra\nACGT\n+\nIIII\n@read2\nTTTT\n+\n!!!!\n")
    recs = list(read_fastq(str(f)))
    assert [r[0] for r in recs] == ["read1", "read2"]
    assert recs[0][1] == "ACGT"
    assert recs[1][2] == "!!!!"


def test_bundled_example_reads_are_readable(wt_fastq):
    reads = list(read_fastq(str(wt_fastq)))
    assert len(reads) > 0
    for _name, seq, qual in reads[:20]:
        assert len(seq) == len(qual), "FASTQ seq/qual length mismatch"
        assert not math.isnan(expected_error_rate(qual))


def test_genbank_records_are_reads(tmp_path):
    from shmtrace.io_utils import read_fastq, sniff_format
    gb = tmp_path / "clone.gb"
    gb.write_text("LOCUS       w2-1-user_added     24 bp ds-DNA     linear\n"
                  "DEFINITION  .\nFEATURES             Location/Qualifiers\nORIGIN\n"
                  "        1 acgtacgtac gtmkacgtac gtac\n//\n")
    assert sniff_format(str(gb)) == "genbank"
    (name, seq, qual), = list(read_fastq(str(gb)))
    assert name == "w2-1-user_added" and seq == "ACGTACGTACGTMKACGTACGTAC"
    assert len(qual) == len(seq)
    # a GenBank file without the extension is sniffed by its LOCUS line
    other = tmp_path / "clone_noext"
    other.write_text(gb.read_text())
    assert sniff_format(str(other)) == "genbank"
