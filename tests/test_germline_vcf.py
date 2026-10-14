import gzip

from shmtrace import germline_vcf

REF = "ACGTACGTACGT"
HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmouse\n"


def _write(tmp_path, body, gz=False):
    f = tmp_path / ("g.vcf.gz" if gz else "g.vcf")
    text = HEAD + body
    if gz:
        with gzip.open(f, "wt") as fh:
            fh.write(text)
    else:
        f.write_text(text)
    return str(f)


def test_reads_carried_substitutions_on_the_amplicon(tmp_path):
    body = ("amp\t2\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\n"          # het: kept
            "amp\t3\t.\tG\tA,T\t30\t.\t.\tGT:DP\t1/2:40\n"      # multi-allelic: both
            "amp\t4\t.\tT\tC\t30\tPASS\t.\tGT\t0/0\n"           # not carried
            "amp\t5\t.\tA\tG\t30\tLowQual\t.\tGT\t1/1\n"         # filtered
            "amp\t6\t.\tCG\tC\t30\tPASS\t.\tGT\t0/1\n"           # deletion: not a substitution
            "chr1\t7\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n")          # another contig
    sites, st = germline_vcf.read_sites(_write(tmp_path, body), "amp - amplicon", REF)
    assert sites == {(1, "T"), (2, "A"), (2, "T")}
    assert st["not_carried"] == 1 and st["filtered"] == 1 and st["not_snv"] == 1 and st["other_contig"] == 1
    assert germline_vcf.describe("g.vcf", st) is None


def test_warns_when_the_vcf_is_not_in_amplicon_coordinates(tmp_path):
    _s, st = germline_vcf.read_sites(_write(tmp_path, "chr12\t113000\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n"),
                                     "amp", REF)
    assert "same FASTA" in germline_vcf.describe("g.vcf", st)
    _s, st = germline_vcf.read_sites(_write(tmp_path, "amp\t1\t.\tT\tG\t30\tPASS\t.\tGT\t0/1\n", gz=True),
                                     "amp", REF)
    assert st["ref_mismatch"] == 1 and "coordinates" in germline_vcf.describe("g.vcf", st)
