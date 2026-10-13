"""LoFreq (Wilm et al. 2012) on exactly the reads SHMtrace analysed, for the benchmark.

    from lofreq_runner import lofreq_calls
    calls = lofreq_calls(reads, reference, lo, hi, preset="map-ont")

`reads` is a list of (name, seq, qual). They are aligned with minimap2 (the same preset SHMtrace
uses), sorted and indexed with samtools, and called with `lofreq call` at its defaults
(Bonferroni-corrected Poisson-binomial test, significance 0.01, default filters), restricted to the
analysed window. Substitutions only. Sites at >= 50% allele frequency are removed, as SHMtrace does,
so every caller is scored on the same somatic definition. Returns {(pos1, ref, alt): af}.
"""
import os
import subprocess
import tempfile

BIN = os.path.expanduser("~/miniforge3/envs/lofreq/bin")
# LoFreq settings benchmarked:
#   lofreq             its defaults (BAQ on, Bonferroni, default filters)
#   lofreq_nobaq       -B: base-alignment quality off. BAQ was built for short reads and on nanopore
#                      reads it suppresses almost every call, so this is LoFreq's realistic long-read use
#   lofreq_permissive  -B, uncorrected p <= 0.05, default filters off: the analogue of Geneious's rule
LOFREQ = {"lofreq": (), "lofreq_nobaq": ("-B",),
          "lofreq_permissive": ("-B", "-a", "0.05", "-b", "1", "--no-default-filter")}


def _which_minimap2():
    for c in (os.path.join(BIN, "minimap2"), "minimap2"):
        try:
            subprocess.run([c, "--version"], capture_output=True, check=True)
            return c
        except (OSError, subprocess.CalledProcessError):
            continue
    raise RuntimeError("minimap2 not found")


def lofreq_calls(reads, reference, lo, hi, preset="map-ont", germline=0.5, extra=()):
    mm2 = _which_minimap2()
    with tempfile.TemporaryDirectory() as d:
        fq, bam = os.path.join(d, "r.fastq"), os.path.join(d, "r.bam")
        with open(fq, "w") as fh:
            for n, s, q in reads:
                fh.write(f"@{n}\n{s}\n+\n{q}\n")
        ref = os.path.join(d, "ref.fa")
        with open(reference) as src, open(ref, "w") as dst:
            dst.write(src.read())
        subprocess.run(f"{mm2} -a -x {preset} --secondary=no {ref} {fq} 2>/dev/null | "
                       f"{BIN}/samtools sort -o {bam} - && {BIN}/samtools index {bam} && "
                       f"{BIN}/samtools faidx {ref}", shell=True, check=True)
        contig = open(ref).readline()[1:].split()[0]
        vcf = os.path.join(d, "out.vcf")
        subprocess.run([f"{BIN}/lofreq", "call", "-f", ref, "-r", f"{contig}:{lo}-{hi}",
                        "-o", vcf, *extra, bam], check=True, capture_output=True)
        out = {}
        for line in open(vcf):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            pos, r, a, info = int(f[1]), f[3], f[4], f[7]
            if len(r) != 1 or len(a) != 1:
                continue
            af = float(next(x for x in info.split(";") if x.startswith("AF="))[3:])
            if af < germline and lo <= pos <= hi:
                out[(pos, r, a)] = af
        return out
