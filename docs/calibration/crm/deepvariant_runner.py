"""DeepVariant (Poplin et al. 2018, germline) and DeepSomatic (Park et al. 2025, tumor-normal somatic) on the reads
SHMtrace analysed.

    from deepvariant_runner import deepvariant_calls, deepsomatic_calls
    calls = deepsomatic_calls(case_reads, normal_reads, reference, lo, hi, setting="deepsomatic")

Reads are aligned exactly as for ClairS (clairs_runner: minimap2 map-ont, samtools) and called with Google's
official Docker images restricted to the analysed window. As for every caller in the benchmark: substitutions
only, FILTER PASS, and sites at >= 50% allele fraction removed. Returns {(pos1, ref, alt): af}.

Settings:
  deepvariant              ONT_R104 model at its defaults (candidate SNVs need >= 12% of reads)
  deepvariant_permissive   candidates from 1% of reads and 2 supporting reads, the lowest the caller is asked to
                           consider
  deepsomatic              ONT model, case against matched normal, at its defaults
  deepsomatic_permissive   candidates from 0.5% of reads and 2 supporting reads
DeepVariant, like Clair3, looks for genotypes; DeepSomatic, like ClairS, compares a case with a normal.
"""
import gzip
import os
import shutil
import subprocess
import tempfile

from clairs_runner import CONTEXT, _bam, _mm2, _ref

DV_IMAGE = "google/deepvariant:1.9.0"
DS_IMAGE = "google/deepsomatic:1.9.0"
DEEPVARIANT = {"deepvariant": (),
               "deepvariant_permissive": ("--make_examples_extra_args=vsc_min_fraction_snps=0.01,vsc_min_count_snps=2",)}
DEEPSOMATIC = {"deepsomatic": (),
               "deepsomatic_permissive": ("--make_examples_extra_args=vsc_min_fraction_snps=0.005,vsc_min_count_snps=2",)}


def _docker(d, image, args):
    cmd = ["docker", "--context", CONTEXT, "run", "--rm", "--platform", "linux/amd64", "-v", f"{d}:{d}", image, *args]
    return subprocess.run(cmd, capture_output=True, text=True)


def _read_vcf(vcf, lo, hi, germline, sample_col=9):
    calls = {}
    for line in gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pos, rb, filt = int(f[1]), f[3], f[6]
        if filt != "PASS" or len(rb) != 1 or not lo <= pos <= hi:
            continue
        fmt = dict(zip(f[8].split(":"), f[sample_col].split(":"), strict=False))
        if fmt.get("GT", "0/0").replace("|", "/") in ("0/0", "./.", "0"):
            continue
        afs = [float(x) for x in fmt.get("VAF", "nan").split(",")]
        for alt, af in zip(f[4].split(","), afs, strict=False):
            if len(alt) == 1 and alt != "*" and af < germline:
                calls[(pos, rb, alt)] = af
    return calls


def deepvariant_calls(reads, reference, lo, hi, setting="deepvariant", germline=0.5, threads=8, keep_vcf=None):
    """`keep_vcf`, a path: also save the full VCF there."""
    mm2 = _mm2()
    with tempfile.TemporaryDirectory(dir=os.path.expanduser("~/Nanoshm")) as d:
        ref, name = _ref(reference, d)
        _bam(reads, ref, os.path.join(d, "s.bam"), mm2, "sample")
        vcf = os.path.join(d, "out.vcf.gz")
        r = _docker(d, DV_IMAGE, ["/opt/deepvariant/bin/run_deepvariant", "--model_type=ONT_R104", f"--ref={ref}",
                                  f"--reads={d}/s.bam", f"--output_vcf={vcf}", f"--regions={name}:{lo}-{hi}",
                                  f"--num_shards={threads}", f"--intermediate_results_dir={d}/tmp",
                                  *DEEPVARIANT[setting]])
        if r.returncode != 0 or not os.path.exists(vcf):
            raise RuntimeError(f"DeepVariant failed ({r.returncode}): {(r.stderr or r.stdout)[-1500:]}")
        if keep_vcf:
            shutil.copy(vcf, keep_vcf)
        return _read_vcf(vcf, lo, hi, germline)


def deepsomatic_calls(tumor_reads, normal_reads, reference, lo, hi, setting="deepsomatic", germline=0.5,
                      threads=8, keep_vcf=None):
    """`keep_vcf`, a path: also save the full VCF there (every candidate, with its FILTER)."""
    mm2 = _mm2()
    with tempfile.TemporaryDirectory(dir=os.path.expanduser("~/Nanoshm")) as d:
        ref, name = _ref(reference, d)
        _bam(tumor_reads, ref, os.path.join(d, "tumor.bam"), mm2, "tumor")
        _bam(normal_reads, ref, os.path.join(d, "normal.bam"), mm2, "normal")
        vcf = os.path.join(d, "out.vcf.gz")
        r = _docker(d, DS_IMAGE, ["run_deepsomatic", "--model_type=ONT", f"--ref={ref}",
                                  f"--reads_tumor={d}/tumor.bam", f"--reads_normal={d}/normal.bam",
                                  "--sample_name_tumor=tumor", "--sample_name_normal=normal",
                                  f"--output_vcf={vcf}", f"--regions={name}:{lo}-{hi}", f"--num_shards={threads}",
                                  f"--intermediate_results_dir={d}/tmp", *DEEPSOMATIC[setting]])
        if r.returncode != 0 or not os.path.exists(vcf):
            raise RuntimeError(f"DeepSomatic failed ({r.returncode}): {(r.stderr or r.stdout)[-1500:]}")
        if keep_vcf:
            shutil.copy(vcf, keep_vcf)
        return _read_vcf(vcf, lo, hi, germline)
