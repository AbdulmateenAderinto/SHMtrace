"""ClairS (Zheng et al. 2026, somatic small variants from long reads) on the reads SHMtrace analysed.

    from clairs_runner import clairs_calls
    calls = clairs_calls(tumor_reads, normal_reads, reference, lo, hi, setting="clairs")

ClairS is the somatic caller closest to SHMtrace's design: it compares a case library ("tumor") with
a matched normal, as SHMtrace compares a sample with its AID-null controls. Here the case is the
library being tested and the normal is the other controls' usable reads, pooled - the same libraries
SHMtrace's null is fitted from. Reads are aligned with minimap2 (map-ont, as SHMtrace), sorted and
indexed with samtools, and called with the official hkubal/clairs Docker image, restricted to the
analysed window. Substitutions only, FILTER PASS, and sites at >= 50% allele fraction are removed,
as for every caller in the benchmark. Returns {(pos1, ref, alt): af}.

Settings:
  clairs             ClairS at its defaults for ONT R10.4.1 SUP 5 kHz reads (snv_min_af 0.05)
  clairs_permissive  snv_min_af 0.005 and minimum coverage 4: the lowest allele fraction ClairS is
                     asked to report, to give it the best chance at mutations carried by few reads
"""
import gzip
import os
import subprocess
import tempfile

BIN = os.path.expanduser("~/miniforge3/envs/lofreq/bin")          # minimap2, samtools
IMAGE = "hkubal/clairs:latest"           # linux/amd64: emulated on Apple Silicon, so slow
CONTEXT = "desktop-linux"
PLATFORM = "ont_r10_dorado_sup_5khz"     # dna_r10.4.1_e8.2_400bps_sup@v5.x, 5 kHz
CLAIRS = {"clairs": (), "clairs_permissive": ("--snv_min_af", "0.005", "--min_coverage", "4")}


def _mm2():
    for c in (os.path.join(BIN, "minimap2"), "minimap2"):
        try:
            subprocess.run([c, "--version"], capture_output=True, check=True)
            return c
        except (OSError, subprocess.CalledProcessError):
            continue
    raise RuntimeError("minimap2 not found")


def _bam(reads, ref, path, mm2, sample):
    fq = path + ".fastq"
    with open(fq, "w") as fh:
        for n, s, q in reads:
            fh.write(f"@{n}\n{s}\n+\n{q}\n")
    subprocess.run(f"{mm2} -a -x map-ont --secondary=no -R '@RG\\tID:{sample}\\tSM:{sample}' {ref} {fq} "
                   f"2>/dev/null | {BIN}/samtools sort -o {path} - && {BIN}/samtools index {path}",
                   shell=True, check=True)


def _ref(reference, d):
    """The first record of `reference`, named by the first word of its header, written to d/ref.fa and indexed."""
    ref = os.path.join(d, "ref.fa")
    name, seq = None, []
    for line in open(reference):
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    with open(ref, "w") as fh:
        fh.write(f">{name}\n{''.join(seq)}\n")
    subprocess.run([f"{BIN}/samtools", "faidx", ref], check=True)
    return ref, name


def clairs_calls(tumor_reads, normal_reads, reference, lo, hi, setting="clairs", germline=0.5, threads=8,
                 keep_vcf=None):
    """`keep_vcf`, a path: also save ClairS's full VCF there (every candidate, with its FILTER)."""
    mm2 = _mm2()
    with tempfile.TemporaryDirectory(dir=os.path.expanduser("~/Nanoshm")) as d:
        ref, name = _ref(reference, d)
        _bam(tumor_reads, ref, os.path.join(d, "tumor.bam"), mm2, "tumor")
        _bam(normal_reads, ref, os.path.join(d, "normal.bam"), mm2, "normal")
        out = os.path.join(d, "out")
        cmd = ["docker", "--context", CONTEXT, "run", "--rm", "--platform", "linux/amd64", "-v", f"{d}:{d}", IMAGE,
               "/opt/bin/run_clairs", "--tumor_bam_fn", f"{d}/tumor.bam", "--normal_bam_fn", f"{d}/normal.bam",
               "--ref_fn", ref, "--threads", str(threads), "--platform", PLATFORM, "--output_dir", out,
               "--region", f"{name}:{lo}-{hi}", "--include_all_ctgs", "--conda_prefix", "/opt/conda/envs/clairs",
               *CLAIRS[setting]]
        r = subprocess.run(cmd, capture_output=True, text=True)
        vcf = next((os.path.join(out, f) for f in ("snv.vcf.gz", "output.vcf.gz") if os.path.exists(os.path.join(out, f))),
                   None)
        if r.returncode != 0 or vcf is None:
            raise RuntimeError(f"ClairS failed ({r.returncode}): {(r.stderr or r.stdout)[-1500:]}")
        if keep_vcf:
            import shutil
            shutil.copy(vcf, keep_vcf)
        calls = {}
        for line in gzip.open(vcf, "rt"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos, rb, alt, filt = int(f[1]), f[3], f[4], f[6]
            if filt != "PASS" or len(rb) != 1 or len(alt) != 1 or not lo <= pos <= hi:
                continue
            fmt = dict(zip(f[8].split(":"), f[9].split(":"), strict=False))
            af = float(fmt.get("AF", "nan"))
            if af >= germline:
                continue
            calls[(pos, rb, alt)] = af
        return calls


# Clair3 (Zheng et al. 2022), the germline caller bundled in the same image, with the model for
# dna_r10.4.1_e8.2_400bps_sup@v5.2: a germline caller looks for genotypes (variants in about half or all
# of the reads), so it is the reference point for what a germline model does with SHM.
#   clair3             its defaults for ONT (snp_min_af 0.08)
#   clair3_permissive  snp_min_af 0.01, min_coverage 4
CLAIR3_MODEL = "/opt/conda/envs/clairs/bin/models/r1041_e82_400bps_sup_v520"
CLAIR3 = {"clair3": (), "clair3_permissive": ("--snp_min_af=0.01", "--min_coverage=4")}


def clair3_calls(reads, reference, lo, hi, setting="clair3", germline=0.5, threads=8):
    mm2 = _mm2()
    with tempfile.TemporaryDirectory(dir=os.path.expanduser("~/Nanoshm")) as d:
        ref, name = _ref(reference, d)
        _bam(reads, ref, os.path.join(d, "s.bam"), mm2, "sample")
        bed = os.path.join(d, "w.bed")
        with open(bed, "w") as fh:
            fh.write(f"{name}\t{lo - 1}\t{hi}\n")
        out = os.path.join(d, "out")
        cmd = ["docker", "--context", CONTEXT, "run", "--rm", "--platform", "linux/amd64", "-v", f"{d}:{d}", IMAGE,
               "/opt/conda/envs/clairs/bin/run_clair3.sh", f"--bam_fn={d}/s.bam", f"--ref_fn={ref}",
               f"--model_path={CLAIR3_MODEL}", f"--threads={threads}", "--platform=ont", f"--output={out}",
               f"--bed_fn={bed}", "--include_all_ctgs", *CLAIR3[setting]]
        r = subprocess.run(cmd, capture_output=True, text=True)
        vcf = os.path.join(out, "merge_output.vcf.gz")
        if r.returncode != 0 or not os.path.exists(vcf):
            raise RuntimeError(f"Clair3 failed ({r.returncode}): {(r.stderr or r.stdout)[-1500:]}")
        calls = {}
        for line in gzip.open(vcf, "rt"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos, rb, filt = int(f[1]), f[3], f[6]
            fmt = dict(zip(f[8].split(":"), f[9].split(":"), strict=False))
            gt = fmt.get("GT", "0/0").replace("|", "/")
            if filt != "PASS" or gt in ("0/0", "./.") or len(rb) != 1 or not lo <= pos <= hi:
                continue
            afs = [float(x) for x in fmt.get("AF", "nan").split(",")]
            for alt, af in zip(f[4].split(","), afs, strict=False):
                if len(alt) == 1 and af < germline:
                    calls[(pos, rb, alt)] = af
        return calls
