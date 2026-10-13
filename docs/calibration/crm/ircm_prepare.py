"""Prepare the IRCM public nanopore runs (PRJNA1481264) for the benchmarks.

    python docs/calibration/crm/ircm_prepare.py

1. Each library is subsampled at random to 6,000 reads (seed 0): sub/<library>.fastq.
2. Each locus gets a reference: the majority consensus of its first Aicda-/- library, started from a
   high-quality read of the median length and iterated (minimap2, `samtools consensus`) until it no
   longer changes. A consensus of knockout reads carries no somatic mutation by construction, which
   is what the reference has to be; there is no published sequence for these amplicons.
Writes references/<locus>.fasta. The data stay in ~/Nanoshm/public_data/PRJNA1481264 (not the repo).
"""
import gzip
import random
import statistics
import subprocess
from pathlib import Path

ROOT = Path.home() / "Nanoshm/public_data/PRJNA1481264"
BIN = Path.home() / "miniforge3/envs/lofreq/bin"
N = 6000
LOCI = {"Igh_Su": "Igh_Su_LPS_AID_KO", "Pim1": "Pim1_LPS_AID_KO"}


def reads(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, _p, q = fh.readline().strip(), fh.readline(), fh.readline().strip()
            yield h[1:].split()[0], s, q


def subsample():
    (ROOT / "sub").mkdir(exist_ok=True)
    for f in sorted((ROOT / "fastq").glob("*.fastq.gz")):
        out = ROOT / "sub" / (f.name[:-3])
        if out.exists():
            continue
        allr = list(reads(f))
        pick = random.Random(0).sample(allr, min(N, len(allr)))
        with open(out, "w") as fh:
            for n, s, q in pick:
                fh.write(f"@{n}\n{s}\n+\n{q}\n")
        print(f.name, len(allr), "->", len(pick), flush=True)


def consensus(locus, lib):
    (ROOT / "references").mkdir(exist_ok=True)
    fq = ROOT / "sub" / f"{lib}.fastq"
    rs = [(n, s, q) for n, s, q in reads_plain(fq)]
    med = statistics.median(len(s) for _n, s, _q in rs)
    meanq = lambda q: sum(ord(c) - 33 for c in q) / len(q)   # noqa: E731
    seed = max((r for r in rs if abs(len(r[1]) - med) <= 2), key=lambda r: meanq(r[2]))
    ref = ROOT / "references" / f"{locus}.fasta"
    ref.write_text(f">{locus}\n{seed[1]}\n")
    for it in range(8):
        bam = ROOT / "references" / f"{locus}.bam"
        subprocess.run(f"minimap2 -a -x map-ont --secondary=no {ref} {fq} 2>/dev/null | "
                       f"{BIN}/samtools sort -o {bam} - && {BIN}/samtools index {bam}", shell=True, check=True)
        out = subprocess.run([f"{BIN}/samtools", "consensus", "-a", "--show-del", "no", "--show-ins", "yes",
                              "-m", "simple", "-c", "0.5", "-f", "fasta", str(bam)],
                             capture_output=True, text=True, check=True).stdout
        seq = "".join(out.split("\n")[1:]).upper().replace("N", "")
        old = "".join(ref.read_text().split("\n")[1:])
        ref.write_text(f">{locus}\n{seq}\n")
        print(locus, "iteration", it, len(seq), "changed" if seq != old else "converged", flush=True)
        if seq == old:
            break


def reads_plain(path):
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, _p, q = fh.readline().strip(), fh.readline(), fh.readline().strip()
            yield h[1:].split()[0], s, q


if __name__ == "__main__":
    subsample()
    for locus, lib in LOCI.items():
        if (ROOT / "sub" / f"{lib}.fastq").exists():
            consensus(locus, lib)
