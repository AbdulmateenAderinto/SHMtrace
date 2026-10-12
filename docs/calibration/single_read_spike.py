"""Spike-in check of the single-read estimate: plant N one-read mutations into a
held-out control and see whether the estimate rises by N.

    python docs/calibration/single_read_spike.py 0,25,50,100 2
"""
import random
import sys
import tempfile
from pathlib import Path

from shmtrace import project, simulate
from shmtrace.io_utils import read_fasta, read_fastq
from shmtrace.project import build_params

ROOT = Path(__file__).resolve().parents[2]
REF = str(ROOT / "src/shmtrace/references/mouse_jh4_intron.fasta")
CON = str(ROOT / "src/shmtrace/references/mouse_jh4_contaminant.fasta")
EX = ROOT / "examples_compare"
KOS = ["KO1", "KO2", "KO3"]
NS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "0,25,50,100").split(",")]
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 2

p = build_params({"region": "jh4_intron"})
ref_name, ref = read_fasta(REF)[0]
lo, hi = p.analysis_window[0] - 1, p.analysis_window[1] - 1
tmp = Path(tempfile.mkdtemp())


def run(sample_path, held, peers):
    proj = {"name": "s", "region": "jh4_intron", "reference": REF, "contaminants": [CON],
            "samples": [{"name": "S", "path": str(sample_path), "group": "S"}]
            + [{"name": c, "path": str(EX / f"{c}.fastq"), "group": "KO", "control": True}
               for c in peers]}
    r = project.run_project(proj)["samples"][0]
    return r


print("held planted  base_est  spiked_est  delta  ci_width  reads_kept")
for held in KOS:
    peers = [c for c in KOS if c != held]
    reads = list(read_fastq(str(EX / f"{held}.fastq")))
    base = run(EX / f"{held}.fastq", held, peers)
    for n in NS:
        for rep in range(REPS if n else 1):
            rng = random.Random(1000 * rep + n)
            sites = simulate.choose_sites(ref, lo, hi, n, rng, read_counts=(1,))
            spiked, planted = simulate.spike(reads, REF, sites, p, rng)
            f = tmp / f"{held}_{n}_{rep}.fastq"
            with open(f, "w") as fh:
                for name, seq, qual in spiked:
                    fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
            r = run(f, held, peers)
            got = sum(planted.values())
            kept = r["reads"] / len(reads)
            d = r["single_read_estimate"] - base["single_read_estimate"]
            print(f"{held} {got:>7} {base['single_read_estimate']:>9} "
                  f"{r['single_read_estimate']:>11} {d:>6.1f} "
                  f"{(r['single_read_density_ci_high'] - r['single_read_density_ci_low']) / 2:>8.2f}"
                  f" {kept:.3f}  expected delta ~{got * kept:.1f}", flush=True)
