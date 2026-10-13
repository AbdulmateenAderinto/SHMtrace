"""Single-read estimate against planted one-read mutations, on the delivered J_H4 control reads.

    python docs/calibration/crm/single_read_spike.py 0,25,50,100 3

As docs/calibration/single_read_spike.py, but each held-out Aicda-/- library is first prepared by
SHMtrace and 500 of its usable reads drawn (seed 0), and the N one-read mutations are planted into
those 500, so every planted mutation is analysed. The other two controls (delivered reads) are the
background. Writes out/single_read_spike.json.
"""
import json
import random
import sys
import tempfile
from pathlib import Path

from shmtrace import pipeline, project, simulate
from shmtrace.io_utils import read_fasta, read_fastq
from shmtrace.project import build_params

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402

DS = DATASETS["jh4"]
REF, CON = DS["reference"], DS["contaminants"]
KOS = ["KO1", "KO2", "KO3"]
NS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "0,25,50,100").split(",")]
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
p = build_params({"region": "jh4_intron"})
ref = read_fasta(REF)[0][1]
lo, hi = p.analysis_window[0] - 1, p.analysis_window[1] - 1
tmp = Path(tempfile.mkdtemp())


def run(path, peers):
    proj = {"name": "s", "region": "jh4_intron", "reference": REF, "contaminants": CON, "depth": 500,
            "samples": [{"name": "S", "path": str(path), "group": "S"}]
            + [{"name": c, "path": DS["samples"][c]["fastq"], "group": "KO", "control": True} for c in peers]}
    return project.run_project(proj)["samples"][0]


rows = []
for held in KOS:
    peers = [c for c in KOS if c != held]
    fq = DS["samples"][held]["fastq"]
    usable, _ = pipeline.prepare_reads(fq, REF, CON, p)
    keep = {a.qname for a in usable}
    reads = random.Random(0).sample([r for r in read_fastq(fq) if r[0].split()[0] in keep], 500)
    base_f = tmp / f"{held}_base.fastq"
    with open(base_f, "w") as fh:
        for n, s, q in reads:
            fh.write(f"@{n}\n{s}\n+\n{q}\n")
    base = run(base_f, peers)
    for n in NS:
        for rep in range(REPS if n else 1):
            rng = random.Random(1000 * rep + n)
            sites = simulate.choose_sites(ref, lo, hi, n, rng, read_counts=(1,))
            spiked, planted = simulate.spike(reads, REF, sites, p, rng)
            f = tmp / f"{held}_{n}_{rep}.fastq"
            with open(f, "w") as fh:
                for name, seq, qual in spiked:
                    fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
            r = run(f, peers)
            row = {"held": held, "planted": sum(planted.values()), "rep": rep, "reads": r["reads"],
                   "base_estimate": base["single_read_estimate"], "estimate": r["single_read_estimate"],
                   "delta": r["single_read_estimate"] - base["single_read_estimate"],
                   "ci_low": r.get("single_read_density_ci_low"), "ci_high": r.get("single_read_density_ci_high")}
            rows.append(row)
            print(row, flush=True)
out = Path(__file__).parent / "out"
out.mkdir(exist_ok=True)
json.dump(rows, open(out / "single_read_spike.json", "w"), indent=1)
