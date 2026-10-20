"""Real hypermutation diluted into AID-null reads: is the estimate linear, and where does it stop being detectable?

    python crm/titration.py [reps]       -> crm/out/titration.json

The spike-ins of Figure 3 plant mutations at chosen sites; this uses the mutations real AID made, with their
clonality, hotspot clustering and multi-mutation reads. For each AID-active library (the donor), a library of 500
usable reads is drawn with a fraction f from the donor and the rest from one AID-null library (the diluent), and
called by `project.run_project` against the other two AID-null libraries. Every read is a usable read of its own
library (SHMtrace's read preparation), so f is the share of analysed reads, not of delivered ones.

  J_H4   donors WT2, WT3, WT4; diluent KO1; controls KO2, KO3
  N114   donors C1, D3, D5;   diluent D1;  controls C3, C4

N114's reference is the parental construct, so it runs with `ancestral_reference` (the opt-in exception to the
50% germline rule); J_H4 runs at its defaults.

Mutation frequency is read-weighted, so its expectation is linear in f: f x (the donor at f = 1) + (1 - f) x (the
diluent at f = 0). Unique-variant density counts each mutation once and is not expected to be linear when clones
share mutations; it is reported beside it.
"""
import json
import os
import random
import sys
import tempfile
import zlib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from datasets import DATASETS  # noqa: E402

from shmtrace import pipeline, project  # noqa: E402
from shmtrace.config import load_preset  # noqa: E402
from shmtrace.io_utils import read_fastq  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
FRACTIONS = (0.0, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.0)
DEPTH = 500
DESIGN = {"jh4": {"region": "jh4_intron", "donors": ("WT2", "WT3", "WT4"), "diluent": "KO1",
                  "controls": ("KO2", "KO3")},
          # the reference is the parental construct, and mutations fixed within a subclone (C1 at 163, D3 at 94
          # and 106) are real: with the default 50% rule they are removed at f = 1 but counted once diluted
          "n114": {"region": "n114_vregion", "donors": ("C1", "D3", "D5"), "diluent": "D1",
                   "controls": ("C3", "C4"), "overrides": {"ancestral_reference": True}}}
METRICS = ("unique_mutations", "unique_density_x1e5", "mutation_frequency_x1e3", "mutated_reads_pct")


def usable(key, name):
    ds = DATASETS[key]
    base = replace(load_preset(ds["preset"]), **ds.get("overrides", {}))
    fq = ds["samples"][name]["fastq"]
    keep = {a.qname for a in pipeline.prepare_reads(fq, ds["reference"], ds.get("contaminants", []), base)[0]}
    return [r for r in read_fastq(fq) if r[0].split()[0] in keep]


def one(job):
    key, donor, f, rep, donor_reads, diluent_reads = job
    d, ds = DESIGN[key], DATASETS[key]
    rng = random.Random(zlib.crc32(f"{key}|{donor}|{f}|{rep}".encode()))
    k = round(f * DEPTH)
    reads = rng.sample(donor_reads, k) + rng.sample(diluent_reads, DEPTH - k)
    rng.shuffle(reads)
    fd, path = tempfile.mkstemp(suffix=".fastq")
    with os.fdopen(fd, "w") as fh:
        for n, s, q in reads:
            fh.write(f"@{n}\n{s}\n+\n{q}\n")
    proj = {"name": "titration", "region": d["region"], "reference": ds["reference"],
            "contaminants": ds.get("contaminants", []), "mode": "recommended", "depth": DEPTH,
            "overrides": d.get("overrides", {}),
            "samples": [{"name": "MIX", "path": path, "group": "MIX"}]
                       + [{"name": c, "path": ds["samples"][c]["fastq"], "group": "CTRL", "control": True}
                          for c in d["controls"]]}
    try:
        out = project.run_project(proj)
    finally:
        os.remove(path)
    row = next(s for s in out["samples"] if s["sample"] == "MIX")
    return {"locus": key, "donor": donor, "fraction": f, "rep": rep, "donor_reads": k,
            **{m: row.get(m) for m in METRICS}}


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    jobs = []
    for key, d in DESIGN.items():
        dil = usable(key, d["diluent"])
        for donor in d["donors"]:
            don = usable(key, donor)
            print(key, donor, len(don), "usable;", d["diluent"], len(dil), flush=True)
            jobs += [(key, donor, f, rep, don, dil) for f in FRACTIONS for rep in range(reps)]
    rows = []
    with ProcessPoolExecutor(6) as ex:
        for i, r in enumerate(ex.map(one, jobs), 1):
            rows.append(r)
            print(f"[{i}/{len(jobs)}] {r['locus']} {r['donor']} f={r['fraction']} rep {r['rep']}: "
                  f"{r['unique_mutations']} unique, freq {r['mutation_frequency_x1e3']}", flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump({"fractions": FRACTIONS, "depth": DEPTH, "reps": reps, "design": DESIGN, "rows": rows},
              open(OUT / "titration.json", "w"), indent=1)


if __name__ == "__main__":
    main()
