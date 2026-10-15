"""ClairS (or DeepSomatic) on the real AID-active libraries of a dataset (J_H4 by default), with its AID-null
libraries pooled as the matched normal.

    python docs/calibration/crm/clairs_real.py                  -> out/clairs_real_jh4.json
    python docs/calibration/crm/clairs_real.py --deepsomatic    -> out/deepsomatic_real_jh4.json
    python docs/calibration/crm/clairs_real.py --deepsomatic --dataset n114   -> out/deepsomatic_real_n114.json

Each Wt library: the 500 usable reads SHMtrace analyses (seed 0), against the Aicda-/- libraries' usable
reads pooled as the normal. Records every candidate ClairS wrote, with its FILTER, QUAL, allele fraction and
read support in the case and in the normal, so the paper can say how many candidates ClairS considered and why
it passed few of them. SHMtrace's call count on the same reads is from the project run (benchmark_real).
"""
import gzip
import json
import os
import random
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from clairs_runner import CLAIRS, clairs_calls  # noqa: E402
from datasets import DATASETS  # noqa: E402
from deepvariant_runner import DEEPSOMATIC, deepsomatic_calls  # noqa: E402

from shmtrace import pipeline  # noqa: E402
from shmtrace.config import load_preset  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

OUT = Path(__file__).parent / "out"


def usable(fq, ds, p):
    u, _ = pipeline.prepare_reads(fq, ds["reference"], ds.get("contaminants", []), p)
    names = {a.qname for a in u}
    return [r for r in read_fastq(fq) if r[0].split()[0] in names]


if __name__ == "__main__":
    somatic = "--deepsomatic" in sys.argv          # DeepSomatic, the other tumor-normal caller
    settings, call, name = ((DEEPSOMATIC, deepsomatic_calls, "deepsomatic") if somatic else
                            (CLAIRS, clairs_calls, "clairs"))
    key = sys.argv[sys.argv.index("--dataset") + 1] if "--dataset" in sys.argv else "jh4"
    ds = DATASETS[key]
    cases = [k for k, s in ds["samples"].items() if not s["control"]]
    controls = [k for k, s in ds["samples"].items() if s["control"]]
    p = load_preset(ds["preset"])
    ref = read_fasta(ds["reference"])[0][1]
    lo, hi = pipeline._window(p, ref)
    normal = [r for k in controls for r in usable(ds["samples"][k]["fastq"], ds, p)]
    out = {"normal_reads": len(normal), "window": [lo + 1, hi + 1], "libraries": {}}
    for wt in cases:
        reads = usable(ds["samples"][wt]["fastq"], ds, p)
        reads = random.Random(0).sample(reads, min(500, len(reads)))
        res = {}
        for setting in settings:
            with tempfile.TemporaryDirectory() as d:
                vcf = os.path.join(d, "c.vcf.gz")
                calls = call(reads, normal, ds["reference"], lo + 1, hi + 1, setting=setting, keep_vcf=vcf)
                cands = []
                for line in gzip.open(vcf, "rt"):
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    fmt = dict(zip(f[8].split(":"), f[9].split(":"), strict=False))
                    ad = fmt.get("AD", "0,0").split(",")
                    cands.append({"pos": int(f[1]), "ref": f[3], "alt": f[4], "qual": float(f[5]), "filter": f[6],
                                  "af": max(float(x) for x in fmt.get("AF", fmt.get("VAF", "nan")).split(",")),
                                  "alt_reads": max(int(x) for x in ad[1:]) if all(x.isdigit() for x in ad[1:]) and
                                  len(ad) > 1 else None,
                                  "normal_af": float(fmt.get("NAF", "nan"))})
            res[setting] = {"pass_somatic": len(calls), "candidates": cands}
            print(wt, setting, f"{len(cands)} candidates, {len(calls)} passed below 50% AF", flush=True)
        out["libraries"][wt] = {"reads": len(reads), **res}
    OUT.mkdir(exist_ok=True)
    json.dump(out, open(OUT / f"{name}_real_{key}.json", "w"), indent=1)
    print(OUT / f"{name}_real_{key}.json")
