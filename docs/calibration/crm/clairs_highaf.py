"""Does ClairS call clean, high-allele-fraction mutations planted into AID-null reads?

    python docs/calibration/crm/clairs_highaf.py                  -> out/clairs_highaf_jh4.json
    python docs/calibration/crm/clairs_highaf.py --deepsomatic    -> out/deepsomatic_highaf_jh4.json

The spike-in benchmark plants mutations into 1-20 of 500 reads (at most 4% allele fraction), below ClairS's
usual range. On real wild-type J_H4 libraries, ClairS also rejected mutations at 7-29%. This separates the two
explanations: 30 AID-weighted substitutions are planted into 25, 50 and 100 of 500 reads (5, 10 and 20%) of each
held-out Aicda-/- library, which otherwise carries few mutations, and ClairS is run with the other two libraries as
its normal. If it calls these, its failure on real SHM comes from the hypermutated reads, not the allele fraction.
"""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from clairs_runner import CLAIRS, clairs_calls  # noqa: E402
from datasets import DATASETS  # noqa: E402
from deepvariant_runner import DEEPSOMATIC, deepsomatic_calls  # noqa: E402

from shmtrace import pipeline, simulate  # noqa: E402
from shmtrace.config import load_preset  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

OUT = Path(__file__).parent / "out"
COUNTS = (25, 50, 100)


def usable(fq, ds, p):
    u, _ = pipeline.prepare_reads(fq, ds["reference"], ds["contaminants"], p)
    names = {a.qname for a in u}
    return [r for r in read_fastq(fq) if r[0].split()[0] in names]


if __name__ == "__main__":
    somatic = "--deepsomatic" in sys.argv          # DeepSomatic, the other tumor-normal caller, on the same plants
    settings, call, name = ((DEEPSOMATIC, deepsomatic_calls, "deepsomatic") if somatic else
                            (CLAIRS, clairs_calls, "clairs"))
    ds = DATASETS["jh4"]
    p = load_preset(ds["preset"])
    ref = read_fasta(ds["reference"])[0][1]
    lo, hi = pipeline._window(p, ref)
    kos = ("KO1", "KO2", "KO3")
    reads = {k: usable(ds["samples"][k]["fastq"], ds, p) for k in kos}
    rows = []
    for k in kos:
        rng = random.Random(7)
        sub = rng.sample(reads[k], 500)
        sites = simulate.choose_sites(ref, lo, hi, 30, rng, read_counts=COUNTS)
        spiked, planted = simulate.spike(sub, ds["reference"], sites, p, rng)
        normal = [r for o in kos if o != k for r in reads[o]]
        for setting in settings:
            calls = call(spiked, normal, ds["reference"], lo + 1, hi + 1, setting=setting)
            called = {(pos - 1, a) for pos, _r, a in calls}
            ev = simulate.evaluate(called, planted)
            rows.append({"control": k, "caller": setting, **ev})
            print(k, setting, ev["recall_by_reads"], "false", ev["false_calls"], flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump({"counts": COUNTS, "rows": rows}, open(OUT / f"{name}_highaf_jh4.json", "w"), indent=1)
    print(OUT / f"{name}_highaf_jh4.json")
