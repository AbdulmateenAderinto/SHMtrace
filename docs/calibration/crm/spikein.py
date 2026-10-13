"""Spike-in recall and false calls for SHMtrace's modes and LoFreq, on held-out AID-null controls.

    python docs/calibration/crm/spikein.py jh4 n114 [--reps 3 --sites 60]

For each control in turn: its usable reads (SHMtrace's read preparation) are drawn to the dataset's
depth, 60 AID-weighted substitutions are planted into 1, 2, 3, 5, 10 or 20 reads each with the
read's own base quality kept (shmtrace.simulate), and the library is called by
  recommended, geneious, v0.1, pbinom   SHMtrace modes; the null from the other controls only
  lofreq                                LoFreq 2.1.5 at its defaults on the same reads
Recall by planted read count, false calls per library and FDR. Writes out/spikein_<dataset>.json.
"""
import argparse
import json
import os
import random
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

from shmtrace import null_model, pipeline, simulate
from shmtrace.config import load_preset
from shmtrace.io_utils import read_fasta, read_fastq

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402
from lofreq_runner import LOFREQ, lofreq_calls  # noqa: E402

OUT = Path(__file__).parent / "out"
PROFILES = ("recommended", "geneious", "v0.1", "pbinom")
PROFILES_RUN, WITH_LOFREQ, TAG = list(PROFILES), True, ""
# extra settings tried on top of a mode: name -> (mode, Params overrides)
VARIANTS = {"recommended_bhall": ("recommended", {"bh_family": "all"})}


def run(key, reps, n_sites, seed=1):
    ds = DATASETS[key]
    base = replace(load_preset(ds["preset"]), **ds.get("overrides", {}))
    ref = read_fasta(ds["reference"])[0][1]
    lo, hi = pipeline._window(base, ref)
    lo0, hi0 = max(lo, base.end_trim_bp), min(hi, len(ref) - 1 - base.end_trim_bp)
    contam = ds.get("contaminants", [])
    controls = [n for n, s in ds["samples"].items() if s["control"]]
    rows = []
    for ctl in controls:
        peers = [ds["samples"][c]["fastq"] for c in controls if c != ctl]
        fq = ds["samples"][ctl]["fastq"]
        usable, _ = pipeline.prepare_reads(fq, ds["reference"], contam, base)
        names = {a.qname for a in usable}
        pool = [r for r in read_fastq(fq) if r[0].split()[0] in names]
        nulls = {}
        rng = random.Random(seed)
        for rep in range(reps):
            sub = rng.sample(pool, min(len(pool), ds["depth"]))
            sites = simulate.choose_sites(ref, lo0, hi0, n_sites, rng)
            spiked, planted = simulate.spike(sub, ds["reference"], sites, base, rng)
            with tempfile.NamedTemporaryFile("w", suffix=".fastq", delete=False) as fh:
                for n, s, q in spiked:
                    fh.write(f"@{n}\n{s}\n+\n{q}\n")
                tmp = fh.name
            try:
                for prof in PROFILES_RUN:
                    mode, extra = VARIANTS.get(prof, (prof, {}))
                    pp = replace(base, **pipeline.mode_overrides(mode), **extra)
                    nr = None
                    if pp.error_model != "geneious":
                        k = (pp.error_model, pp.phase_mode)
                        if k not in nulls:
                            nulls[k] = null_model.fit(peers, ds["reference"], contam, pp)[0]
                        nr = nulls[k]
                    res = pipeline.run(tmp, ds["reference"], contam, replace(pp, downsample_to=None), nr)
                    called = {(v.pos - 1, v.alt) for v in res["variants"]}
                    rows.append({"control": ctl, "rep": rep, "caller": prof, "reads": len(res["_aligned"]),
                                 **simulate.evaluate(called, planted)})
                    print(key, ctl, rep, prof, rows[-1]["recall"], rows[-1]["false_calls"], flush=True)
                    if prof == "recommended" and WITH_LOFREQ:
                        keep = {a.qname for a in res["_aligned"]}
                        reads = [r for r in spiked if r[0].split()[0] in keep]
                        for lname, extra in LOFREQ.items():
                            lf = lofreq_calls(reads, ds["reference"], lo, hi,
                                              preset=base.minimap2_preset, extra=extra)
                            called = {(pos - 1, a) for pos, _r, a in lf}
                            rows.append({"control": ctl, "rep": rep, "caller": lname, "reads": len(reads),
                                         **simulate.evaluate(called, planted)})
                            print(key, ctl, rep, lname, rows[-1]["recall"], rows[-1]["false_calls"],
                                  flush=True)
            finally:
                os.remove(tmp)
    OUT.mkdir(exist_ok=True)
    json.dump({"dataset": key, "reps": reps, "sites": n_sites, "rows": rows},
              open(OUT / f"spikein_{key}{TAG}.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="+")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--sites", type=int, default=60)
    ap.add_argument("--profiles", default=",".join(PROFILES), help="modes or VARIANTS names")
    ap.add_argument("--no-lofreq", action="store_true")
    ap.add_argument("--tag", default="", help="suffix for the output file")
    a = ap.parse_args()
    PROFILES_RUN = a.profiles.split(",")
    WITH_LOFREQ = not a.no_lofreq
    TAG = a.tag
    for k in a.datasets:
        run(k, a.reps, a.sites)
