"""SHMtrace against Geneious's rule and LoFreq on real libraries with matched AID-null controls.

    python docs/calibration/crm/benchmark_real.py jh4 n114 su_lps su_ch12 pim1

For each dataset, every library is prepared once by SHMtrace (read gates, competitive alignment,
paralogue filter) and the same analysed reads (a random draw at the dataset's depth) are called by:
  recommended  SHMtrace's control-calibrated test (each control called against the other controls)
  geneious     Geneious's approximate per-candidate p <= 0.05, uncorrected (SHMtrace's geneious mode)
  lofreq       LoFreq 2.1.5 `call` at its defaults on minimap2 alignments of the same reads
Every caller is restricted to the analysed window, substitutions only, sites at >= 50% removed.
Densities are distinct calls per analysed base (x 1e5) with one denominator per library.
Writes docs/calibration/crm/out/benchmark_real_<dataset>.json.
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

from shmtrace import batch, null_model, pipeline
from shmtrace.config import load_preset
from shmtrace.io_utils import read_fasta, read_fastq

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402
from lofreq_runner import LOFREQ, lofreq_calls  # noqa: E402

OUT = Path(__file__).parent / "out"
MODES = ("recommended", "geneious")


def window_bases(aligned, lo, hi):
    return sum(1 for a in aligned for rpos, base, _q in a.pairs if base and lo - 1 <= rpos <= hi - 1)


def run(key):
    ds = DATASETS[key]
    base = load_preset(ds["preset"])
    base = replace(base, **ds.get("overrides", {}))
    ref_seq = read_fasta(ds["reference"])[0][1]
    ref_name = read_fasta(ds["reference"])[0][0]
    lo, hi = pipeline._window(base, ref_seq)
    if base.end_trim_bp:
        lo, hi = max(lo, base.end_trim_bp + 1), min(hi, len(ref_seq) - base.end_trim_bp)
    samples = {n: s["fastq"] for n, s in ds["samples"].items()}
    controls = {n for n, s in ds["samples"].items() if s["control"]}
    records = {}
    rows = []
    null_cache = {}
    for mode in MODES:
        p = replace(base, **pipeline.mode_overrides(mode))
        b = batch.prepare_batch(samples, ds["reference"], ds.get("contaminants", []), p,
                                depth=ds["depth"], controls=controls)
        for name in samples:
            null_rates = None
            if p.error_model != "geneious":
                peers = sorted(c for c in controls if c != name)
                k = tuple(peers)
                if k not in null_cache:
                    null_cache[k] = null_model.fit_aligned([b["reads"][c] for c in peers], ref_seq, p,
                                                           lo, hi, calibrate=p.error_model == "qcal",
                                                           ref_name=ref_name)[0]
                null_rates = null_cache[k]
            res = pipeline.analyse(b["reads"][name], ds["reference"], p, b["stats"][name],
                                   null_rates=null_rates, depth=b["depth"])
            analysed = res["_aligned"]
            bases = window_bases(analysed, lo, hi)
            calls = {(v.pos, v.ref, v.alt): v.supporting_reads for v in res["variants"]}
            rows.append(dict(dataset=key, sample=name, control=name in controls, mode=mode,
                             reads=len(analysed), bases=bases, calls=len(calls),
                             density_x1e5=1e5 * len(calls) / bases if bases else None,
                             sites=[list(c) for c in calls]))
            print(key, mode, name, len(analysed), len(calls), flush=True)
            if mode == "recommended":       # LoFreq on the very same reads
                if name not in records:
                    records[name] = {r[0].split()[0]: r for r in read_fastq(samples[name])}
                reads = [records[name][a.qname] for a in analysed if a.qname in records[name]]
                for lname, extra in LOFREQ.items():
                    lf = lofreq_calls(reads, ds["reference"], lo, hi, preset=base.minimap2_preset,
                                      extra=extra)
                    rows.append(dict(dataset=key, sample=name, control=name in controls, mode=lname,
                                     reads=len(reads), bases=bases, calls=len(lf),
                                     density_x1e5=1e5 * len(lf) / bases if bases else None,
                                     sites=[list(c) for c in lf]))
                    print(key, lname, name, len(reads), len(lf), flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump({"dataset": key, "window": [lo, hi], "depth": ds["depth"], "rows": rows},
              open(OUT / f"benchmark_real_{key}.json", "w"), indent=1)


if __name__ == "__main__":
    for k in sys.argv[1:]:
        run(k)
