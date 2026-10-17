"""The eight nanopore reads drawn in Figure S12, over the whole J_H4 amplicon.

    python docs/calibration/crm/scicanvas/s12_reads.py   -> prints the READS array for figS12_platforms.js

The per-read raster in out/intuition_jh4.json is clipped to the 529-bp window shared with Sanger, so a read drawn
from it carries no mutation past position 629 and the panel would imply SHM stops there. It does not: in the
full-length run (figures/jh4_shmtrace_full) a quarter of wild-type calls lie beyond that window. This re-runs the
wild-type library Wt3 over the full window (positions 90-1,230, as Figure 6A) and takes eight reads of four
clones, with every called mutation each read carries.

Sequencing errors are not taken from the reads: a read's own mismatches that are not called cannot be separated
from uncalled mutations. They are drawn at the rate the AID-null controls show over the same window, per base and
per quality, which is what the figure's caption states.
"""
import json
import random
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

CRM = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CRM))
from datasets import DATASETS  # noqa: E402

from shmtrace import batch, clones, null_model, pipeline  # noqa: E402
from shmtrace.config import load_preset  # noqa: E402
from shmtrace.io_utils import read_fasta  # noqa: E402

CASE, QMIN, SEED = "WT3", 20, 450
WANT = [(0, 2), (1, 2), (3, 2), (2, 1), (None, 1)]   # clone -> how many reads to draw


def main():
    ds = DATASETS["jh4"]
    p = replace(load_preset(ds["preset"]), **ds.get("overrides", {}))
    p = replace(p, **pipeline.mode_overrides("recommended"))
    p = replace(p, analysis_window=[90, 1230])       # the reliably covered amplicon, as Figure 6A
    ref_name, ref = read_fasta(ds["reference"])[0]
    lo, hi = pipeline._window(p, ref)
    samples = {n: s["fastq"] for n, s in ds["samples"].items()}
    controls = sorted(n for n, s in ds["samples"].items() if s["control"])
    b = batch.prepare_batch(samples, ds["reference"], ds.get("contaminants", []), p, depth=ds["depth"],
                            controls=set(controls))
    peers = tuple(c for c in controls if c != CASE)
    null = null_model.fit_aligned([b["reads"][c] for c in peers], ref, p, lo, hi,
                                  calibrate=p.error_model == "qcal", ref_name=ref_name)[0]
    res = pipeline.analyse(b["reads"][CASE], ds["reference"], p, b["stats"][CASE], null_rates=null,
                           depth=b["depth"])
    aligned = res["_aligned"]
    called = {(v.pos - 1, v.alt) for v in res["variants"]}
    cl = (clones.call_clones(aligned, res["variants"], p, lo, hi) or {}).get("clone_of", {})
    # germline and parental differences, which nearly every read carries, are not sequencing error: drop the
    # sites SHMtrace drops, those where half or more of a library's reads differ from the reference
    germ = set()
    for lib in b["reads"].values():
        depth, diff = Counter(), Counter()
        for a in lib:
            for rpos, base, _q in a.pairs:
                if base and lo - 1 <= rpos <= hi - 1:
                    depth[rpos] += 1
                    diff[rpos] += base != ref[rpos]
        germ |= {rpos for rpos, d in depth.items() if d and diff[rpos] / d >= 0.5}
    # the error rate the controls show over the same window, and the share of it at Q >= 20
    qb, qm = Counter(), Counter()
    for c in controls:
        for a in b["reads"][c]:
            for rpos, base, q in a.pairs:
                if base and lo - 1 <= rpos <= hi - 1 and rpos not in germ:
                    qb[int(q)] += 1
                    qm[int(q)] += base != ref[rpos]
    rate = sum(qm.values()) / sum(qb.values())
    hq = sum(v for q, v in qm.items() if q >= QMIN) / sum(qm.values())
    reads = []
    for a in aligned:
        muts = sorted(rpos + 1 for rpos, base, _q in a.pairs
                      if base and int(_q) >= QMIN and (rpos, base) in called and rpos not in germ)
        reads.append({"clone": cl.get(a.qname), "start": a.r_st + 1, "end": a.r_en, "muts": muts})
    pick = []
    for clone, n in WANT:
        rs = [r for r in reads if r["clone"] == clone and r["start"] <= 60 and r["muts"]]
        pick += (rs or [r for r in reads if r["clone"] == clone])[:n]
    span = sum(r["end"] - r["start"] + 1 for r in pick)
    want_n, want_hq = round(span * rate), round(span * rate * hq)
    for seed in range(SEED, SEED + 5000):           # a draw at the expected counts, clear of the mutations
        rng = random.Random(seed)
        errs = []
        for r in pick:
            k = sum(1 for _ in range(r["end"] - r["start"] + 1) if rng.random() < rate)
            errs.append(sorted((rng.randint(r["start"], r["end"]), int(rng.random() < hq)) for _ in range(k)))
        n, h = sum(len(e) for e in errs), sum(x[1] for e in errs for x in e)
        if n == want_n and h == want_hq and not any(abs(e[0] - m) < 12 for e, r in zip(errs, pick, strict=True)
                                                    for e in e for m in r["muts"]):
            break
    print(f"// window {lo}-{hi}; {len(aligned)} reads of {CASE}; control error {rate:.5f}/base, "
          f"{100 * hq:.0f}% at Q >= {QMIN}; {want_n} errors drawn over these reads, {want_hq} at Q >= {QMIN} "
          f"(seed {seed})")
    print("  const READS = " + json.dumps([[r["start"], r["end"], r["muts"], e]
                                           for r, e in zip(pick, errs, strict=True)], separators=(",", " ")) + ";")
    beyond = sum(1 for r in pick for m in r["muts"] if m > 629)
    print(f"// mutations drawn: {sum(len(r['muts']) for r in pick)}, of which {beyond} lie beyond position 629")


if __name__ == "__main__":
    main()
