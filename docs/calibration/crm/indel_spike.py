"""Indel recall and false calls on held-out nanopore controls (J_H4 and N114), as the CRM spike-ins.

    python docs/calibration/crm/indel_spike.py jh4 n114 [--draws 8]

For each AID-null control in turn: its usable reads are drawn to the dataset's depth (seed 0), 10 indels
of mixed kinds (1-bp deletions in and outside homopolymers, 1- and 2-bp insertions, 2- and 3-bp
deletions, 10-40-bp deletions) are planted into 1, 2, 3, 5, 10 or 20 reads each with the reads' own base
qualities kept (shmtrace.simulate.spike_indels), and the library is called against an indel null fitted
and calibrated on the other controls only, at the defaults. As docs/calibration/indel_spike.py, on the
delivered reads and the CRM datasets. Writes out/indel_spike_<dataset>.json.
"""
import argparse
import json
import random
import sys
import tempfile
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from shmtrace import batch, indels, pipeline, simulate
from shmtrace.config import load_preset
from shmtrace.io_utils import read_fasta, read_fastq
from shmtrace.variants import homopolymer_run

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402

OUT = Path(__file__).parent / "out"
SUPPORT = (1, 2, 3, 5, 10, 20)
KINDS = ["del1_hp", "del1", "ins1", "del2", "ins2", "del3", "dellong"]


def pick_sites(ref, lo, hi, rng, n=10):
    """Spaced sites of mixed kinds: (pos0, "ins"/"del", bases, reads, kind)."""
    sites, used = [], []
    for _ in range(10000):
        if len(sites) == n:
            break
        kind = rng.choice(KINDS)
        L = {"del1_hp": 1, "del1": 1, "ins1": 1, "del2": 2, "ins2": 2, "del3": 3,
             "dellong": rng.randint(10, 40)}[kind]
        pos = rng.randint(lo + 5, hi - L - 5)
        if (kind == "del1_hp") != (homopolymer_run(ref, pos) >= 3):
            continue
        if any(abs(pos - u) < 10 + max(L, ul) for u, ul in used):
            continue
        site = (pos, "ins", "".join(rng.choice("ACGT") for _ in range(L))) if kind.startswith("ins") \
            else (pos, "del", ref[pos:pos + L])
        sites.append((*site, rng.choice(SUPPORT), kind))
        used.append((pos, L))
    return sites


def run(key, draws):
    ds = DATASETS[key]
    p = replace(load_preset(ds["preset"]), **ds.get("overrides", {}), exclude_indels=False)
    ref = read_fasta(ds["reference"])[0][1]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)
    contam = ds.get("contaminants", [])
    controls = [n for n, s in ds["samples"].items() if s["control"]]
    tmp = Path(tempfile.mkdtemp())

    def pile(path, name):
        b = batch.prepare_batch({name: str(path)}, ds["reference"], contam, p, depth=ds["depth"])
        return indels.Pileup(b["reads"][name], ref, lo, hi, p.indel_min_quality)

    rows = []
    for held in controls:
        peers = [pile(ds["samples"][c]["fastq"], c) for c in controls if c != held]
        null = indels.IndelNull(peers, ref)
        null.qscale = indels.calibrate(peers, ref, p)[0]
        fq = ds["samples"][held]["fastq"]
        usable, _ = pipeline.prepare_reads(fq, ds["reference"], contam, p)
        names = {a.qname for a in usable}
        pool = [r for r in read_fastq(fq) if r[0].split()[0] in names]
        reads = random.Random(0).sample(pool, min(len(pool), ds["depth"]))
        for draw in range(draws):
            rng = random.Random(100 * draw + controls.index(held))
            sites = pick_sites(ref, lo, hi, rng)
            spiked, planted = simulate.spike_indels(reads, ds["reference"], [s[:4] for s in sites], p, rng)
            f = tmp / f"{held}_{draw}.fastq"
            with open(f, "w") as fh:
                for n, s, q in spiked:
                    fh.write(f"@{n}\n{s}\n+\n{q}\n")
            called, _ = indels.call(pile(f, "S"), null, p)
            got = {(v.pos - 1, v.kind, v.seq) for v in called}
            kind_of = {indels.normalise(ref, pos, "D" if k == "del" else "I",
                                        len(b) if k == "del" else b): label
                       for pos, k, b, _n, label in sites}
            found = [{"kind": kind_of[k], "reads": n, "found": k in got} for k, n in planted.items() if n]
            false = [{"kind": v.kind, "length": v.length, "reads": v.supporting_reads}
                     for v in called if (v.pos - 1, v.kind, v.seq) not in planted]
            rows.append({"control": held, "draw": draw, "planted": found, "false": false})
            print(key, held, draw, f"found {sum(x['found'] for x in found)}/{len(found)}, false {len(false)}",
                  flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump({"dataset": key, "draws": draws, "rows": rows}, open(OUT / f"indel_spike_{key}.json", "w"),
              indent=1)


def summary(key):
    d = json.load(open(OUT / f"indel_spike_{key}.json"))
    by, kind = defaultdict(lambda: [0, 0]), defaultdict(lambda: [0, 0])
    for r in d["rows"]:
        for x in r["planted"]:
            by[x["reads"]][0] += x["found"]
            by[x["reads"]][1] += 1
            if x["reads"] >= 3:
                kind[x["kind"]][0] += x["found"]
                kind[x["kind"]][1] += 1
    fp = sum(len(r["false"]) for r in d["rows"]) / len(d["rows"])
    tp = sum(x["found"] for r in d["rows"] for x in r["planted"])
    return {"by_reads": {k: v[0] / v[1] for k, v in sorted(by.items())},
            "by_kind_ge3": {k: v[0] / v[1] for k, v in sorted(kind.items())},
            "ge3": sum(v[0] for k, v in by.items() if k >= 3) / sum(v[1] for k, v in by.items() if k >= 3),
            "fp": fp, "fdr": sum(len(r["false"]) for r in d["rows"]) / max(1, tp + sum(len(r["false"]) for r in d["rows"]))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="+")
    ap.add_argument("--draws", type=int, default=8)
    a = ap.parse_args()
    for k in a.datasets:
        run(k, a.draws)
        print(k, json.dumps(summary(k), indent=1))
