"""In-silico spike-in validation: known mutations planted in real control reads.

There is no ground truth in a real SHM library, so sensitivity and false
discovery cannot be read off the WT samples. This module makes ground truth:

1. take reads from a mutation-free control (Aicda-/-, empty vector),
2. choose sites the way AID would (WRC/GYW hotspots weighted up, G:C
   transitions favoured) and write each chosen substitution into a known number
   of reads, keeping the sequencer's own base quality at that position,
3. run the caller exactly as on real data (alignment, phasing, test, FDR),
4. score: recall by number of supporting reads, and false calls, which are
   everything called that was not planted.

Because the background comes from real nanopore reads, the false-call side is
realistic; the planted side is idealised (a planted base never has a
quality-score mismatch or an alignment artefact of its own). Use a control that
is *not* part of the null, so the test is not graded on its own training data.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
from collections import defaultdict
from dataclasses import replace

from . import motifs, null_model, pipeline
from .align import build_index
from .config import Params, load_preset
from .io_utils import read_fasta, read_fastq

BASES = ("A", "C", "G", "T")
_COMP = str.maketrans("ACGT", "TGCA")
DEFAULT_READ_COUNTS = (1, 2, 3, 5, 10, 20)


def _weight(ref: str, pos0: int) -> float:
    if motifs.is_site("WRC", ref, pos0) or motifs.is_site("GYW", ref, pos0):
        return 4.0
    if motifs.is_site("WA", ref, pos0) or motifs.is_site("TW", ref, pos0):
        return 2.0
    return 1.0


def _alt_for(refbase: str, rng: random.Random) -> str:
    """G:C transitions half the time at G/C (the AID replication signature)."""
    if refbase == "C" and rng.random() < 0.5:
        return "T"
    if refbase == "G" and rng.random() < 0.5:
        return "A"
    return rng.choice([b for b in BASES if b != refbase])


def choose_sites(ref: str, lo: int, hi: int, n: int, rng: random.Random,
                 read_counts=DEFAULT_READ_COUNTS) -> list[tuple[int, str, int]]:
    """Pick n distinct positions in [lo, hi] (0-based), AID-weighted.

    Returns (pos0, alt, reads_to_carry) with read counts cycled through
    `read_counts` so every support level is equally represented.
    """
    pool = [i for i in range(lo, hi + 1) if ref[i] in BASES]
    weights = [_weight(ref, i) for i in pool]
    chosen: list[int] = []
    while len(chosen) < min(n, len(pool)):
        i = rng.choices(range(len(pool)), weights=weights)[0]
        if weights[i] > 0:
            chosen.append(pool[i])
            weights[i] = 0.0
    chosen.sort()
    return [(pos, _alt_for(ref[pos], rng), read_counts[k % len(read_counts)])
            for k, pos in enumerate(chosen)]


def _ref_to_query(hit, qlen: int) -> dict[int, int]:
    """Map reference positions to positions on the read as sequenced."""
    out = {}
    # mappy's CIGAR walks the reference forward; for a reverse hit that is the
    # read's reverse complement, which starts qlen - q_en into it
    q = hit.q_st if hit.strand == 1 else qlen - hit.q_en
    r = hit.r_st
    for length, op in hit.cigar:
        if op in (0, 7, 8):
            for k in range(length):
                out[r + k] = q + k
            q += length
            r += length
        elif op in (1, 4):
            q += length
        elif op in (2, 3):
            r += length
    return out


def spike(reads: list[tuple[str, str, str]], reference: str, sites, p: Params,
          rng: random.Random) -> tuple[list[tuple[str, str, str]], dict]:
    """Write each planted substitution into its number of reads.

    A read can only carry a site it covers with an aligned (not deleted) base.
    Reverse-strand reads get the complement at the mirrored position, so every
    read stays in its sequenced orientation. Returns the new reads and how many
    reads actually received each site.
    """
    aligner, _t, _c = build_index([reference], [], p.minimap2_preset)
    seqs = [list(s) for _n, s, _q in reads]
    covering: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for idx, (_name, seq, _qual) in enumerate(reads):
        hit = next((h for h in aligner.map(seq) if h.is_primary), None)
        if hit is None:
            continue
        L = len(seq)
        r2q = _ref_to_query(hit, L)
        for pos, _alt, _k in sites:
            if pos in r2q:
                qpos = r2q[pos]
                covering[pos].append((idx, qpos if hit.strand == 1 else L - 1 - qpos,
                                      hit.strand))
    planted = {}
    for pos, alt, k in sites:
        cands = covering.get(pos, [])
        take = rng.sample(cands, min(k, len(cands)))
        for idx, qpos, strand in take:
            seqs[idx][qpos] = alt if strand == 1 else alt.translate(_COMP)
        planted[(pos, alt)] = len(take)
    out = [(n, "".join(seqs[i]), q) for i, (n, _s, q) in enumerate(reads)]
    return out, planted


def evaluate(called: set, planted: dict) -> dict:
    """Recall per planted read count, and false calls."""
    by_k: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for (pos, alt), k in planted.items():
        if k == 0:
            continue
        by_k[k][1] += 1
        by_k[k][0] += (pos, alt) in called
    truth = {key for key, k in planted.items() if k > 0}
    tp = len(called & truth)
    fp = len(called - truth)
    return {"recall_by_reads": {k: (v[0], v[1]) for k, v in sorted(by_k.items())},
            "true_calls": tp, "false_calls": fp, "planted": len(truth),
            "recall": tp / len(truth) if truth else None,
            "fdr": fp / (tp + fp) if (tp + fp) else 0.0}


PROFILES = ("geneious", "recommended", "v0.1", "qcal", "qcal+filter", "strand",
            "strand+filter", "pbinom", "recommended+filter")


def run_validation(control: str, nulls: list[str], reference: str, contaminants: list[str],
                   p: Params, profiles: list[str], n_sites: int = 60, reps: int = 3,
                   seed: int = 1, depth: int | None = None) -> dict:
    ref = read_fasta(reference)[0][1]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)
    reads = list(read_fastq(control))
    null = None
    if nulls:
        null, _s = null_model.fit(nulls, reference, contaminants, p)
    rng = random.Random(seed)
    depth = depth or p.downsample_to
    rows = []
    for rep in range(reps):
        sub = rng.sample(reads, min(len(reads), int((depth or len(reads)) * 1.3)))
        sites = choose_sites(ref, lo, hi, n_sites, rng)
        spiked, planted = spike(sub, reference, sites, p, rng)
        tmp = f"/tmp/nanoshm_spike_{os.getpid()}_{rep}.fastq"
        with open(tmp, "w") as fh:
            for name, seq, qual in spiked:
                fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
        try:
            for prof in profiles:
                pp = replace(p, **pipeline.mode_overrides(prof))
                nr = None if pp.error_model == "geneious" else null
                # all reads go in: planted reads must not be lost to downsampling
                res = pipeline.run(tmp, reference, contaminants,
                                   replace(pp, downsample_to=None), nr)
                called = {(v.pos - 1, v.alt) for v in res["variants"]}
                ev = evaluate(called, planted)
                rows.append({"rep": rep, "profile": prof,
                             "reads": len(res["_aligned"]), **ev})
        finally:
            os.remove(tmp)
    return {"control": os.path.basename(control),
            "nulls": [os.path.basename(f) for f in nulls],
            "sites_per_rep": n_sites, "reps": reps, "rows": rows,
            "summary": summarise(rows)}


def summarise(rows: list[dict]) -> dict:
    out = {}
    for prof in dict.fromkeys(r["profile"] for r in rows):
        rs = [r for r in rows if r["profile"] == prof]
        byk: dict[int, list[int]] = defaultdict(lambda: [0, 0])
        for r in rs:
            for k, (hit, n) in r["recall_by_reads"].items():
                byk[k][0] += hit
                byk[k][1] += n
        tp = sum(r["true_calls"] for r in rs)
        fp = sum(r["false_calls"] for r in rs)
        out[prof] = {"recall_by_reads": {k: round(h / n, 3) if n else None
                                         for k, (h, n) in sorted(byk.items())},
                     "recall": round(tp / sum(r["planted"] for r in rs), 3),
                     "false_calls_per_sample": round(fp / len(rs), 1),
                     "fdr": round(fp / (tp + fp), 3) if tp + fp else 0.0}
    return out


def report_markdown(res: dict) -> str:
    ks = sorted({k for s in res["summary"].values() for k in s["recall_by_reads"]})
    out = ["# Spike-in validation", "",
           f"Control reads: {res['control']}; null fitted from: "
           f"{', '.join(res['nulls']) or 'none'}. {res['sites_per_rep']} planted sites "
           f"x {res['reps']} replicate draws.", "",
           "Recall by number of reads carrying the planted mutation, overall recall, "
           "false calls (anything called that was not planted) per sample, and FDR.", "",
           "| Profile | " + " | ".join(f"{k} read{'s' if k > 1 else ''}" for k in ks)
           + " | Recall | False calls | FDR |",
           "|---|" + "---:|" * (len(ks) + 3)]
    for prof, s in res["summary"].items():
        cells = [f"{s['recall_by_reads'].get(k, 0):.2f}" for k in ks]
        out.append(f"| {prof} | " + " | ".join(cells)
                   + f" | {s['recall']:.2f} | {s['false_calls_per_sample']} | {s['fdr']:.2f} |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description="Spike-in validation of nanoshm calling profiles")
    ap.add_argument("--control", required=True, help="mutation-free reads to plant mutations in")
    ap.add_argument("--null", action="append", default=[],
                    help="other mutation-free controls for the empirical null (repeatable)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--contaminant", action="append", default=[])
    ap.add_argument("--preset", required=True)
    ap.add_argument("--profiles", default="geneious,recommended,strand,strand+filter",
                    help=f"comma-separated, from {', '.join(PROFILES)}")
    ap.add_argument("--sites", type=int, default=60)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--depth", type=int)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="validation")
    args = ap.parse_args()
    p = load_preset(args.preset)
    res = run_validation(args.control, args.null, args.reference, args.contaminant, p,
                         args.profiles.split(","), args.sites, args.reps, args.seed, args.depth)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "validation.json"), "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    with open(os.path.join(args.out, "validation.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rep", "profile", "reads", "planted", "true_calls", "false_calls",
                    "recall", "fdr"])
        for r in res["rows"]:
            w.writerow([r["rep"], r["profile"], r["reads"], r["planted"], r["true_calls"],
                        r["false_calls"], r["recall"], r["fdr"]])
    md = report_markdown(res)
    with open(os.path.join(args.out, "report.md"), "w") as fh:
        fh.write(md)
    print(md)


if __name__ == "__main__":
    main()
