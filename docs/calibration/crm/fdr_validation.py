"""Validate the whole calling procedure on fresh simulations: threshold, margin and the noisy-control rule.

    python crm/fdr_validation.py [cells]        -> crm/out/fdr_validation.json

Issues 3 and 4 of the review ask for the same thing from two directions. Issue 3: a threshold chosen so that
the estimated false calls are at most 5% of the discoveries is not a finite-sample guarantee, the decoys come
from a few control libraries, and the same controls help choose the safety factor, so the *whole* procedure
needs testing on fresh data, reporting the per-library false-discovery proportion and its spread and the
chance of any discovery at all in an all-null library, not only the pooled false/(true+false). Issue 4: the
noisy-control rule must be shown with all controls against the exclusion, side by side, and tested under
heterogeneous nulls.

Both are answered by running `project.run_project` itself -- the full procedure, including threshold
selection and the noisy-decoy rule -- over fresh libraries built from real AID-null reads, and varying the
things the review names: how many controls there are, how heterogeneous they are, the depth, the mutation
burden, and whether an artifact is shared with the controls or absent from them.

What each library is
  case        built from KO1's reads with a known set of planted mutations (none, for an all-null library)
  controls    built from KO2 and KO3, each split in two, so up to four controls exist that share no molecule
              with the case; the two halves of one library share its artifacts, as two controls from one run do

The noise conditions (`NOISE`) are the review's list:
  none              the controls as delivered
  one_noisy         one control carries extra errors at recurrent positions: a heterogeneous null
  case_only         the case carries recurrent artifacts at positions no control has: the failure mode the
                    control-based threshold cannot see, and the one worth measuring honestly
  shared            the case and every control carry the same recurrent artifacts: what the controls exist for

Every cell is run twice, with the noisy-control rule on and off (`exclude_noisy_decoys`), on the identical
libraries, which is the side-by-side issue 4 asks for.

Reported per run: discoveries, true and false calls, the false-discovery proportion of that library (every
call is false in an all-null library, so its FDP is 1 whenever it discovers anything), the threshold the
procedure chose, which control it set aside, and the decoy noise it measured.
"""
import json
import os
import random
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from datasets import DATASETS  # noqa: E402

from shmtrace import lineage_sim, pipeline, project, simulate  # noqa: E402
from shmtrace.align import build_index  # noqa: E402
from shmtrace.lineage_sim import _COMP  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

CRM = Path(__file__).resolve().parent
OUT = CRM / "out"
REF = Path(__file__).resolve().parents[3] / "src/shmtrace/references"
REGION = "jh4_intron"
REFERENCE = str(REF / "mouse_jh4_intron.fasta")
CONTAM = [str(REF / "mouse_jh4_contaminant.fasta")]

CASE_SRC = "KO1"                      # the case library's donor reads
CTL_SRC = ("KO2", "KO3")              # split in two each: four controls that share no molecule with the case
NOISE = ("none", "one_noisy", "case_only", "shared")
ARTIFACT_SITES = 6                    # recurrent positions an artifact occupies
ARTIFACT_SHARE = 0.06                 # share of reads carrying it: above the error floor, below the 50% mask


def _sources():
    s = DATASETS["jh4"]["samples"]
    return {k: v["fastq"] for k, v in s.items()}


def split_controls(rng):
    """KO2 and KO3, each split into two halves: four control read sets, none sharing a molecule with KO1."""
    src = _sources()
    out = {}
    for name in CTL_SRC:
        reads = list(read_fastq(src[name]))
        rng.shuffle(reads)
        half = len(reads) // 2
        out[f"{name}a"], out[f"{name}b"] = reads[:half], reads[half:]
    return out


def add_artifact(reads, ref, aligner, positions, share, rng):
    """Write a recurrent substitution into `share` of the reads at each of `positions`.

    The positions are *reference* coordinates, so each read is aligned and the base is written at the query
    index that reference position maps to, exactly as lineage_sim plants a mutation. Writing at a read index
    instead would smear the artifact across the amplicon, because nanopore reads do not start together. The
    changed base keeps the read's own quality, as a real recurrent artifact does: it passes the quality
    filter, which is what makes it dangerous.
    """
    alt_of = {pos: rng.choice([b for b in "ACGT" if b != ref[pos]]) for pos in positions}
    out = []
    for name, seq, qual in reads:
        hit = next((h for h in aligner.map(seq) if h.is_primary), None)
        if hit is None:
            out.append((name, seq, qual))
            continue
        r2q = simulate._ref_to_query(hit, len(seq))
        fs = list(seq if hit.strand == 1 else seq.translate(_COMP)[::-1])
        fq = qual if hit.strand == 1 else qual[::-1]
        for pos in positions:
            if pos in r2q and rng.random() <= share:
                fs[r2q[pos]] = alt_of[pos]
        out.append((name, "".join(fs), fq))
    return out


def write_fastq(path, reads):
    with open(path, "w") as fh:
        for name, seq, qual in reads:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def run_cell(signal, n_controls, noise, exclude, depth, seed, workdir, estimate="mean"):
    """One simulated experiment, run end to end by the real pipeline."""
    rng = random.Random(seed)
    p = project.build_params({"region": REGION, "mode": "recommended"})
    ref = read_fasta(REFERENCE)[0][1]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)

    # the case library: KO1's reads with a known set of planted mutations ("null" plants none)
    case_reads = list(read_fastq(_sources()[CASE_SRC]))
    truth = lineage_sim.design(len(case_reads), ref, lo, hi, lineage_sim.SCENARIOS[signal], rng,
                               None, p.junction_length, (p.junction_ref_start or 53) - 1)
    lib, carried = lineage_sim.build_library(case_reads, REFERENCE, truth, p, rng)

    ctl = split_controls(rng)
    names = sorted(ctl)[:n_controls]
    sites = rng.sample(range(lo, hi + 1), ARTIFACT_SITES)
    aligner, _t, _c = build_index([REFERENCE], [], p.minimap2_preset)
    if noise == "case_only":
        lib = add_artifact(lib, ref, aligner, sites, ARTIFACT_SHARE, rng)
    elif noise == "shared":
        lib = add_artifact(lib, ref, aligner, sites, ARTIFACT_SHARE, rng)
        ctl = {k: (add_artifact(v, ref, aligner, sites, ARTIFACT_SHARE, rng) if k in names else v)
               for k, v in ctl.items()}
    elif noise == "one_noisy":
        first = names[0]
        ctl = {k: (add_artifact(v, ref, aligner, sites, ARTIFACT_SHARE, rng) if k == first else v)
               for k, v in ctl.items()}

    paths = {}
    case_path = os.path.join(workdir, f"case_{seed}.fastq")
    write_fastq(case_path, lib)
    for n in names:
        paths[n] = os.path.join(workdir, f"{n}_{seed}.fastq")
        write_fastq(paths[n], ctl[n])
    proj = {"name": "fdr", "region": REGION, "reference": REFERENCE, "contaminants": CONTAM,
            "mode": "recommended", "depth": depth,
            "overrides": {"exclude_noisy_decoys": exclude, "control_decoy_estimate": estimate},
            "samples": [{"name": "CASE", "path": case_path, "group": "CASE"}]
                       + [{"name": n, "path": paths[n], "group": "CTRL", "control": True} for n in names]}
    try:
        out = project.run_project(proj)
    finally:
        for f in [case_path, *paths.values()]:
            if os.path.exists(f):
                os.remove(f)

    row = next(s for s in out["samples"] if s["sample"] == "CASE")
    res = out["_details"]["CASE"]
    analysed = [r.qname for r in res["_aligned"]]
    support = Counter(m for n in analysed for m in carried.get(n, {}).get("muts", ()) if lo <= m[0] <= hi)
    truth_set = set(support)
    called = {(v.pos - 1, v.alt) for v in res["variants"]}
    tp, fp = len(called & truth_set), len(called - truth_set)
    cal = out.get("calibration") or {}
    s_call = row.get("summary", {}).get("calling", {}) if isinstance(row.get("summary"), dict) else {}
    return {
        "signal": signal, "controls": n_controls, "noise": noise, "exclude_noisy": exclude,
        "depth": depth, "seed": seed, "control_names": names, "estimate": estimate,
        "planted": len(truth_set), "discoveries": len(called), "true_calls": tp, "false_calls": fp,
        # the false-discovery proportion of this one library: 1 whenever an all-null library discovers anything
        "fdp": round(fp / len(called), 4) if called else 0.0,
        "any_discovery": bool(called),
        "threshold": s_call.get("control_fdr_threshold") or res.get("summary", {}).get(
            "calling", {}).get("control_fdr_threshold"),
        "decoys_set_aside": cal.get("decoys_set_aside") or [],
        "decoy_noise_x1e5": cal.get("decoy_noise_x1e5") or {},
        "uneven_control_pair": cal.get("uneven_control_pair"),
    }


# (signal, controls, noise, depth): the grid, kept to what the review asks about
def grid(reps=5):
    cells = []
    for signal in ("null", "polyclonal"):
        for noise in NOISE:
            for exclude in (True, False):
                for rep in range(reps):
                    cells.append((signal, 3, noise, exclude, 500, 7000 + rep))
    # the rule needs a majority to compare against: with two controls it cannot fire
    for exclude in (True, False):
        for rep in range(reps):
            cells.append(("null", 2, "one_noisy", exclude, 500, 7000 + rep))
    # depth
    for depth in (250, 1000):
        for exclude in (True, False):
            for rep in range(reps):
                cells.append(("null", 3, "one_noisy", exclude, depth, 7000 + rep))
    return cells


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cells = grid(reps)
    print(f"{len(cells)} runs", flush=True)
    rows = []
    workdir = tempfile.mkdtemp(prefix="shmtrace_fdr_")
    for i, (signal, nc, noise, exclude, depth, seed) in enumerate(cells, 1):
        r = run_cell(signal, nc, noise, exclude, depth, seed, workdir)
        rows.append(r)
        print(f"[{i}/{len(cells)}] {signal} ctl={nc} {noise} excl={exclude} d={depth} "
              f"-> {r['discoveries']} calls, {r['false_calls']} false, FDP {r['fdp']}, "
              f"set aside {r['decoys_set_aside']}", flush=True)
        OUT.mkdir(exist_ok=True)
        json.dump({"runs": rows}, open(OUT / "fdr_validation.json", "w"), indent=1)
    print("written", OUT / "fdr_validation.json")


if __name__ == "__main__":
    main()
