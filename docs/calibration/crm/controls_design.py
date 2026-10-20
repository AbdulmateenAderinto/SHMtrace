"""How many AID-null controls, and which: recall and false calls as the control set changes.

    python crm/controls_design.py jh4 [reps]     -> crm/out/controls_design_jh4.json
    python crm/controls_design.py su [reps]      -> crm/out/controls_design_su.json

jh4  The whole procedure (`project.run_project`: margin, decoys, threshold, noisy-control rule) on fresh libraries
     built from real J_H4 AID-null reads, as in fdr_validation.py: the case is KO1 carrying planted polyclonal
     lineages (or nothing, for an all-null library); the controls are halves of KO2 and KO3, so up to four exist
     and none shares a molecule with the case. Varied: the number of controls (1-4), and one control carrying
     recurrent artifacts, with the noisy-control rule on and off. Halves of one library are technical partitions
     of one animal, not new animals, and share that library's artifacts.

su   The same locus from two experiments of PRJNA1481264: Sμ in LPS-stimulated B cells and in CH12 cells, each with
     three AID-null libraries. Each LPS AID-null library in turn receives 60 planted substitutions (1-20 reads
     each, three draws) and is called against: its own experiment's other two controls; two or three CH12
     controls; or all five pooled. The SRA reads carry no run identifiers, so this compares a matched with an
     unmatched experiment, not verified sequencing runs.
"""
import json
import random
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from datasets import DATASETS  # noqa: E402

from shmtrace import null_model, pipeline, simulate  # noqa: E402
from shmtrace.config import load_preset  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
WORKERS = 6


# ---------------- J_H4: number of controls and one noisy control, through run_project
def _jh4_cell(cell):
    import fdr_validation as fv
    signal, n, noise, exclude, seed = cell
    workdir = tempfile.mkdtemp(prefix="shmtrace_ctl_")
    return fv.run_cell(signal, n, noise, exclude, 500, seed, workdir)


def jh4(reps):
    cells = [(signal, n, "none", True, 8000 + rep)
             for signal in ("polyclonal", "null") for n in (1, 2, 3, 4) for rep in range(reps)]
    cells += [(signal, n, "one_noisy", exclude, 8000 + rep)
              for signal in ("polyclonal", "null") for n in (3, 4) for exclude in (True, False) for rep in range(reps)]
    print(f"{len(cells)} runs", flush=True)
    rows = []
    with ProcessPoolExecutor(WORKERS) as ex:
        for i, r in enumerate(ex.map(_jh4_cell, cells), 1):
            rows.append(r)
            print(f"[{i}/{len(cells)}] {r['signal']} ctl={r['controls']} {r['noise']} excl={r['exclude_noisy']} -> "
                  f"{r['true_calls']}/{r['planted']} true, {r['false_calls']} false", flush=True)
    json.dump({"runs": rows}, open(OUT / "controls_design_jh4.json", "w"), indent=1)


# ---------------- Sμ: controls from the matched experiment against another experiment's
ARMS = {"matched (2)": ("lps", 2), "other experiment (2)": ("ch12", 2), "other experiment (3)": ("ch12", 3),
        "pooled (5)": ("both", 5)}


def _su_case(job):
    case, reps, n_sites, seed = job
    lps, ch12 = DATASETS["su_lps"], DATASETS["su_ch12"]
    base = replace(load_preset(lps["preset"]), **lps.get("overrides", {}))
    pp = replace(base, **pipeline.mode_overrides("recommended"))
    ref_path = lps["reference"]
    ref = read_fasta(ref_path)[0][1]
    lo, hi = pipeline._window(base, ref)
    lo0, hi0 = max(lo, base.end_trim_bp), min(hi, len(ref) - 1 - base.end_trim_bp)
    lps_ctl = [s["fastq"] for k, s in lps["samples"].items() if s["control"] and k != case]
    ch12_ctl = [s["fastq"] for s in ch12["samples"].values() if s["control"]]
    peers = {"lps": lps_ctl, "ch12": ch12_ctl, "both": lps_ctl + ch12_ctl}
    fq = lps["samples"][case]["fastq"]
    usable, _ = pipeline.prepare_reads(fq, ref_path, [], base)
    names = {a.qname for a in usable}
    pool = [r for r in read_fastq(fq) if r[0].split()[0] in names]
    nulls = {arm: null_model.fit(peers[src][:n], ref_path, [], pp)[0] for arm, (src, n) in ARMS.items()}
    rng, rows = random.Random(seed), []
    for rep in range(reps):
        sub = rng.sample(pool, min(len(pool), lps["depth"]))
        sites = simulate.choose_sites(ref, lo0, hi0, n_sites, rng)
        spiked, planted = simulate.spike(sub, ref_path, sites, base, rng)
        with tempfile.NamedTemporaryFile("w", suffix=".fastq", delete=False) as fh:
            for n, s, q in spiked:
                fh.write(f"@{n}\n{s}\n+\n{q}\n")
            tmp = fh.name
        try:
            for arm in ARMS:
                res = pipeline.run(tmp, ref_path, [], replace(pp, downsample_to=None), nulls[arm])
                called = {(v.pos - 1, v.alt) for v in res["variants"]}
                rows.append({"case": case, "rep": rep, "arm": arm, **simulate.evaluate(called, planted)})
                print(case, rep, arm, rows[-1]["recall"], rows[-1]["false_calls"], flush=True)
        finally:
            Path(tmp).unlink()
    return rows


def su(reps, n_sites=60):
    cases = [k for k, s in DATASETS["su_lps"]["samples"].items() if s["control"]]
    rows = []
    with ProcessPoolExecutor(len(cases)) as ex:
        for part in ex.map(_su_case, [(c, reps, n_sites, 1) for c in cases]):
            rows += part
    json.dump({"arms": list(ARMS), "reps": reps, "sites": n_sites, "rows": rows},
              open(OUT / "controls_design_su.json", "w"), indent=1)


if __name__ == "__main__":
    part, reps = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else (10 if sys.argv[1] == "jh4" else 3)
    OUT.mkdir(exist_ok=True)
    {"jh4": jh4, "su": su}[part](reps)
