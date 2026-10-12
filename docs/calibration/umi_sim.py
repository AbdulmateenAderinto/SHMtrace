"""Does UMI consensus quantify SHM better? Simulated UMI libraries on J꜀4 and N114.

    python docs/calibration/umi_sim.py --out umi_sim [--jobs 8]

For each amplicon, each control is held out in turn: its reads become the copies of
simulated molecules (lineage_sim.build_umi_library, ~6 copies per molecule, each copy
a different real read, UMIs written with nanopore-like errors). The other two controls
become UMI libraries with nothing planted, the null. The same simulated reads are then
analysed twice: with UMI consensus (umi="flank") and without (every read on its own).
Both are scored against the molecules' true mutations.

J꜀4 reads begin at their V(D)J junction, not at the reference start, so its UMI is at
the 3' end only; N114 has UMIs at both ends.
"""
import argparse
import json
import os
import random
import statistics
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from quant_accuracy import AMPLICONS  # noqa: E402

from shmtrace import lineage_sim as L  # noqa: E402
from shmtrace import pipeline, project  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

ENDS = {"JH4": "right", "N114": "both"}
SCEN = {"JH4": ["polyclonal", "low", "null"], "N114": ["subclonal", "low", "null"]}


def _write(lib, path):
    with open(path, "w") as fh:
        for n, s, q in lib:
            fh.write(f"@{n}\n{s}\n+\n{q}\n")


def usable(path, a, p, ref, ends):
    """The control's reads that pass the pipeline's filters and reach the UMI end(s),
    so copies are not lost to filtering after molecules are assigned."""
    aligned, _ = pipeline.prepare_reads(path, a["reference"], a["contaminants"], p)
    ok = {r.qname for r in aligned
          if (ends == "left" or r.r_en >= len(ref) - 3) and (ends == "right" or r.r_st <= 3)}
    return [x for x in read_fastq(path) if x[0] in ok]


def task(args):
    amp, scen, k, seed = args
    a = AMPLICONS[amp]
    ends = ENDS[amp]
    ctl = a["controls"][k]
    peers = [c for c in a["controls"] if c != ctl]
    p = project.build_params({"region": a["region"], "mode": "recommended"})
    ref = read_fasta(a["reference"])[0][1]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)
    rng = random.Random(seed)
    tmp = tempfile.mkdtemp(prefix="umisim_")
    reads = usable(ctl, a, p, ref, ends)
    n_mol = len(reads) // 6
    truth = L.design(n_mol, ref, lo, hi, L.SCENARIOS[scen], rng, None, p.junction_length,
                     (p.junction_ref_start or 53) - 1)
    lib, carried = L.build_umi_library(reads, a["reference"], truth, p, rng, ends=ends)
    sim = os.path.join(tmp, "sim.fastq")
    _write(lib, sim)
    nulls = []
    for i, c in enumerate(peers):
        pr = usable(c, a, p, ref, ends)
        t0 = L.design(len(pr) // 6, ref, lo, hi, L.SCENARIOS["null"], rng, None,
                      p.junction_length, (p.junction_ref_start or 53) - 1)
        nl, _ = L.build_umi_library(pr, a["reference"], t0, p, rng, ends=ends)
        f = os.path.join(tmp, f"null{i}.fastq")
        _write(nl, f)
        nulls.append(f)
    out = {"amplicon": amp, "scenario": scen, "control": k}
    for mode in ("umi", "reads"):
        over = {"umi": "flank", "umi_ends": ends} if mode == "umi" else {}
        proj = {"name": "s", "region": a["region"], "mode": "recommended",
                "reference": a["reference"], "contaminants": a["contaminants"],
                "overrides": over,
                "samples": [{"name": "SIM", "path": sim, "group": "SIM"}]
                + [{"name": f"C{i}", "path": f, "group": "C", "control": True}
                   for i, f in enumerate(nulls)]}
        r = project.run_project(proj)
        row = next(s for s in r["samples"] if s["sample"] == "SIM")
        res = r["_details"]["SIM"]
        if mode == "umi":
            sc = L.score_umi(res, row, carried, lo, hi)
            sc.update({k2: v for k2, v in res["summary"]["reads"].items() if k2.startswith("umi_")})
        else:
            # truth per molecule even here: what biology is being measured
            sc = L.score_umi_reads(res, row, carried, lo, hi)
        out[mode] = sc
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="umi_sim")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tasks = [(amp, scen, k, 100 + 7 * k + 31 * j)
             for amp in ("JH4", "N114") for j, scen in enumerate(SCEN[amp])
             for k in ((0,) if a.quick else (0, 1, 2))]
    with ProcessPoolExecutor(a.jobs) as ex:
        runs = list(ex.map(task, tasks))
    json.dump(runs, open(os.path.join(a.out, "runs.json"), "w"), indent=1, default=str)
    lines = ["| Amplicon | Scenario | Analysis | Units | False calls | Recall | Recall, 1-molecule "
             "mutations | Frequency est/true | Unique density est/true | Group purity |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for amp in ("JH4", "N114"):
        for scen in SCEN[amp]:
            rs = [r for r in runs if r["amplicon"] == amp and r["scenario"] == scen]
            for mode in ("reads", "umi"):
                xs = [r[mode] for r in rs]

                def m(key, xs=xs):
                    v = [x.get(key) for x in xs if x.get(key) is not None]
                    return statistics.mean(v) if v else None

                def ratio(e, t, xs=xs):
                    v = [x[e] / x[t] for x in xs if x.get(e) is not None and x.get(t)]
                    return f"{statistics.mean(v):.2f}" if v else "–"
                f = lambda v, d=2: "–" if v is None else f"{v:.{d}f}"  # noqa: E731
                units = m("molecules_analysed") if mode == "umi" else m("reads_analysed")
                lines.append(f"| {amp} | {scen} | {'UMI consensus' if mode == 'umi' else 'every read'} | "
                             f"{f(units, 0)} | {f(m('false_calls'), 1)} | {f(m('recall'))} | "
                             f"{f(m('recall_single_molecule'))} | "
                             f"{ratio('est_read_weighted_x1e3', 'true_read_weighted_x1e3')} | "
                             f"{ratio('est_unique_density_x1e5', 'true_unique_density_x1e5')} | "
                             f"{f(m('group_purity'))} |")
    md = "\n".join(lines) + "\n"
    open(os.path.join(a.out, "report.md"), "w").write(md)
    print(md)


if __name__ == "__main__":
    main()
