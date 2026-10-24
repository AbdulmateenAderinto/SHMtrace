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

    python docs/calibration/umi_sim.py --grouping --out umi_sim

tests umi.group alone, at library sizes the amplicon controls cannot reach: N molecules,
each with 1 + Poisson(5) copies of an 18-nt UMI (one end) or two (both ends), read at
Q15 with 2% indels, 5% of copies cut 1-2 bases short at the read end. UMIs are uniform,
or T-rich (35% T at every position, as the Takara SMARTer UMIs VDJtrace measured). For
the directional step at +2 and +1 edits and off: how many groups of 3 or more reads mix
molecules, and how often one molecule is split over several such groups.
"""
import argparse
import json
import os
import random
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from quant_accuracy import AMPLICONS  # noqa: E402

from shmtrace import lineage_sim as L  # noqa: E402
from shmtrace import pipeline, project, umi  # noqa: E402
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
    amp, scen, k, seed, dirs, controls = args
    a = dict(AMPLICONS[amp], controls=controls) if controls else AMPLICONS[amp]
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
    for mode in ["reads"] + [f"umi{d}" for d in dirs]:
        over = ({"umi": "flank", "umi_ends": ends, "umi_directional": int(mode[3:])}
                if mode != "reads" else {})
        proj = {"name": "s", "region": a["region"], "mode": "recommended",
                "reference": a["reference"], "contaminants": a["contaminants"],
                "overrides": over,
                "samples": [{"name": "SIM", "path": sim, "group": "SIM"}]
                + [{"name": f"C{i}", "path": f, "group": "C", "control": True}
                   for i, f in enumerate(nulls)]}
        r = project.run_project(proj)
        row = next(s for s in r["samples"] if s["sample"] == "SIM")
        res = r["_details"]["SIM"]
        if mode != "reads":
            sc = L.score_umi(res, row, carried, lo, hi)
            mem = [Counter(carried[n]["molecule"] for n in getattr(r, "members", [])
                           if n in carried).most_common(1) for r in res["_aligned"]]
            mols = Counter(m[0][0] for m in mem if m)
            sc["molecules_split_pct"] = (round(100 * sum(v > 1 for v in mols.values())
                                               / len(mols), 2) if mols else None)
            sc.update({k2: v for k2, v in res["summary"]["reads"].items() if k2.startswith("umi_")})
        else:
            # truth per molecule even here: what biology is being measured
            sc = L.score_umi_reads(res, row, carried, lo, hi)
        out[mode] = sc
    return out


def grouping_task(args):
    n_mol, ends, comp, d, seed = args
    rng = random.Random(seed)
    w = [0.217, 0.217, 0.217, 0.35] if comp == "T-rich" else [0.25] * 4

    def one():
        return "".join(rng.choices("ACGT", w, k=18))

    keys, mol = {}, {}
    for m in range(n_mol):
        umis = [one() for _ in range(2 if ends == "both" else 1)]
        for c in range(1 + L._poisson(5, rng)):
            parts = []
            for j, u in enumerate(umis):
                x, _q = L._noisy(u, "0" * 18, rng, 0.02)        # "0" is Q15
                if rng.random() < 0.05:                          # cut short at the read end
                    cut = rng.choice((1, 2))
                    x = x[cut:] if j == 0 else x[:-cut]
                parts.append(x)
            keys[f"m{m}c{c}"] = "".join(parts)
            mol[f"m{m}c{c}"] = m
    gid = umi.group(keys, 3 * len(umis), d)
    groups = defaultdict(list)
    for name, g in gid.items():
        groups[g].append(mol[name])
    big = [ms for ms in groups.values() if len(ms) >= 3]
    major = [Counter(ms).most_common(1)[0] for ms in big]
    per_mol = Counter(m for m, _k in major)
    true3 = sum(n >= 3 for n in Counter(mol.values()).values())
    return {"molecules": n_mol, "ends": ends, "umi": comp, "directional": d,
            "molecules_3plus_copies": true3, "groups_3plus_reads": len(big),
            "groups_mixing_molecules_pct": 100 * sum(k < len(ms) for ms, (_m, k)
                                                     in zip(big, major, strict=True)) / len(big),
            "group_purity": statistics.mean(k / len(ms) for ms, (_m, k)
                                            in zip(big, major, strict=True)),
            "molecules_split_pct": 100 * sum(n > 1 for n in per_mol.values()) / len(per_mol),
            "reads_in_groups_3plus_pct": 100 * sum(map(len, big)) / len(keys)}


def grouping(a):
    tasks = [(n, ends, comp, d, 1000 + n + 7 * k)
             for n in ((150, 1000) if a.quick else (150, 1000, 3000))
             for ends in ("right", "both") for comp in ("uniform", "T-rich")
             for d in (2, 1, 0) for k in range(1 if a.quick else 3)]
    with ProcessPoolExecutor(a.jobs) as ex:
        runs = list(ex.map(grouping_task, tasks))
    json.dump(runs, open(os.path.join(a.out, "grouping.json"), "w"), indent=1)
    lines = ["| Molecules | UMI | Composition | Directional | Groups ≥ 3 reads | Mixing molecules "
             "| Purity | Molecules split | Reads in groups ≥ 3 |",
             "|---:|---|---|---|---:|---:|---:|---:|---:|"]
    for n in sorted({r["molecules"] for r in runs}):
        for ends in ("right", "both"):
            for comp in ("uniform", "T-rich"):
                for d in (2, 1, 0):
                    rs = [r for r in runs if (r["molecules"], r["ends"], r["umi"],
                                              r["directional"]) == (n, ends, comp, d)]

                    def m(k, rs=rs):
                        return statistics.mean(r[k] for r in rs)
                    lines.append(
                        f"| {n} | {'18 nt' if ends == 'right' else '2 × 18 nt'} | {comp} | "
                        f"{'+' + str(d) if d else 'off'} | {m('groups_3plus_reads'):.0f} | "
                        f"{m('groups_mixing_molecules_pct'):.1f}% | {m('group_purity'):.4f} | "
                        f"{m('molecules_split_pct'):.1f}% | {m('reads_in_groups_3plus_pct'):.1f}% |")
    md = "\n".join(lines) + "\n"
    open(os.path.join(a.out, "grouping.md"), "w").write(md)
    print(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="umi_sim")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--grouping", action="store_true", help="umi.group alone, at scale")
    ap.add_argument("--directional", default="0",
                    help="umi_directional settings to compare, e.g. 2,1,0")
    ap.add_argument("--amplicons", default="JH4,N114")
    ap.add_argument("--jh4-controls", nargs=3,
                    help="three J꜀4 control FASTQs in place of the lab's (e.g. examples_compare/KO*)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.grouping:
        grouping(a)
        return
    dirs = [int(x) for x in a.directional.split(",")]
    amps = a.amplicons.split(",")
    tasks = [(amp, scen, k, 100 + 7 * k + 31 * j, dirs,
              a.jh4_controls if amp == "JH4" else None)
             for amp in amps for j, scen in enumerate(SCEN[amp])
             for k in ((0,) if a.quick else (0, 1, 2))]
    with ProcessPoolExecutor(a.jobs) as ex:
        runs = list(ex.map(task, tasks))
    json.dump(runs, open(os.path.join(a.out, "runs.json"), "w"), indent=1, default=str)
    lines = ["| Amplicon | Scenario | Analysis | Units | False calls | Recall | Recall, 1-molecule "
             "mutations | Frequency est/true | Unique density est/true | Group purity | "
             "Molecules split |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for amp in amps:
        for scen in SCEN[amp]:
            rs = [r for r in runs if r["amplicon"] == amp and r["scenario"] == scen]
            for mode in ["reads"] + [f"umi{d}" for d in dirs]:
                xs = [r[mode] for r in rs]

                def m(key, xs=xs):
                    v = [x.get(key) for x in xs if x.get(key) is not None]
                    return statistics.mean(v) if v else None

                def ratio(e, t, xs=xs):
                    v = [x[e] / x[t] for x in xs if x.get(e) is not None and x.get(t)]
                    return f"{statistics.mean(v):.2f}" if v else "–"
                f = lambda v, d=2: "–" if v is None else f"{v:.{d}f}"  # noqa: E731
                units = m("reads_analysed") if mode == "reads" else m("molecules_analysed")
                d = mode[3:]
                label = ("every read" if mode == "reads" else "UMI consensus" if len(dirs) == 1
                         else f"UMI, directional {'+' + d if d != '0' else 'off'}")
                lines.append(f"| {amp} | {scen} | {label} | "
                             f"{f(units, 0)} | {f(m('false_calls'), 1)} | {f(m('recall'))} | "
                             f"{f(m('recall_single_molecule'))} | "
                             f"{ratio('est_read_weighted_x1e3', 'true_read_weighted_x1e3')} | "
                             f"{ratio('est_unique_density_x1e5', 'true_unique_density_x1e5')} | "
                             f"{f(m('group_purity'))} | {f(m('molecules_split_pct'), 1)} |")
    md = "\n".join(lines) + "\n"
    open(os.path.join(a.out, "report.md"), "w").write(md)
    print(md)


if __name__ == "__main__":
    main()
