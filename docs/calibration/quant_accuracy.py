"""How accurately does nanoshm quantify SHM? Truth-known libraries on both amplicons.

    python docs/calibration/quant_accuracy.py --out quant_accuracy [--jobs 8] [--quick]

Libraries are built by lineage_sim from the deposit's own mutation-free controls (J꜀4
Aicda-/- KO1-3, 5,000 reads each; N114 empty-vector C3, C4, D1, 1,240-1,865 reads),
each control held out in turn with the other two as the null, and analysed at 125, 250,
500 and 1,000 reads. For every library the reported numbers are compared with the truth:

  mutation frequency  mutations per 10^3 bp, read-weighted
  unique density      distinct mutations per base sequenced x 10^5 (the manuscript's measure)
  combined density    unique density with the single-read estimate added (SINGLE_READ.md)

Writes runs.json and report.md (est/true ratios by amplicon, scenario, depth and mode).
"""
import argparse
import json
import os
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from nanoshm import lineage_sim

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "src/nanoshm/references"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data")
AMPLICONS = {
    "JH4": {"region": "jh4_intron", "reference": str(REF / "mouse_jh4_intron.fasta"),
            "contaminants": [str(REF / "mouse_jh4_contaminant.fasta")],
            "controls": [str(DEP / f"02_JH4_ONT/Original_FASTQ/DHJGKB_{i}_AID_KO_{i - 4}.fastq")
                         for i in (5, 6, 7)],
            "scenarios": ["polyclonal", "expanded", "low", "null"]},
    "N114": {"region": "n114_vregion", "reference": str(REF / "n114_vregion.fasta"),
             "contaminants": [],
             "controls": [str(DEP / f"03_N114_ONT/Processed_FASTQ_Depth_Subsets/{n}.fastq")
                          for n in ("75K3B3_1_C3_Empty_March_2024_1_N114",
                                    "75K3B3_3_C4_Empty_March_2024_N114",
                                    "75K3B3_4_D1_Empty_March_2024_N114")],
             "scenarios": ["subclonal", "low", "null"]},
}
DEPTHS = (125, 250, 500, 1000)
MODES = ("recommended", "geneious")
KEYS = ("mutation_frequency", "unique_density", "combined_density")


def task(args):
    amp, scen, depth, mode, k, seed = args
    a = AMPLICONS[amp]
    ctl = a["controls"][k]
    peers = [c for c in a["controls"] if c != ctl]
    r = lineage_sim.run_one(ctl, peers, a["reference"], a["contaminants"], a["region"],
                            lineage_sim.SCENARIOS[scen], seed, mode, depth=depth)
    r.update(amplicon=amp, scenario=scen)
    r.pop("roc", None)
    r.pop("lineages", None)
    return r


def ratio(est, true):
    if est is None or true is None or true == 0:
        return None
    return est / true


def report(runs):
    pairs = {"mutation_frequency": ("est_read_weighted_x1e3", "true_read_weighted_x1e3"),
             "unique_density": ("est_unique_density_x1e5", "true_unique_density_x1e5"),
             "combined_density": ("est_combined_density_x1e5", "true_unique_density_x1e5")}
    groups = defaultdict(list)
    for r in runs:
        groups[(r["amplicon"], r["scenario"], r["depth"], r["mode"])].append(r)
    out = ["# Quantification accuracy", "",
           "Estimate / truth, mean of the three held-out controls (1.00 = exact). For the "
           "null scenario, the estimates themselves (truth is 0).", "",
           "| Amplicon | Scenario | Depth | Mode | Mutation frequency | Unique density | "
           "Combined density | False calls | Recall, 2+ reads |",
           "|---|---|---:|---|---:|---:|---:|---:|---:|"]
    for key in sorted(groups, key=lambda k: (k[0], k[1], k[2], k[3])):
        rs = groups[key]
        cells = []
        for m in KEYS:
            e, t = pairs[m]
            if key[1] == "null":
                v = [r.get(e) for r in rs if r.get(e) is not None]
                cells.append(f"{statistics.mean(v):.2f} (est)" if v else "–")
            else:
                v = [ratio(r.get(e), r.get(t)) for r in rs]
                v = [x for x in v if x is not None]
                cells.append(f"{statistics.mean(v):.2f}" if v else "–")
        fc = statistics.mean(r["false_calls"] for r in rs)
        rc = [r["recall_2plus"] for r in rs if r.get("recall_2plus") is not None]
        out.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | " + " | ".join(cells)
                   + f" | {fc:.1f} | {statistics.mean(rc):.2f} |" if rc else
                   f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | " + " | ".join(cells)
                   + f" | {fc:.1f} | – |")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="quant_accuracy")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--quick", action="store_true", help="500 reads, one control, for a smoke test")
    ap.add_argument("--amplicons", default="JH4,N114")
    ap.add_argument("--modes", default=",".join(MODES))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tasks = []
    for amp in a.amplicons.split(","):
        for scen in AMPLICONS[amp]["scenarios"]:
            for depth in ((500,) if a.quick else DEPTHS):
                for mode in a.modes.split(","):
                    for k in ((0,) if a.quick else (0, 1, 2)):
                        tasks.append((amp, scen, depth, mode, k, 7 + 31 * k + depth))
    with ProcessPoolExecutor(a.jobs) as ex:
        runs = list(ex.map(task, tasks))
    json.dump(runs, open(os.path.join(a.out, "runs.json"), "w"), indent=1, default=str)
    md = report(runs)
    open(os.path.join(a.out, "report.md"), "w").write(md)
    print(md)


if __name__ == "__main__":
    main()
