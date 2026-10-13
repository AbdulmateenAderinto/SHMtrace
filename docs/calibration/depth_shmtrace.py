"""Depth series for Supplementary S6/S7 with SHMtrace, from the deposited depth subsets.

    python docs/calibration/depth_shmtrace.py {jh4|n114}

The manuscript's down-sampling design is kept: one independently drawn, non-nested subsample
per library at 30, 60, 125, 250, 500 and 1,000 reads (the deposit's *_subN.fastq files). Each
depth is one SHMtrace project (recommended mode), so the controls are calibrated at the same
depth as the cases. Writes figures/depth_shmtrace/<assay>.json:
{depth: {sample: {"reads", "calls", "analysed_bp", "density"}}}.
"""
import json
import sys
from pathlib import Path

from shmtrace import project

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "src/shmtrace/references"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data")
DEPTHS = [30, 60, 125, 250, 500, 1000]
ASSAYS = {
    "jh4": dict(region="jh4_intron", reference=REF / "mouse_jh4_intron.fasta",
                contaminants=[str(REF / "mouse_jh4_contaminant.fasta")],
                dir=DEP / "02_JH4_ONT/Processed_FASTQ_Depth_Subsets", suffix="_C57BL6",
                libs={"WT2": ("DHJGKB_2_Wildtype_2", "Wt"), "WT3": ("DHJGKB_3_Wildtype_3", "Wt"),
                      "WT4": ("DHJGKB_4_Wildtype_4", "Wt"), "KO1": ("DHJGKB_5_AID_KO_1", "Aicda-/-"),
                      "KO2": ("DHJGKB_6_AID_KO_2", "Aicda-/-"), "KO3": ("DHJGKB_7_AID_KO_3", "Aicda-/-")},
                control="Aicda-/-"),
    "n114": dict(region="n114_vregion", reference=REF / "n114_vregion.fasta", contaminants=[],
                 dir=DEP / "03_N114_ONT/Processed_FASTQ_Depth_Subsets", suffix="",
                 libs={"C1": ("75K3B3_6_C1_AID_10.06.24_N114", "pMX-Aicda"),
                       "D3": ("75K3B3_7_D3_AID_10.01.24_N114", "pMX-Aicda"),
                       "D5": ("75K3B3_10_D5_AID_6.28.24__N114", "pMX-Aicda"),
                       "C3": ("75K3B3_1_C3_Empty_March_2024_1_N114", "pMX-Empty"),
                       "C4": ("75K3B3_3_C4_Empty_March_2024_N114", "pMX-Empty"),
                       "D1": ("75K3B3_4_D1_Empty_March_2024_N114", "pMX-Empty")},
                 control="pMX-Empty"),
}


def main():
    name = sys.argv[1]
    a = ASSAYS[name]
    out = {}
    for d in DEPTHS:
        proj = {"name": f"{name} depth {d}", "region": a["region"], "reference": str(a["reference"]),
                "contaminants": a["contaminants"], "mode": "recommended", "depth": d,
                "samples": [{"name": n, "path": str(a["dir"] / f"{f}{a['suffix']}_sub{d}.fastq"),
                             "group": g, "control": g == a["control"]} for n, (f, g) in a["libs"].items()]}
        res = project.run_project(proj)
        out[d] = {s["sample"]: {"group": s["group"], "reads": s["reads"], "calls": s["unique_mutations"],
                                "analysed_bp": s["analysed_bp"], "density": s["unique_density_x1e5"]}
                  for s in res["samples"]}
        print(d, {k: v["calls"] for k, v in out[d].items()}, flush=True)
    od = Path(__file__).parent / "figures/depth_shmtrace"
    od.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(od / f"{name}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
