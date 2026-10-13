"""N114 re-analysed with the current SHMtrace pipeline, from the deposited reads.

    python docs/calibration/n114_shmtrace.py

Six libraries (pMX-Aicda subclones C1, D3, D5; pMX-Empty C3, C4, D1, the controls), all
reads of each deposited FASTQ, analysed at 500 reads per library (the manuscript's depth)
in the recommended (quality-calibrated) mode with the bundled N114 reference. Writes the
standard project figures and tables to figures/n114_shmtrace/ and prints the comparison
with the manuscript's Geneious-based numbers (Figure 4: density 30.6 vs 1.69, 18-fold).
"""
import json
from pathlib import Path

from shmtrace import project, report

ROOT = Path(__file__).resolve().parents[2]
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/03_N114_ONT/"
                     "Processed_FASTQ_Depth_Subsets")
LIBS = {"C1": ("75K3B3_6_C1_AID_10.06.24_N114", "pMX-Aicda"),
        "D3": ("75K3B3_7_D3_AID_10.01.24_N114", "pMX-Aicda"),
        "D5": ("75K3B3_10_D5_AID_6.28.24__N114", "pMX-Aicda"),
        "C3": ("75K3B3_1_C3_Empty_March_2024_1_N114", "pMX-Empty"),
        "C4": ("75K3B3_3_C4_Empty_March_2024_N114", "pMX-Empty"),
        "D1": ("75K3B3_4_D1_Empty_March_2024_N114", "pMX-Empty")}
OUT = Path(__file__).parent / "figures/n114_shmtrace"


def main():
    proj = {"name": "N114 V region, SHMtrace", "region": "n114_vregion",
            "reference": str(ROOT / "src/shmtrace/references/n114_vregion.fasta"),
            "contaminants": [], "mode": "recommended", "depth": 500,
            "samples": [{"name": n, "path": str(DEP / f"{f}.fastq"), "group": g,
                         "control": g == "pMX-Empty"} for n, (f, g) in LIBS.items()]}
    res = project.run_project(proj)
    OUT.mkdir(parents=True, exist_ok=True)
    report.write_all(res, str(OUT))
    json.dump(project.public(res), open(OUT / "project.json", "w"), default=str)
    # the calls in the shape of a Geneious annotation export, for the manuscript's figure
    # scripts (Supplementary S4 reads these; positions are on the 655-nt N114 reference)
    import csv
    calls = Path(__file__).parent / "manuscript_corrected/n114_calls"
    calls.mkdir(exist_ok=True)
    for name in LIBS:
        with open(calls / f"{name} SHMtrace Annotations.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Name", "Minimum", "Polymorphism Type", "Reference Nucleotide(s)",
                        "Variant Nucleotide(s)", "Variant Frequency", "Variant P-Value (approximate)",
                        "Coverage", "Supporting reads"])
            for v in res["variants"]:
                if v["sample"] != name:
                    continue
                kind = "transition" if v["ref"] + v["alt"] in ("AG", "GA", "CT", "TC") else "transversion"
                w.writerow([f"{v['ref']}->{v['alt']}", v["pos"], f"SNP ({kind})", v["ref"], v["alt"],
                            f"{100 * v['frequency']:.2f}%", "0", v["depth"], v["reads"]])
    for s in res["samples"]:
        print(f"{s['sample']:3} {s['group']:10} reads {s['reads']:4} calls {s['unique_mutations']:4} "
              f"unique density {s['unique_density_x1e5']:7.2f}  combined {s.get('combined_density_x1e5') or 0:7.2f}"
              f"  freq/1e3 {s['mutation_frequency_x1e3']:.3f}")
    for t in res["animal_tests"]:
        print(t["metric"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in t.items()
                            if k in ("fold", "fold_ci", "p_negbinomial", "p_welch_log", "p_permutation",
                                     "fold_is_lower_bound")})


if __name__ == "__main__":
    main()
