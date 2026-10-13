"""N114 V region analysed with SHMtrace end to end, from the delivered reads.

    python docs/calibration/n114_shmtrace.py [--input delivered|subsets]

Six libraries (pMX-Aicda subclones C1, D3, D5; empty-vector C3, C4, D1, the controls) in the
recommended (quality-calibrated) mode with the bundled N114 reference, 500 reads per library.

  delivered (default) the reads as delivered (Original_FASTQ/): SHMtrace applies its own length
            and quality gates, alignment, a random draw of 500 usable reads and the calling.
  subsets   the deposit's 500-read subsets from the earlier preprocessing script; for comparison.

Writes figures/n114_shmtrace[_subsets]/ (report figures, project.json, table_s1.json) and, for
the delivered run, Geneious-shaped call CSVs in manuscript_corrected/n114_calls/.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

from shmtrace import project, report

sys.path.insert(0, str(Path(__file__).parent))
from read_accounting import account  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "src/shmtrace/references/n114_vregion.fasta"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/03_N114_ONT")
LIBS = {"C1": ("75K3B3_6_C1_AID_10.06.24", "pMX-Aicda"),
        "D3": ("75K3B3_7_D3_AID_10.01.24", "pMX-Aicda"),
        "D5": ("75K3B3_10_D5_AID_6.28.24_", "pMX-Aicda"),
        "C3": ("75K3B3_1_C3_Empty_March_2024_1", "pMX-Empty"),
        "C4": ("75K3B3_3_C4_Empty_March_2024", "pMX-Empty"),
        "D1": ("75K3B3_4_D1_Empty_March_2024", "pMX-Empty")}


def fastq(f, source):
    if source == "delivered":
        return DEP / "Original_FASTQ" / f"{f}.fastq"
    return DEP / "Processed_FASTQ_Depth_Subsets" / f"{f}_N114_sub500.fastq"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", choices=("delivered", "subsets"), default="delivered")
    a = ap.parse_args()
    out = Path(__file__).parent / "figures" / ("n114_shmtrace" + ("_subsets" if a.input == "subsets" else ""))
    proj = {"name": f"N114 V region, SHMtrace ({a.input} reads)", "region": "n114_vregion",
            "reference": str(REF), "contaminants": [], "mode": "recommended", "depth": 500,
            "samples": [{"name": n, "path": str(fastq(f, a.input)), "group": g, "control": g == "pMX-Empty"}
                        for n, (f, g) in LIBS.items()]}
    res = project.run_project(proj)
    out.mkdir(parents=True, exist_ok=True)
    report.write_all(res, str(out))
    json.dump(project.public(res), open(out / "project.json", "w"), default=str)
    acct = {}
    for n, (f, _g) in LIBS.items():
        acct[n] = account(fastq(f, a.input), REF, [], "n114_vregion")
        sm = res["_details"][n]["summary"]
        acct[n].update(candidates_tested=sm["calling"].get("candidates_tested"),
                       calls=sm["calling"].get("unique_variants"))
    json.dump(acct, open(out / "table_s1.json", "w"), indent=1)
    if a.input == "delivered":
        calls = Path(__file__).parent / "manuscript_corrected/n114_calls"
        calls.mkdir(exist_ok=True)
        for n in LIBS:
            with open(calls / f"{n} SHMtrace Annotations.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Name", "Minimum", "Polymorphism Type", "Reference Nucleotide(s)",
                            "Variant Nucleotide(s)", "Variant Frequency", "Variant P-Value (approximate)",
                            "Coverage", "Supporting reads"])
                for v in res["variants"]:
                    if v["sample"] != n:
                        continue
                    kind = "transition" if v["ref"] + v["alt"] in ("AG", "GA", "CT", "TC") else "transversion"
                    w.writerow([f"{v['ref']}->{v['alt']}", v["pos"], f"SNP ({kind})", v["ref"], v["alt"],
                                f"{100 * v['frequency']:.2f}%", "0", v["depth"], v["reads"]])
    for s in res["samples"]:
        print(f"{s['sample']:3} {s['group']:10} reads {s['reads']:4} calls {s['unique_mutations']:4} "
              f"density {s['unique_density_x1e5']:7.2f}  freq/1e3 {s['mutation_frequency_x1e3']:.3f}")
    print("accounting", json.dumps(acct))
    print("calibration", res.get("calibration"))


if __name__ == "__main__":
    main()
