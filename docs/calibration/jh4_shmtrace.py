"""JH4 ONT analysed with SHMtrace end to end, from the delivered reads.

    python docs/calibration/jh4_shmtrace.py [--full-length] [--input delivered|subsets]

Six libraries (Wt mice 2, 3, 4; Aicda-/- mice 1, 2, 3, the controls) in the recommended
(quality-calibrated) mode with the bundled JH4 intron reference and its paralogue reference.

  delivered (default) the reads as delivered by the sequencing provider (Original_FASTQ/).
            SHMtrace does all of the work: length and quality gates, competitive alignment
            against the paralogue, the paralogue-allele filter, a random draw of 500 usable
            reads, and calling against the controls. This is the manuscript's analysis.
  subsets   the deposit's 500-read subsets made by the earlier preprocessing script
            (Processed_FASTQ_Depth_Subsets/*_C57BL6_sub500.fastq); kept for comparison.

--full-length analyses positions 90-1,230 instead of the 529-bp window (Figure 2C).

Writes figures/jh4_shmtrace[_subsets][_full]/ (report figures, project.json, table_s1.json with
SHMtrace's read accounting per library) and, for the delivered 529-bp run, the calls in the
shape of a Geneious export in manuscript_corrected/jh4_calls/ for the figure scripts.
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
REF = ROOT / "src/shmtrace/references"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT")
LIBS = {"WT2": ("DHJGKB_2_Wildtype_2", "Wt"), "WT3": ("DHJGKB_3_Wildtype_3", "Wt"),
        "WT4": ("DHJGKB_4_Wildtype_4", "Wt"), "KO1": ("DHJGKB_5_AID_KO_1", "Aicda-/-"),
        "KO2": ("DHJGKB_6_AID_KO_2", "Aicda-/-"), "KO3": ("DHJGKB_7_AID_KO_3", "Aicda-/-")}


def fastq(f, source):
    if source == "delivered":
        return DEP / "Original_FASTQ" / f"{f}.fastq"
    return DEP / "Processed_FASTQ_Depth_Subsets" / f"{f}_C57BL6_sub500.fastq"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", choices=("delivered", "subsets"), default="delivered")
    ap.add_argument("--full-length", action="store_true")
    a = ap.parse_args()
    name = "jh4_shmtrace" + ("_subsets" if a.input == "subsets" else "") + ("_full" if a.full_length else "")
    out = Path(__file__).parent / "figures" / name
    overrides = {"analysis_window": [90, 1230]} if a.full_length else {}
    contaminants = [str(REF / "mouse_jh4_contaminant.fasta")]
    proj = {"name": f"JH4 intron, SHMtrace ({a.input} reads)", "region": "jh4_intron",
            "reference": str(REF / "mouse_jh4_intron.fasta"), "contaminants": contaminants,
            "mode": "recommended", "depth": 500, "overrides": overrides,
            "samples": [{"name": n, "path": str(fastq(f, a.input)), "group": g, "control": g == "Aicda-/-"}
                        for n, (f, g) in LIBS.items()]}
    res = project.run_project(proj)
    out.mkdir(parents=True, exist_ok=True)
    report.write_all(res, str(out))
    json.dump(project.public(res), open(out / "project.json", "w"), default=str)
    # SHMtrace's read accounting per library, and candidates tested
    acct = {}
    for n, (f, _g) in LIBS.items():
        acct[n] = account(fastq(f, a.input), REF / "mouse_jh4_intron.fasta", contaminants, "jh4_intron",
                          overrides=overrides)
        sm = res["_details"][n]["summary"]
        acct[n].update(candidates_tested=sm["calling"].get("candidates_tested"),
                       calls=sm["calling"].get("unique_variants"))
    json.dump(acct, open(out / "table_s1.json", "w"), indent=1)
    if a.input == "delivered" and not a.full_length:
        calls = Path(__file__).parent / "manuscript_corrected/jh4_calls"
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
        print(f"{s['sample']:3} {s['group']:9} reads {s['reads']:4} calls {s['unique_mutations']:4} "
              f"density {s['unique_density_x1e5']:7.2f}  freq/1e3 {s['mutation_frequency_x1e3']:.3f}")
    print("accounting", json.dumps(acct))
    print("calibration", res.get("calibration"))


if __name__ == "__main__":
    main()
