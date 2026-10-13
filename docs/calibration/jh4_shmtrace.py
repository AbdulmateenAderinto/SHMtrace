"""JH4 ONT re-analysed with the current SHMtrace pipeline, from the deposited reads.

    python docs/calibration/jh4_shmtrace.py [--input selected|raw]

Six libraries (Wt mice 2, 3, 4; Aicda-/- mice 1, 2, 3, the controls), analysed at 500 reads
per library (the manuscript's depth) in the recommended (quality-calibrated) mode with the
bundled JH4 intron reference and its paralogue (contaminant) reference.

  selected  (default) the deposit's realised 500-read subsets of the target-selected reads,
            Processed_FASTQ_Depth_Subsets/*_C57BL6_sub500.fastq: the inputs Methods 2.7
            names, as for N114, so SHMtrace replaces only the Geneious calling.
  raw       the delivered reads, Original_FASTQ/: SHMtrace end to end, its own paralogue
            filter doing the target selection. A sensitivity check.

--full-length analyses positions 90-1230 instead of the 529-bp window (Figure 2C).

Writes figures/jh4_shmtrace[_raw][_full]/ (report figures, project.json) and, for `selected`,
Geneious-shaped call CSVs in manuscript_corrected/jh4_calls/ for the figure scripts.
"""
import argparse
import csv
import json
from pathlib import Path

from shmtrace import project, report

ROOT = Path(__file__).resolve().parents[2]
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT")
LIBS = {"WT2": ("DHJGKB_2_Wildtype_2", "Wt"), "WT3": ("DHJGKB_3_Wildtype_3", "Wt"),
        "WT4": ("DHJGKB_4_Wildtype_4", "Wt"), "KO1": ("DHJGKB_5_AID_KO_1", "Aicda-/-"),
        "KO2": ("DHJGKB_6_AID_KO_2", "Aicda-/-"), "KO3": ("DHJGKB_7_AID_KO_3", "Aicda-/-")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", choices=("selected", "raw"), default="selected")
    ap.add_argument("--full-length", action="store_true",
                    help="analyse positions 90-1230 (the full-length map, Figure 2C) "
                         "instead of the 529-bp Sanger-comparison window")
    a = ap.parse_args()
    if a.input == "selected":
        path = lambda f: DEP / "Processed_FASTQ_Depth_Subsets" / f"{f}_C57BL6_sub500.fastq"  # noqa: E731
        out = Path(__file__).parent / "figures/jh4_shmtrace"
    else:
        path = lambda f: DEP / "Original_FASTQ" / f"{f}.fastq"  # noqa: E731
        out = Path(__file__).parent / "figures/jh4_shmtrace_raw"
    if a.full_length:
        out = out.with_name(out.name + "_full")
    ref = ROOT / "src/shmtrace/references"
    proj = {"name": f"JH4 intron, SHMtrace ({a.input} reads)", "region": "jh4_intron",
            "reference": str(ref / "mouse_jh4_intron.fasta"),
            "contaminants": [str(ref / "mouse_jh4_contaminant.fasta")],
            "mode": "recommended", "depth": 500,
            "overrides": {"analysis_window": [90, 1230]} if a.full_length else {},
            "samples": [{"name": n, "path": str(path(f)), "group": g, "control": g == "Aicda-/-"}
                        for n, (f, g) in LIBS.items()]}
    res = project.run_project(proj)
    out.mkdir(parents=True, exist_ok=True)
    report.write_all(res, str(out))
    json.dump(project.public(res), open(out / "project.json", "w"), default=str)
    # read accounting and candidates tested per library (Supplementary Table S1)
    acct = {}
    for name, d in res["_details"].items():
        sm = d["summary"]
        acct[name] = {"reads_in": sm["reads"].get("reads_in"), "passed": sm["reads"].get("passed_prefilter"),
                      "aligned": sm["reads"].get("aligned_to_target"),
                      "candidates_tested": sm["calling"].get("candidates_tested"),
                      "calls": sm["calling"].get("unique_variants")}
    json.dump(acct, open(out / "table_s1.json", "w"), indent=1)
    print(acct)
    if a.input == "selected" and not a.full_length:
        # the calls in the shape of a Geneious annotation export, for the figure scripts
        calls = Path(__file__).parent / "manuscript_corrected/jh4_calls"
        calls.mkdir(exist_ok=True)
        for name in LIBS:
            with open(calls / f"{name} SHMtrace Annotations.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Name", "Minimum", "Polymorphism Type", "Reference Nucleotide(s)",
                            "Variant Nucleotide(s)", "Variant Frequency",
                            "Variant P-Value (approximate)", "Coverage", "Supporting reads"])
                for v in res["variants"]:
                    if v["sample"] != name:
                        continue
                    kind = "transition" if v["ref"] + v["alt"] in ("AG", "GA", "CT", "TC") else "transversion"
                    w.writerow([f"{v['ref']}->{v['alt']}", v["pos"], f"SNP ({kind})", v["ref"], v["alt"],
                                f"{100 * v['frequency']:.2f}%", "0", v["depth"], v["reads"]])
    for s in res["samples"]:
        print(f"{s['sample']:3} {s['group']:9} reads {s['reads']:4} calls {s['unique_mutations']:4} "
              f"unique density {s['unique_density_x1e5']:7.2f}  combined {s.get('combined_density_x1e5') or 0:7.2f}"
              f"  freq/1e3 {s['mutation_frequency_x1e3']:.3f}")
    for t in res["animal_tests"]:
        print(t["metric"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in t.items()
                            if k in ("fold", "fold_ci", "p_negbinomial", "p_welch_log", "p_permutation",
                                     "fold_is_lower_bound")})


if __name__ == "__main__":
    main()
