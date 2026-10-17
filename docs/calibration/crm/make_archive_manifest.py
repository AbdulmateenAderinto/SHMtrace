"""Record the environment and the command behind every displayed result.

    python docs/calibration/crm/make_archive_manifest.py   -> ARCHIVE.md, ARCHIVE_ENVIRONMENT.txt

A reader should be able to see, without reading the prose, which command produced each figure, table and
quoted number, what it read and what it wrote. This writes two files beside the manuscript sources:

  ARCHIVE_ENVIRONMENT.txt  the interpreter, the installed packages, the external tools with their versions,
                           and the git commit the manuscript was built from
  ARCHIVE.md               one row per analysis: command, inputs, outputs, and what it feeds

Paths are given relative to docs/calibration, and outputs under out/ are intermediates the figure and
manuscript builders read. Nothing here is generated from prose: the rows are declared below and checked
against the filesystem, so a stale row shows up as "missing" rather than passing silently.
"""
import importlib.metadata as md
import subprocess
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
CAL = CRM.parent
PKGS = ["numpy", "scipy", "pysam", "mappy", "biopython", "matplotlib", "openpyxl", "python-docx", "edlib"]
TOOLS = [("minimap2", ["minimap2", "--version"]), ("samtools", ["samtools", "--version"]),
         ("lofreq", ["lofreq", "version"])]

# (what it produces, command, inputs, outputs). Commands are run from docs/calibration.
ANALYSES = [
    ("J_H4 and N114 nanopore projects (Figures 2, 4-6; Tables S1, S3)",
     "python crm/jh4_extras.py", "the deposit's FASTQ files; src/shmtrace/references/",
     "figures/jh4_shmtrace/, figures/n114_shmtrace/, figures/jh4_shmtrace_full/"),
    ("Caller comparison on real libraries (Figure 4; Table S2)",
     "python crm/benchmark_real.py jh4 n114 su_lps su_ch12 pim1 ot1 ot2",
     "the same FASTQ files; LoFreq, ClairS, Clair3, DeepVariant, DeepSomatic",
     "crm/out/benchmark_real_*.json"),
    ("Planted mutations (Figure 3A-C, 3H; Table S2)",
     "python crm/spikein.py jh4 n114 su_lps pim1", "the AID-null libraries of each locus",
     "crm/out/spikein_*.json"),
    ("Simulated lineages (Figure 3D-F)",
     "shmtrace-simulate-lineages --region jh4_intron --control KO1.fastq --control KO2.fastq "
     "--control KO3.fastq --reference mouse_jh4_intron.fasta --contaminant mouse_jh4_contaminant.fasta "
     "--targeting-model HH_S5F.csv --out lineage_sim",
     "the three J_H4 AID-null libraries; HH_S5F.csv (SHazaM, not redistributed)",
     "crm/out/lineage_sim_jh4/"),
    ("The high-allele-fraction rule on a mutated dominant clone (Results, STAR Methods)",
     "shmtrace-simulate-lineages ... --scenarios expanded,expanded_mutated_40,expanded_mutated_60,"
     "expanded_mutated_80 --modes recommended --out mask_sim",
     "the same three AID-null libraries", "crm/out/mask_sim_jh4/"),
    ("Depth series with repeated draws (Figure S5)", "python crm/depth_draws.py jh4 5; python crm/depth_draws.py n114 5",
     "the J_H4 and N114 FASTQ files", "crm/out/depth_draws_jh4.json, crm/out/depth_draws_n114.json"),
    ("Power with the pilot CV's uncertainty (Figure S3B)", "python crm/power_s3.py",
     "figures/jh4_shmtrace/project.json (the wild-type pilot)", "crm/out/power_s3.json"),
    ("Group-test calibration (Figure 3G)", "python replicates_calibration.py",
     "simulated counts", "replicates_calibration/runs.json"),
    ("Sanger of the same mice (Figure 5)", "python crm/sanger_cohort2.py",
     "the deposit's chromatograms; samplesheets/", "figures/sanger_vs_nanopore/"),
    ("Published Sanger study, Maea (Figure S8)", "python crm/sanger_cohort2.py cohort2",
     "samplesheets/cohort2_barbulescu2024.csv and its chromatograms", "crm/out/sanger_cohort2.json"),
    ("Published Sanger study, Mkln1 (Figure S9)", "python crm/sanger_cohort2.py mkln",
     "samplesheets/jh4_mkln_barbulescu2025.csv and its chromatograms", "crm/out/sanger_mkln.json"),
    ("External Sanger study, SW_HEL (Figure S10)", "python external/swhel_samhd1.py <genbank-clone-dir>",
     "GenBank MH198782-MH200576", "external/out/swhel_samhd1.json"),
    ("The two published studies metric by metric (Table S4)", "python crm/sanger_compare.py",
     "the two runs above; crm/published/*.json", "(printed rows, read by the manuscript builder)"),
    ("Public Illumina and nanopore datasets (Figure 7)", "python crm/public_illumina.py",
     "PRJNA973112, PRJNA918596, PRJNA820976, PRJNA170048, Karst et al.", "crm/out/public/"),
    ("Figure S12's nanopore reads", "python crm/scicanvas/s12_reads.py",
     "the J_H4 libraries", "(printed READS array, pasted into the builder)"),
    ("Figures", "python crm/figs/fig2_call.py (and fig3..fig7, figS2, figS3, figS6, figS8, figS9)",
     "the outputs above", "crm/figures/Figure*.png, .pdf, .svg, .tif"),
    ("SciCanvas figures", "node crm/scicanvas/build.js crm/scicanvas/<builder>.js crm/figures/<name>",
     "the outputs above", "crm/figures/Figure1.pdf, and FigureS1_workflow, S10, S12, S13, S14, GraphicalAbstract"),
    ("Mutation trees (Figures 6D and S11)", "python mutation_tree_figure.py",
     "the J_H4 FASTQ files", "figures/mutation_trees/"),
    ("Table S2 workbook", "python crm/make_tables.py", "crm/out/benchmark_real_*.json, spikein_*.json",
     "crm/out/Table_S2_benchmark.xlsx"),
    ("The manuscript and supplement", "python crm/manuscript/build_crm.py crm/build",
     "everything above; crm/manuscript/*.py", "crm/build/*.docx"),
    ("The reviewer package", "python crm/make_review_package.py", "crm/build/*.docx, crm/figures/",
     "crm/build/SHMtrace_review_package.zip"),
]


def environment() -> str:
    out = [f"python {sys.version.split()[0]}"]
    for p in PKGS:
        try:
            out.append(f"{p} {md.version(p)}")
        except md.PackageNotFoundError:
            out.append(f"{p} (not installed)")
    for name, cmd in TOOLS:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
            out.append(f"{name} {(r.stdout or r.stderr).splitlines()[0].strip()}")
        except Exception:
            out.append(f"{name} (not on PATH)")
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=CRM, capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=CRM, capture_output=True, text=True).stdout
        out.append(f"git {sha}{' (uncommitted changes present)' if dirty.strip() else ''}")
    except Exception:
        pass
    return "\n".join(out)


def main():
    env = environment()
    (CRM / "ARCHIVE_ENVIRONMENT.txt").write_text(env + "\n")
    rows = []
    for what, cmd, ins, outs in ANALYSES:
        # an output is "present" when its first path exists; a stale row is visible rather than silent
        first = outs.split(",")[0].strip()
        if first.startswith("("):                      # the command prints its result rather than writing a file
            state = "printed"
        elif "*" in first:
            state = "present" if list(CAL.glob(first)) else "missing"
        else:
            state = "present" if (CAL / first.rstrip("/")).exists() else "missing"
        rows.append((what, cmd, ins, outs, state))
    lines = ["# Archive: what produced each result", "",
             "Commands are run from `docs/calibration`. Paths are relative to it. `out/` holds intermediates that",
             "the figure and manuscript builders read; the figures and the manuscript are rebuilt from them, so a",
             "number in the text and the same number in a figure always come from one file.", "",
             "The environment this was last recorded in is in `ARCHIVE_ENVIRONMENT.txt`:", "", "```", env, "```", "",
             "| Produces | Command | Inputs | Outputs | Output present |", "|---|---|---|---|---|"]
    lines += [f"| {w} | `{c}` | {i} | `{o}` | {s} |" for w, c, i, o, s in rows]
    lines += ["", "Not redistributed: the SHazaM HH_S5F model (AGPL) and the Karst et al. rRNA reads (GPL-3) are",
              "fetched from their own sources; every other input is in the data deposit or a public BioProject."]
    (CRM / "ARCHIVE.md").write_text("\n".join(lines) + "\n")
    missing = [w for w, _c, _i, _o, s in rows if s == "missing"]
    print(f"ARCHIVE.md: {len(rows)} analyses, {len(missing)} with outputs not on this machine")
    for m in missing:
        print("  missing:", m)


if __name__ == "__main__":
    main()
