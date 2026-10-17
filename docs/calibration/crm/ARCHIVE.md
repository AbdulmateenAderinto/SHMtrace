# Archive: what produced each result

Commands are run from `docs/calibration`. Paths are relative to it. `out/` holds intermediates that
the figure and manuscript builders read; the figures and the manuscript are rebuilt from them, so a
number in the text and the same number in a figure always come from one file.

The environment this was last recorded in is in `ARCHIVE_ENVIRONMENT.txt`:

```
python 3.13.12
numpy 2.5.4
scipy 1.18.1
pysam 0.24.1
mappy 2.31
biopython 1.88
matplotlib 3.11.2
openpyxl 3.1.5
python-docx 1.2.0
edlib (not installed)
minimap2 2.30-r1287
samtools samtools 1.21
lofreq (not on PATH)
git 8d30d61416260f7f899fdcae7c473aa91f4e2ed6 (uncommitted changes present)
```

| Produces | Command | Inputs | Outputs | Output present |
|---|---|---|---|---|
| J_H4 and N114 nanopore projects (Figures 2, 4-6; Tables S1, S3) | `python crm/jh4_extras.py` | the deposit's FASTQ files; src/shmtrace/references/ | `figures/jh4_shmtrace/, figures/n114_shmtrace/, figures/jh4_shmtrace_full/` | present |
| Caller comparison on real libraries (Figure 4; Table S2) | `python crm/benchmark_real.py jh4 n114 su_lps su_ch12 pim1 ot1 ot2` | the same FASTQ files; LoFreq, ClairS, Clair3, DeepVariant, DeepSomatic | `crm/out/benchmark_real_*.json` | present |
| Planted mutations (Figure 3A-C, 3H; Table S2) | `python crm/spikein.py jh4 n114 su_lps pim1` | the AID-null libraries of each locus | `crm/out/spikein_*.json` | present |
| Simulated lineages (Figure 3D-F) | `shmtrace-simulate-lineages --region jh4_intron --control KO1.fastq --control KO2.fastq --control KO3.fastq --reference mouse_jh4_intron.fasta --contaminant mouse_jh4_contaminant.fasta --targeting-model HH_S5F.csv --out lineage_sim` | the three J_H4 AID-null libraries; HH_S5F.csv (SHazaM, not redistributed) | `crm/out/lineage_sim_jh4/` | present |
| The high-allele-fraction rule on a mutated dominant clone (Results, STAR Methods) | `shmtrace-simulate-lineages ... --scenarios expanded,expanded_mutated_40,expanded_mutated_60,expanded_mutated_80 --modes recommended --out mask_sim` | the same three AID-null libraries | `crm/out/mask_sim_jh4/` | present |
| Depth series with repeated draws (Figure S5) | `python crm/depth_draws.py jh4 5; python crm/depth_draws.py n114 5` | the J_H4 and N114 FASTQ files | `crm/out/depth_draws_jh4.json, crm/out/depth_draws_n114.json` | present |
| Power with the pilot CV's uncertainty (Figure S3B) | `python crm/power_s3.py` | figures/jh4_shmtrace/project.json (the wild-type pilot) | `crm/out/power_s3.json` | present |
| Group-test calibration (Figure 3G) | `python replicates_calibration.py` | simulated counts | `replicates_calibration/runs.json` | present |
| Sanger of the same mice (Figure 5) | `python crm/sanger_cohort2.py` | the deposit's chromatograms; samplesheets/ | `figures/sanger_vs_nanopore/` | present |
| Published Sanger study, Maea (Figure S8) | `python crm/sanger_cohort2.py cohort2` | samplesheets/cohort2_barbulescu2024.csv and its chromatograms | `crm/out/sanger_cohort2.json` | present |
| Published Sanger study, Mkln1 (Figure S9) | `python crm/sanger_cohort2.py mkln` | samplesheets/jh4_mkln_barbulescu2025.csv and its chromatograms | `crm/out/sanger_mkln.json` | present |
| External Sanger study, SW_HEL (Figure S10) | `python external/swhel_samhd1.py <genbank-clone-dir>` | GenBank MH198782-MH200576 | `external/out/swhel_samhd1.json` | present |
| The two published studies metric by metric (Table S4) | `python crm/sanger_compare.py` | the two runs above; crm/published/*.json | `(printed rows, read by the manuscript builder)` | printed |
| Public Illumina and nanopore datasets (Figure 7) | `python crm/public_illumina.py` | PRJNA973112, PRJNA918596, PRJNA820976, PRJNA170048, Karst et al. | `crm/out/public/` | present |
| Figure S12's nanopore reads | `python crm/scicanvas/s12_reads.py` | the J_H4 libraries | `(printed READS array, pasted into the builder)` | printed |
| Figures | `python crm/figs/fig2_call.py (and fig3..fig7, figS2, figS3, figS6, figS8, figS9)` | the outputs above | `crm/figures/Figure*.png, .pdf, .svg, .tif` | present |
| SciCanvas figures | `node crm/scicanvas/build.js crm/scicanvas/<builder>.js crm/figures/<name>` | the outputs above | `crm/figures/Figure1.pdf, and FigureS1_workflow, S10, S12, S13, S14, GraphicalAbstract` | present |
| Mutation trees (Figures 6D and S11) | `python mutation_tree_figure.py` | the J_H4 FASTQ files | `figures/mutation_trees/` | present |
| Table S2 workbook | `python crm/make_tables.py` | crm/out/benchmark_real_*.json, spikein_*.json | `crm/out/Table_S2_benchmark.xlsx` | present |
| The manuscript and supplement | `python crm/manuscript/build_crm.py crm/build` | everything above; crm/manuscript/*.py | `crm/build/*.docx` | present |
| The reviewer package | `python crm/make_review_package.py` | crm/build/*.docx, crm/figures/ | `crm/build/SHMtrace_review_package.zip` | present |

Not redistributed: the SHazaM HH_S5F model (AGPL) and the Karst et al. rRNA reads (GPL-3) are
fetched from their own sources; every other input is in the data deposit or a public BioProject.
