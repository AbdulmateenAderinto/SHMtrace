"""Assemble the Ver21 submission package and the Ver11 data deposit beside the earlier ones.

    python build_submission_v21.py          (needs python-docx and Pillow; LibreOffice for PDFs)

Nothing earlier is modified: the Ver13 package (cover letter, highlights) and the Ver10 deposit are
read and copied.
  ONT_SHM_JIM_Submission_Ver21_<date>/   manuscript, cover letter, highlights, figures, graphical
                                         abstract, code, data, supporting notes
  ONT_SHM_DataDeposit_Ver11_<date>/      Ver10 with the SHMtrace analyses replaced by the runs on the
                                         delivered reads (the manuscript's), the earlier subset runs
                                         kept for comparison, SHMtrace's references and read accounting
Run revise_manuscript_v21.py first, and commit, so the source archive matches the analysis.
"""
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import date
from pathlib import Path

from docx import Document
from PIL import Image

HERE = Path(__file__).resolve().parent
CAL = HERE.parent
REPO = CAL.parents[1]
FIG = HERE / "figures"
RUNS_DIR = CAL / "figures"
SV = RUNS_DIR / "sanger_vs_nanopore"
SC = HERE / "figure1_scicanvas/out"
SUB = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4")
V13 = SUB / "ONT_SHM_JIM_Submission_Ver13_2026-08-01"
DEP10 = SUB / "ONT_SHM_DataDeposit_Ver10_2026-10-12"
MS = SUB / "ONT_SHM_Manuscript_JIM_Ver21_shmtrace.docx"
SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
STAMP = date.today().isoformat()
PKG = SUB / f"ONT_SHM_JIM_Submission_Ver21_{STAMP}"
DEP = SUB / f"ONT_SHM_DataDeposit_Ver11_{STAMP}"
TITLE = ("Cloning-Free Profiling of Somatic Hypermutation with Oxford Nanopore R10.4.1/Kit 14 Sequencing: "
         "Proof of Concept in the Mouse JH4 Intron and an Immunoglobulin V Region")

# manuscript figure -> source stem; SciCanvas figures (1, S11) also have a 1,000-dpi TIFF
FIGURES = {
    "Figure1": SC / "Figure1_overview",
    "Figure2": FIG / "Figure2_JH4_combined",
    "Figure3": FIG / "Figure3_concordance",
    "Figure4": FIG / "Figure4_hybridoma_combined",
    "FigureS1": FIG / "FigureS_spectrum_529_vs_fulllength",
    "FigureS2": FIG / "FigureS5_individual_Sanger_matrices",
    "FigureS3": FIG / "FigureS6_individual_JH4_ONT_matrices",
    "FigureS4": FIG / "FigureS7_individual_N114_matrices",
    "FigureS5": FIG / "FigureS7_recurrence_weighted",
    "FigureS6": FIG / "FigureS3_depth_discovery_density",
    "FigureS7": FIG / "FigureS4_saturation",
    "FigureS8": SV / "figS_sanger_window_original",
    "FigureS9": SV / "recommended_manuscript/fig3_agreement_position_hotspot",
    "FigureS10": RUNS_DIR / "mutation_trees/mutation_trees_wt3",
    "FigureS11": SC / "Figure_workflow_detailed",
}
# SHMtrace analyses: run folder -> deposit name. The manuscript's runs, on the delivered reads ...
RUNS = {"jh4_shmtrace": "JH4_529bp_window", "jh4_shmtrace_full": "JH4_full_length_90-1230",
        "n114_shmtrace": "N114_V_region", "depth_shmtrace": "Depth_series",
        "sanger_vs_nanopore/recommended_manuscript": "Sanger_vs_ONT", "mutation_trees": "Mutation_trees"}
# ... and the same analyses on the earlier preprocessing script's 500-read subsets, for comparison
SUBSET_RUNS = {"jh4_shmtrace_subsets": "JH4_529bp_window", "jh4_shmtrace_full_subsets": "JH4_full_length_90-1230",
               "n114_shmtrace_subsets": "N114_V_region", "depth_shmtrace_subsets": "Depth_series",
               "sanger_vs_nanopore/recommended_subsets": "Sanger_vs_ONT"}
ONT_ROWS = {"Wt 2": ("jh4", "WT2"), "Wt 3": ("jh4", "WT3"), "Wt 4": ("jh4", "WT4"),
            "Aicda-null 1": ("jh4", "KO1"), "Aicda-null 2": ("jh4", "KO2"), "Aicda-null 3": ("jh4", "KO3"),
            "C1": ("n114", "C1"), "D3": ("n114", "D3"), "D5": ("n114", "D5"),
            "C3": ("n114", "C3"), "C4": ("n114", "C4"), "D1": ("n114", "D1")}
SANGER_ROWS = {"Wt 2": "Wt 2", "Wt 3": "Wt 3", "Wt 4": "Wt 4", "Aicda-null 1": "Aicda−/− 1",
               "Aicda-null 2": "Aicda−/− 2", "Aicda-null 3": "Aicda−/− 3"}

COVER = {
    0: f"{date.today():%B} {date.today().day}th, {date.today().year}",
    2: f"Please consider our manuscript by Abdulmateen Aderinto & Alberto Martin, entitled \"{TITLE},\" for "
       "publication in the Journal of Immunological Methods.",
    5: "We therefore developed a cloning-free nanopore workflow for SHM, implemented it as SHMtrace, an "
       "open-source analysis tool that takes the reads as delivered by the sequencing provider, removes "
       "off-target, chimeric and low-quality reads, and tests every candidate substitution against matched "
       "negative-control libraries with false-discovery control. We benchmarked it directly against "
       "cloned-amplicon Sanger sequencing in the same animals, re-analysed from the original chromatograms.",
    6: "The two platforms agreed. Both recovered the expected AID-dependent difference between wild-type and "
       "Aicda−/− mice (at least 137-fold by nanopore, ~38-fold by Sanger), and their pooled mutation spectra "
       "correlated closely across the twelve substitution types (Pearson r = 0.93). An independent monoclonal "
       "system, an AID-expressing hybridoma line, gave at least a 92-fold difference over empty-vector "
       "controls, none of which had a call, with the expected G:C-biased signature, indicating that the "
       "framework is not specific to a single locus. Residual background in AID-deficient controls was "
       "~0.5 × 10−5 per analysed base.",
    7: "Nanopore also yields information the cloning-based assay cannot. Because each read spans the entire "
       "amplicon, we could analyse ~1,140 bp rather than the ~530 bp window accessible to Sanger sequencing, "
       "and retain co-occurring mutations linked on individual reads, from which clonal mutation trees can be "
       "drawn.",
}
HIGHLIGHTS = {
    1: "Control-calibrated ONT calling profiles full-length SHM amplicons without cloning",
    2: "Matched biological controls define a Q53-equivalent residual background",
    3: "ONT reproduces the paired Sanger JH4 substitution spectrum (r = 0.93)",
    4: "At 500 reads, variant discovery reaches 91–94% of its maximum in this dataset",
}


def set_paragraph(p, text):
    runs = p.runs
    runs[0].text = text
    for r in runs[1:]:
        r._r.getparent().remove(r._r)


def to_pdf(docx, outdir):
    subprocess.run([SOFFICE, "--headless", "--convert-to", "pdf", str(docx), "--outdir", str(outdir)],
                   check=True, capture_output=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def accounting():
    """SHMtrace's per-library read accounting (table_s1.json) for the manuscript runs."""
    return {run: json.load(open(RUNS_DIR / f"{run}_shmtrace/table_s1.json")) for run in ("jh4", "n114")}


def build_deposit():
    if DEP.exists():
        shutil.rmtree(DEP)
    shutil.copytree(DEP10, DEP, ignore=shutil.ignore_patterns(".DS_Store"))
    raw = DEP / "Raw_Data"
    man = raw / "04_Manifests_and_Integrity"
    out = raw / "05_SHMtrace_Analysis"
    shutil.rmtree(out)
    for src, name in RUNS.items():
        shutil.copytree(RUNS_DIR / src, out / name, ignore=shutil.ignore_patterns(".DS_Store"))
    for src, name in SUBSET_RUNS.items():
        shutil.copytree(RUNS_DIR / src, out / "Comparison_preprocessed_subsets" / name,
                        ignore=shutil.ignore_patterns(".DS_Store"))
    shutil.copytree(HERE / "jh4_calls", out / "JH4_529bp_window/calls_by_sample")
    shutil.copytree(HERE / "n114_calls", out / "N114_V_region/calls_by_sample")
    refs = out / "SHMtrace_references"
    refs.mkdir()
    for f in ("mouse_jh4_intron.fasta", "mouse_jh4_contaminant.fasta", "n114_vregion.fasta"):
        shutil.copy(REPO / "src/shmtrace/references" / f, refs / f)
    for f in ("jh4_intron.yaml", "n114_vregion.yaml"):
        shutil.copy(REPO / "src/shmtrace/presets" / f, refs / f)
    # read accounting, one row per ONT library
    acct = accounting()
    keys = ["delivered", "too_long", "low_quality", "not_target", "chimeric", "short_alignment",
            "paralogue_alleles", "usable", "analysed", "candidates_tested", "calls"]
    with open(out / "READ_ACCOUNTING.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assay", "library"] + keys)
        for run, libs in acct.items():
            for lib, a in libs.items():
                w.writerow([run.upper(), lib] + [a[k] for k in keys])
    (out / "README.md").write_text(DEPOSIT_ANALYSIS_README)

    # sample manifest: SHMtrace's usable/analysed counts; the earlier preprocessing column is renamed
    sanger = {k: v for k, v in json.load(open(HERE / "sanger_calls/corrected_sanger.json")).items()
              if isinstance(v, dict)}
    rows = list(csv.DictReader(open(man / "SAMPLE_MANIFEST.csv", encoding="utf-8")))
    head = [c for c in rows[0] if c not in ("records_passing_preprocessing", "primary_analysis_depth", "notes")]
    head += ["shmtrace_usable", "shmtrace_analysed", "earlier_preprocessing_passing", "notes"]
    with open(man / "SAMPLE_MANIFEST.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, head)
        w.writeheader()
        for r in rows:
            name = r["manuscript_sample"]
            if r["assay"].startswith("ONT"):
                run, lib = ONT_ROWS[name]
                a = acct[run][lib]
                usable, analysed = a["usable"], a["analysed"]
            else:
                usable = analysed = sanger[SANGER_ROWS[name]]["clones"]
            w.writerow({**{c: r[c] for c in head if c in r}, "shmtrace_usable": usable,
                        "shmtrace_analysed": analysed,
                        "earlier_preprocessing_passing": r["records_passing_preprocessing"],
                        "notes": r["notes"]})
    (man / "VER11_NOTES.md").write_text(VER11_NOTES)
    (raw / "README.md").write_text(RAW_README)
    (DEP / "README.md").write_text(DEPOSIT_README)
    (man / "VER10_ADDITIONS_MANIFEST.csv").unlink()      # listed the Ver10 analysis folder, now replaced
    added = sorted(p for p in out.rglob("*") if p.is_file())
    with open(man / "VER11_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in added:
            w.writerow([str(p.relative_to(raw)), p.stat().st_size, sha256(p)])
    skip = {"SHA256SUMS.txt", "FILE_MANIFEST.csv"}
    files = sorted(p for p in raw.rglob("*") if p.is_file() and not (p.parent == man and p.name in skip))
    digests = {p: sha256(p) for p in files}
    with open(man / "FILE_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["relative_path", "size_bytes", "sha256", "top_level_category"])
        for p in files:
            rel = p.relative_to(raw)
            w.writerow([str(rel), p.stat().st_size, digests[p], rel.parts[0] if len(rel.parts) > 1 else ""])
    (man / "SHA256SUMS.txt").write_text("".join(f"{digests[p]}  {p.relative_to(raw)}\n" for p in files))
    print(f"deposit: {DEP.name}, {len(files)} checksummed files, {len(added)} in 05_SHMtrace_Analysis")


def build_package():
    if PKG.exists():
        shutil.rmtree(PKG)
    for d in ("1_Manuscript", "2_Cover_Letter", "3_Highlights", "4_Figures/PNG", "4_Figures/TIFF",
              "4_Figures/Vector_PDF", "4_Figures/Editable_SciCanvas", "5_Graphical_Abstract", "6_Code", "7_Data",
              "8_Supporting"):
        (PKG / d).mkdir(parents=True)
    ms = PKG / "1_Manuscript/ONT_SHM_Manuscript_JIM_Ver21.docx"
    shutil.copy(MS, ms)
    to_pdf(ms, ms.parent)
    (ms.parent / "README.md").write_text(
        "# Manuscript\n\n`ONT_SHM_Manuscript_JIM_Ver21.docx` and its PDF. Passages changed are highlighted "
        "yellow; remove the highlighting before submission. See `8_Supporting/CHANGES_Ver21.md`.\n")

    cover = Document(V13 / "2_Cover_Letter/ONT_SHM_CoverLetter_JIM.docx")
    for i, text in COVER.items():
        set_paragraph(cover.paragraphs[i], text)
    cl = PKG / "2_Cover_Letter/ONT_SHM_CoverLetter_JIM.docx"
    cover.save(cl)
    to_pdf(cl, cl.parent)
    hl = Document(V13 / "3_Highlights/ONT_SHM_Highlights.docx")
    for i, text in HIGHLIGHTS.items():
        if len(text) > 85:
            raise RuntimeError(f"highlight over 85 characters: {text}")
        set_paragraph(hl.paragraphs[i], text)
    hp = PKG / "3_Highlights/ONT_SHM_Highlights.docx"
    hl.save(hp)
    to_pdf(hp, hp.parent)

    rows = []
    for name, stem in FIGURES.items():
        sci = stem.parent == SC
        png = stem.parent / f"{stem.name}_preview.png" if sci else stem.with_suffix(".png")
        shutil.copy(png, PKG / f"4_Figures/PNG/{name}.png")
        if sci:          # line art: SciCanvas's own 1,000-dpi render, and the editable files
            shutil.copy(stem.parent / f"{stem.name}_1000dpi.tif", PKG / f"4_Figures/TIFF/{name}.tif")
            for ext in (".svg", ".scifig"):
                shutil.copy(stem.with_suffix(ext), PKG / f"4_Figures/Editable_SciCanvas/{name}{ext}")
            dpi = 1000
        else:
            im = Image.open(png).convert("RGB")
            im.save(PKG / f"4_Figures/TIFF/{name}.tif", dpi=(600, 600), compression="tiff_lzw")
            dpi = round(im.info.get("dpi", (600, 600))[0])
        pdf = stem.with_suffix(".pdf")
        if pdf.exists():
            shutil.copy(pdf, PKG / f"4_Figures/Vector_PDF/{name}.pdf")
        w, h = Image.open(png).size
        rows.append(f"| {name} | `{png.relative_to(REPO)}` | {w} × {h} | {dpi} |")
    (PKG / "4_Figures/README.md").write_text(
        "# Figures\n\nPNG, LZW-compressed TIFF and vector PDF for every figure in Ver21. Figure 1 and "
        "Supplementary Figure S11 are line art drawn in SciCanvas: their TIFFs are 1,000 dpi at 174 mm, and the "
        "editable SVG and .scifig files are in `Editable_SciCanvas/`. The other TIFFs are 600 dpi. Each figure is "
        "regenerated by the script named in `6_Code/README.md` from the deposited data.\n\n"
        "| Figure | Source in the SHMtrace repository | Pixels | TIFF dpi |\n|---|---|---|---|\n"
        + "\n".join(rows) + "\n")
    shutil.copytree(FIG / "graphical_abstract", PKG / "5_Graphical_Abstract", dirs_exist_ok=True)

    code = PKG / "6_Code"
    shutil.copytree(HERE / "scripts", code / "figure_scripts",
                    ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
    sc = code / "figure_scripts/figure1_scicanvas"
    sc.mkdir()
    for f in ("build.js", "figure1_overview.js", "figure1_detailed.js"):
        shutil.copy(HERE / "figure1_scicanvas" / f, sc)
    (code / "analysis_scripts").mkdir()
    for f in ("jh4_shmtrace.py", "n114_shmtrace.py", "read_accounting.py", "depth_shmtrace.py",
              "sanger_vs_nanopore.py", "sanger_window_figure.py", "mutation_tree_figure.py", "organize_jh4_ab1.py",
              "sanger_lineages.py", "samplesheets/make_jh4_samplesheets.py"):
        shutil.copy(CAL / f, code / "analysis_scripts")
    for f in ("manuscript_stats.py", "export_sanger_calls.py", "revise_manuscript_sanger.py",
              "revise_manuscript_v18.py", "revise_manuscript_v19.py", "revise_manuscript_v20.py",
              "revise_manuscript_v21.py", "build_submission_v21.py"):
        shutil.copy(HERE / f, code / "analysis_scripts")
    src = code / "SHMtrace_v0.2_source.zip"
    subprocess.run(["git", "-C", str(REPO), "archive", "--format=zip", "-o", str(src), "HEAD",
                    "src", "pyproject.toml", "README.md", "LICENSE", "tests"], check=True)
    (code / "README.md").write_text(CODE_README)

    data = PKG / "7_Data"
    for src_run, name in RUNS.items():
        d = RUNS_DIR / src_run
        for f in ("samples.csv", "variants.csv", "animal_tests.csv", "table_s1.json", "stats.json",
                  "jh4.json", "n114.json", "trees.json"):
            if (d / f).exists():
                shutil.copy(d / f, data / f"{name}__{f}")
    shutil.copy(DEP / "Raw_Data/05_SHMtrace_Analysis/READ_ACCOUNTING.csv", data)
    shutil.copy(FIG / "Individual_mutation_matrices_counts.csv", data)
    (data / "README.md").write_text(
        "# Data underlying the figures\n\nPer-run SHMtrace tables (`<run>__samples.csv`, `__variants.csv`, "
        "`__animal_tests.csv`, `__table_s1.json`) for the JH4 529-bp window, the JH4 full-length interval, N114, "
        "the depth series, the Sanger–ONT comparison and the mutation trees, all from the reads as delivered. "
        "`READ_ACCOUNTING.csv` gives what SHMtrace did with every delivered read (Supplementary Table S1 and "
        "Figure S11E). The complete project files, the reads and the earlier subset runs are in the data "
        "deposit (`Raw_Data/05_SHMtrace_Analysis/`).\n")
    for f in ("CHANGES_Ver21.md", "CHANGES_Ver18.md"):
        shutil.copy(HERE / f, PKG / "8_Supporting" / f)
    (PKG / "README.md").write_text(PACKAGE_README)
    print(f"package: {PKG.name}")


PACKAGE_README = f"""# Submission package, Ver21 ({STAMP})

**{TITLE}** (Aderinto and Martin). Target: Journal of Immunological Methods.

| Folder | Contents |
|---|---|
| `1_Manuscript/` | Ver21 manuscript (docx, PDF); changed passages highlighted |
| `2_Cover_Letter/` | Cover letter, updated to the Ver21 results |
| `3_Highlights/` | Highlights (≤ 85 characters each) |
| `4_Figures/` | Figures 1–4 and S1–S11: PNG, TIFF, vector PDF; editable SciCanvas files for Figure 1 and S11 |
| `5_Graphical_Abstract/` | Graphical abstract (PNG, TIFF, PDF) |
| `6_Code/` | SHMtrace v0.2 source, the analysis drivers, and the figure scripts |
| `7_Data/` | The tables behind every figure, and SHMtrace's read accounting |
| `8_Supporting/` | `CHANGES_Ver21.md` (Ver19–21) and `CHANGES_Ver18.md` |

Every ONT result is from SHMtrace run end to end on the reads as delivered by the sequencing provider.
Reads, chromatograms and full analysis outputs: `ONT_SHM_DataDeposit_Ver11_{STAMP}`.
"""

DEPOSIT_README = f"""# Raw and intermediate sequencing data, Ver11 ({STAMP})

Companion data deposit for **"{TITLE}"** (Aderinto and Martin, Department of Immunology,
University of Toronto). Intended for SRA/BioProject and Zenodo; the resulting accessions
should be cited in the manuscript's Data and code availability section.

## Contents

| Directory | Contents |
|---|---|
| `Raw_Data/01_JH4_Sanger/` | Sanger clone records (GenBank), the original `.ab1` chromatograms of the six manuscript mice, the clone sample sheet, and SHMtrace's Sanger calls |
| `Raw_Data/02_JH4_ONT/` | J<sub>H</sub>4 ONT FASTQs as delivered (`Original_FASTQ/`, the analysis input), and the earlier preprocessing script's subsets and Geneious exports, kept for comparison |
| `Raw_Data/03_N114_ONT/` | N114 V-region equivalents |
| `Raw_Data/04_Manifests_and_Integrity/` | Sample and file manifests, validation, `VER11_NOTES.md`, `VER11_ADDITIONS_MANIFEST.csv`, SHA-256 checksums |
| `Raw_Data/05_SHMtrace_Analysis/` | Every SHMtrace analysis behind the Ver21 figures and tables, its references and settings, and the per-library read accounting (see its README) |

## What changed from Ver10

- SHMtrace now analyses the delivered FASTQs end to end: length and expected-error gates, competitive
  alignment against paralogue references, the paralogue-allele filter, a random draw of 500 usable reads
  (seed 0), and calling against the negative controls at a false-discovery rate of 0.05. Ver10's analyses,
  made on the earlier preprocessing script's 500-read subsets, are kept in
  `05_SHMtrace_Analysis/Comparison_preprocessed_subsets/`.
- `SAMPLE_MANIFEST.csv` gives SHMtrace's usable and analysed counts per library (and the Sanger clones
  analysed); the earlier preprocessing's counts are kept as `earlier_preprocessing_passing`.
- `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over every file under `Raw_Data/`.

## Integrity

```bash
cd Raw_Data && shasum -a 256 -c 04_Manifests_and_Integrity/SHA256SUMS.txt
```
"""

RAW_README = """# Raw and intermediate sequencing data

The sequencing data and analyses that underpin the submitted manuscript, organised by assay.

## Directory map

| Directory | Contents |
|---|---|
| `01_JH4_Sanger/Original_Chromatograms_AB1/` | The `.ab1` chromatograms of the six manuscript mice: the Sanger analysis input. |
| `01_JH4_Sanger/Original_Clone_Sequences/` | The hand-collated GenBank clone records, grouped as originally filed. |
| `01_JH4_Sanger/SHMtrace_Sanger_calls/` | SHMtrace's Sanger-mode calls (base quality ≥ 20) per mouse. |
| `01_JH4_Sanger/Alignment_and_Variant_Exports/` | Earlier Geneious alignments and variant exports, kept for comparison. |
| `02_JH4_ONT/Original_FASTQ/` | ONT FASTQ files as delivered for the three Wt and three *Aicda*−/− JH4 libraries: the analysis input. |
| `02_JH4_ONT/Reference_Sequences/` | References used by the earlier preprocessing; SHMtrace's are in `05_SHMtrace_Analysis/SHMtrace_references/`. |
| `02_JH4_ONT/Processed_FASTQ_Depth_Subsets/`, `Alignment_and_Variant_Exports/`, `Pipeline_Scripts_and_Summaries/` | The earlier preprocessing script, its filtered read sets and depth subsets, and Geneious exports made from them. Not used for the reported results; kept for comparison. |
| `03_N114_ONT/` | The same for the three pMX-Aicda and three empty-vector N114 subclones. |
| `04_Manifests_and_Integrity/` | Sample and file manifests, source scope, validation and SHA-256 checksums. |
| `05_SHMtrace_Analysis/` | Every reported analysis (SHMtrace v0.2), with references, settings and read accounting. |

## How the reported results were made

SHMtrace v0.2 read each delivered FASTQ and (1) removed reads longer than 2,000 bp (JH4) or 900 bp (N114)
and reads averaging more than 0.01 expected errors per base; (2) aligned the rest with minimap2 (map-ont)
to the target and its competing references, removing reads aligning best elsewhere, chimeric reads and
alignments shorter than 800 bp (JH4) or 400 bp (N114); (3) removed JH4 reads carrying ten or more
paralogue-specific bases; (4) drew 500 usable reads per library at random (seed 0); and (5) tested every
candidate substitution against the error rate measured in the negative controls (Aicda−/− for JH4, empty
vector for N114) at a Benjamini–Hochberg false-discovery rate of 0.05, removing sites at ≥ 50%. A unique
variant is one substitution at one reference position, counted once per sample. Windows: JH4 positions
101–629 (529 bp; the Sanger comparison) or 90–1,230 (full length); N114 positions 21–635 (615 bp).

## Scope

Only libraries that underpin the manuscript are included. Files from separate MKLN1 and MAEA experiments,
archival Sanger folders and alternate N114 libraries were not analysed and are not included; see
`04_Manifests_and_Integrity/SOURCE_SCOPE_AND_EXCLUSIONS.md`. `PACKAGING_AUDIT.md` describes the original
(July 2026) packaging; `VER11_NOTES.md` what changed since.

## Integrity

`04_Manifests_and_Integrity/FILE_MANIFEST.csv` lists every file with its size and SHA-256 digest, and
`SHA256SUMS.txt` gives the same digests for `shasum -a 256 -c`. The two lists omit themselves. Raw FASTQ,
chromatograms and GenBank records are included without modification.
"""

VER11_NOTES = f"""# Ver11 notes ({STAMP})

- The reported ONT analyses are SHMtrace runs on `Original_FASTQ/` (the reads as delivered), not on
  `Processed_FASTQ_Depth_Subsets/`. The subsets, the earlier preprocessing scripts and the Geneious
  exports remain for comparison; `PACKAGING_AUDIT.md` (July 2026) describes them as they were packaged.
- `SAMPLE_MANIFEST.csv`: `shmtrace_usable` and `shmtrace_analysed` are SHMtrace's counts (Supplementary
  Table S1); for Sanger they are the clones analysed from the chromatograms. `earlier_preprocessing_passing`
  is the former `records_passing_preprocessing` column (the earlier script's counts, and for Sanger the
  GenBank records as originally filed).
- `05_SHMtrace_Analysis/READ_ACCOUNTING.csv`: the fate of every delivered read per library.
- `VER10_ADDITIONS_MANIFEST.csv` is replaced by `VER11_ADDITIONS_MANIFEST.csv` (the analysis folder).
"""

DEPOSIT_ANALYSIS_README = """# SHMtrace analyses (SHMtrace v0.2)

Each folder is one SHMtrace project: `project.json` (every result), `report.html`, CSV tables and
figures; `table_s1.json` gives SHMtrace's read accounting and the candidates tested per library.
All manuscript analyses start from the reads as delivered (`0x_*/Original_FASTQ/`).

| Folder | Analysis | Used for |
|---|---|---|
| `JH4_529bp_window/` | JH4 positions 101–629; calls per sample in `calls_by_sample/` | Figures 2A, 3; S3; Tables S1, S2 |
| `JH4_full_length_90-1230/` | JH4 positions 90–1,230, otherwise identical | Figure 2C; S1 |
| `N114_V_region/` | N114 positions 21–635; calls per sample in `calls_by_sample/` | Figure 4; S4, S5 |
| `Depth_series/` | one random draw of usable reads per library at 30–1,000 reads | S6, S7 |
| `Sanger_vs_ONT/` | the chromatograms and the JH4 calls | Figure 3 statistics, S9 |
| `Mutation_trees/` | full-length reads grouped by V(D)J junction and ordered into trees | S10 |
| `SHMtrace_references/` | target and competing references, and the presets (settings) used | all |
| `READ_ACCOUNTING.csv` | per library: delivered, removed at each step, usable, analysed, candidates, calls | Table S1; S11E |
| `Comparison_preprocessed_subsets/` | the same analyses on the earlier script's 500-read subsets | not reported; comparison only |

Calling: recommended mode (`error_model: qcal`), negative controls Aicda−/− (JH4) and empty vector
(N114); each control is called against the other two. `calibration` in `project.json` records the
leave-one-control-out factor (5 for every run) and its trace. Reproduce with
`6_Code/analysis_scripts/` of the submission package.
"""

CODE_README = """# Code

- `SHMtrace_v0.2_source.zip`: the analysis software (Python package `shmtrace`), with its tests.
  Install with `pip install .` from the unpacked folder (Python ≥ 3.11).
- `analysis_scripts/`: the drivers that produced every result from the data deposit's delivered reads.
  - `jh4_shmtrace.py` (and `--full-length`), `n114_shmtrace.py`: the manuscript runs; `--input subsets`
    repeats them on the earlier 500-read subsets. `read_accounting.py` tallies each read's fate.
  - `depth_shmtrace.py`; `mutation_tree_figure.py` (Supplementary Figure S10)
  - `sanger_vs_nanopore.py --ont-mode recommended --sanger-dir <AB1 folders> --ont-dir <JH4 Original_FASTQ>`
  - `organize_jh4_ab1.py --original`, `make_jh4_samplesheets.py`, `export_sanger_calls.py`,
    `sanger_lineages.py`, `sanger_window_figure.py`: Sanger clones
  - `manuscript_stats.py`: every statistic quoted in the text
  - `revise_manuscript_*.py`: the manuscript revisions Ver16 → Ver21; `build_submission_v21.py`: this package
- `figure_scripts/`: one script per figure; `config.py` holds the paths. `figure1_scicanvas/` draws
  Figure 1 and Supplementary Figure S11 with SciCanvas (`node build.js figure1_overview.js out/Figure1_overview`;
  needs a SciCanvas checkout, `SCICANVAS=<path>`, default `~/SciCanvas`).

Paths at the top of each driver point at the data deposit; edit them to your copy.
"""

if __name__ == "__main__":
    build_deposit()
    build_package()
