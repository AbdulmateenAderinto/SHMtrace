# How to trial SHM-webtool

Open **Terminal** (Cmd-Space → "Terminal") and paste each block.

## 1. One-time setup
```bash
cd "/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/Paper #1 - Nanopore sequencing/June 2nd analyses/SHM-webtool"
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
(Next time you only need the `cd` and `source .venv/bin/activate` lines.)

## 2. Start the web app
```bash
uvicorn webapp.main:app --reload
```
Open **http://127.0.0.1:8000** in your browser. There are two tabs at the top:
**Single sample** and **Compare groups**. Press **Ctrl-C** in Terminal to stop.

---

## Test case A — Single sample
In the **Single sample** tab, choose these files (all in the `examples/` folder):

| Field | File |
|---|---|
| Reads (FASTQ) | `examples/WT3_JH4_sub500.fastq` |
| Reference (FASTA) | `examples/JH4_reference.fasta` |
| Contaminant (FASTA) | `examples/JH4_contaminant.fasta` |
| Assay preset | `jh4_intron` |

Click **Run analysis**. Expect ~**320 unique variants**, a G:C/A:T-mixed spectrum,
and the spectrum / position-map / per-read figures.

Then try the negative control `examples/KO3_JH4_sub500.fastq` (same reference /
preset) → should drop to ~**25–30 variants**. That WT-vs-KO contrast is the
AID-dependence signal.

## Test case B — Compare groups (the SHMTool-style analysis)
In the **Compare groups** tab, files are in the `examples_compare/` folder:

| Field | Files (select multiple) |
|---|---|
| CONTROL reads | `WT2.fastq`, `WT3.fastq`, `WT4.fastq` |
| CASE reads | `KO1.fastq`, `KO2.fastq`, `KO3.fastq` |
| Reference | `JH4_reference.fasta` |
| Contaminant | `JH4_contaminant.fasta` |
| Preset | `jh4_intron` |
| Subregion | (leave blank, or try `200-400`) |

Click **Run comparison**. Expect the wild-type control to show ~**7× higher**
overall mutation frequency than the Aicda-/- case, with **\*\*\*** significance
on overall, A:T, and G:C mutations and collapsed AID-hotspot targeting in the
knockout — plus three comparison figures.

> To select several files at once in the file picker, hold **Cmd** and click each.

---

## Command-line alternative (no browser)
```bash
# single sample
python -m nanoshm.pipeline --fastq examples/WT3_JH4_sub500.fastq \
  --reference examples/JH4_reference.fasta --contaminant examples/JH4_contaminant.fasta \
  --preset presets/jh4_intron.yaml --out results/WT3 --plots

# case vs control
python -m nanoshm.compare \
  --control examples_compare/WT2.fastq examples_compare/WT3.fastq examples_compare/WT4.fastq \
  --case    examples_compare/KO1.fastq examples_compare/KO2.fastq examples_compare/KO3.fastq \
  --reference examples_compare/JH4_reference.fasta \
  --contaminant examples_compare/JH4_contaminant.fasta \
  --preset presets/jh4_intron.yaml --plots --out compare_results
```

## Using your own data
Swap in any FASTQ + a reference FASTA for that amplicon (+ optional contaminant
FASTA). Pick the `jh4_intron` or `n114_vregion` preset, or adjust the sliders /
window. Remember: this is the uncalibrated prototype, so read the magnitudes as
illustrative until the Phase-1 validation in `ROADMAP.md` is done.
