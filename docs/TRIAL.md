# Trying SHMtrace on the bundled example

A ten-minute check that the install works and produces the expected numbers.
[`QUICKSTART.md`](QUICKSTART.md) is the one-page guide to using it on your own
data; [`USER_GUIDE.md`](USER_GUIDE.md) is the full account.

## 1. One-time setup

```bash
cd ~/Nanoshm/nanoshm_v0.2
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[app,bam,sanger]"
```

Next time, only `cd` and `source .venv/bin/activate`.

## 2. The app

```bash
shmtrace-app
```

A window opens. If no window appears, run `shmtrace-app --no-window` and open
the address it prints.

In **1 · Region and data type** leave the defaults (`jh4_intron`, `ont`,
`recommended`). In **2 · Reference** choose the built-in **Jh4 – Base**; the
matching contaminants are added for you. In **3 · Samples** drag in the six
FASTQ files from `examples_compare/`, then press **Controls by name**, which
marks the three KO libraries as controls. Name the project and press
**Analyse**.

It takes about a minute. Expect:

| | |
|---|---|
| Wild type | ~57 × 10⁻⁵ unique mutations per base |
| *Aicda*-/- | ~0.6 × 10⁻⁵ |
| Separation | ~100-fold, P ≈ 0.002 (Welch on log densities, n = 3 vs 3) |

The **Results** tab draws the density, frequency, spectrum, positional map and
hotspot figures. **Tables** lists every call. **Export** writes CSV, Excel,
figures (PNG, SVG, PDF), a standalone HTML report, or a ZIP of all of it.

Per-sample numbers to expect, at the common depth of 497 reads:

| Sample | Unique mutations | Density ×10⁻⁵ |
|---|---:|---:|
| WT2 | 146 | 55.7 |
| WT3 | 334 | 129.3 |
| WT4 | 68 | 26.0 |
| KO1 | 0 | 0.0 |
| KO2 | 3 | 1.1 |
| KO3 | 1 | 0.4 |

## 3. The command line

The same project from a sample sheet:

```bash
shmtrace-benchmark --sheet examples_compare/samples.csv \
  --reference examples_compare/JH4_reference.fasta \
  --contaminant examples_compare/JH4_contaminant.fasta \
  --preset jh4_intron --out bench
```

One sample on its own, with the manuscript's Geneious rule for comparison:

```bash
SHMtrace --fastq examples/WT3_JH4_sub500.fastq \
  --reference examples/JH4_reference.fasta \
  --contaminant examples/JH4_contaminant.fasta \
  --preset jh4_intron --mode geneious --out results/WT3 --plots
```

## 4. Your own data

Swap in your FASTQ files and a reference FASTA for your amplicon. Three things
matter more than any setting:

- **Mark your controls.** A mutation-free sample from the same flow cell (an
  AID-knockout animal, an empty-vector line, germline DNA) is what measures the
  sequencing background. Two or more is better than one; with none, every call
  rests on base quality alone and the background rises.
- **Check the window.** The preset's analysis window must lie inside the region
  your reads actually cover.
- **Read the warnings column.** Low depth, a missing control, and calls that
  cluster in a few bases are all reported there.

`docs/calibration/` records how the defaults were chosen and how they behave on
independent data.
