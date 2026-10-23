# SHMtrace demo dataset: J_H4 intron in wild-type and *Aicda*⁻/⁻ mice

This is real data from the SHMtrace paper. The same six mice were read two ways: nanopore reads of the J_H4 intron,
and Sanger clones of the same amplicon. Running it reproduces the paper's headline J_H4 numbers exactly. It takes about
half a minute on a laptop.

## What it shows

| SHMtrace feature | Where you see it | Paper |
|---|---|---|
| Read preparation: length and error gates, competitive alignment against the paralog, matched depth | `qc_waterfall.png`, `samples.csv` | Figure S1E |
| Control-calibrated calling: AID-null libraries as the error model, with an empirical decoy threshold | calls in `variants.csv`: hundreds per Wt mouse, 0–2 per control | Figures 2 and 4 |
| Animal-level statistics, one value per mouse | `animal_tests.csv`, `density_by_group.png` | Figure 4 |
| AID signature: spectrum and WRC/GYW hotspot use | `spectrum_by_group.png`, `hotspot_targeting.png` | Figure 5 |
| Whole-amplicon positional profile | `positional_map.png` | Figure 6A |
| Clones from the V(D)J junction, and read- vs clone-weighted frequency | `clones.csv`, `clone_map.png` | Figure 6C |
| Sanger route: Q≥20 counting, background reported beside the counts | `sanger/` run outputs | Figure 5 |
| Sanger bookkeeping checks (flags, never silent changes) | `qc_flags.csv`, `excluded.csv` | STAR Methods |
| Agreement between platforms in the same mice | `shmtrace-concordance` | Figure 5 |

## Contents

```
nanopore/   WT2, WT3, WT4 (wild type) and KO1, KO2, KO3 (Aicda-/-): every read delivered (4,558-5,000 each), gzipped
            samples.csv   the sample sheet: file, sample, group, control
sanger/     traces/       184 .ab1 chromatograms filed to the same six mice, in their sequencing-order folders
            samples.csv   the sample sheet the paper used, including two rows excluded with their reasons
reference/  mouse_jh4_intron.fasta and mouse_jh4_contaminant.fasta (the paralog used for competitive alignment)
expected/   the summary tables a correct run produces, to compare against
```

## Run it

**Install SHMtrace** first; see the main README. Then, from this folder:

```bash
# nanopore: about 20 s
shmtrace-sheet run nanopore/samples.csv --reference reference/mouse_jh4_intron.fasta \
  --contaminant reference/mouse_jh4_contaminant.fasta --region jh4_intron -o results/nanopore

# Sanger: about 10 s
shmtrace-sheet run sanger/samples.csv --reference reference/mouse_jh4_intron.fasta \
  --contaminant reference/mouse_jh4_contaminant.fasta --region jh4_intron --platform sanger -o results/sanger

# do the two platforms agree in the same mice?
shmtrace-concordance results/sanger/project.json results/nanopore/project.json \
  --pair WT2=WT2 --pair WT3=WT3 --pair WT4=WT4 --labels Sanger,Nanopore
```

Open `results/nanopore/report.html` for every result on one page.

**In the desktop app** (`shmtrace-app`):
1. Choose region `jh4_intron`, platform `ont` and mode `recommended`.
2. Choose the built-in J_H4 reference. Its paralog is added automatically.
3. Drop in the six files from `nanopore/` and press **Controls by name**, so that KO1–KO3 become controls.
4. Press **Analyse**.

## What you should get

Nanopore, at 500 reads per mouse:

| Mouse | Group | Calls | Unique variants per 10⁵ bases |
|---|---|---:|---:|
| WT2 | Wt | 221 | 83.7 |
| WT3 | Wt | 344 | 133.2 |
| WT4 | Wt | 73 | 27.8 |
| KO1 | *Aicda*⁻/⁻ | 0 | 0 |
| KO2 | *Aicda*⁻/⁻ | 2 | 0.76 |
| KO3 | *Aicda*⁻/⁻ | 0 | 0 |

The Wt/*Aicda*⁻/⁻ separation in unique-variant density is 141-fold (Welch test on log density, p = 0.0027, n = 3 mice per
group). Controls with no call enter the test at one event, as stated in the table.

Sanger: 134 clones are analysed: 37, 40 and 26 from the Wt mice and 11, 10 and 10 from the *Aicda*⁻/⁻ mice. They carry
135, 230 and 106 mutations in the Wt mice and 1, 1 and 2 in the controls.

The bookkeeping checks raise 33 flags without changing anything. They include:
- an identical mutated clone filed under both WT2 and WT3;
- tubes whose labels name a different mouse;
- two traces the sample sheet excludes, with its reasons.

This is what SHMtrace is for: problems are shown to you rather than silently resolved.

The two platforms give similar substitution spectra (`shmtrace-concordance`). Their absolute densities are compared
within a platform, not between: Sanger counts every Q≥20 change in a few dozen single clones, while nanopore counts
control-calibrated calls among 500 reads.

Your numbers should match the files in `expected/` exactly. The analysis is deterministic (fixed-seed sampling).

## Notes

- The S5F targeting index needs the HH_S5F mutability model, which is not bundled (AGPL licence). Supply your own copy
  if you want it. The WRC/GYW and WA/TW hotspot shares need nothing extra.
- To try your own data, replace the files and the sample sheets, and keep at least two AID-null controls from the same
  sequencing run. The main README lists what to check before believing a number.

## Source and citation

These data are from the Martin laboratory (University of Toronto). The nanopore libraries will be deposited with the
paper, and the Sanger traces are from the mice analysed in Barbulescu et al. (2025, *J. Immunol.*). If you use this
dataset, please cite the SHMtrace paper (citation to be added on publication).
