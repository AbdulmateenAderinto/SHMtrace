# Reproducing the Geneious calls

Run 11 October 2026 on the six bundled J꜀4 libraries (`examples_compare/`, 500 reads each, run DHJGKB). Reproduce with:

```bash
shmtrace-benchmark --sheet examples_compare/samples.csv \
  --reference examples_compare/JH4_reference.fasta \
  --contaminant examples_compare/JH4_contaminant.fasta \
  --preset jh4_intron --out bench
```

## The Geneious rule

The Geneious manual describes its approximate variant p-value as the probability of sequencing error producing the variant bases, computed by averaging the quality of the bases that carry the variant and taking a binomial tail. Its worked example (three Q20 reads, two carrying the variant) gives 3C2 × 0.01² × 0.99 + 0.01³ ≈ 0.0003. So, for a candidate seen in k of n reads whose variant bases average quality Q:

    p = P(X ≥ k),  X ~ Binomial(n, 10^(−Q/10))

The per-base error is **not** divided by three. The manuscript kept candidates with p ≤ 0.05, uncorrected.

A practical consequence: at 500 reads, a variant seen in one read passes only if that base is about Q40 or better. Geneious's filter is therefore mostly a filter on the base quality of single-read calls.

`--mode geneious` implements exactly this, with phasing off and no multiple-testing correction.

## Agreement with the manuscript (Table S1)

| Library | Geneious calls | `--mode geneious` | Geneious density ×10⁻⁵ | SHMtrace density ×10⁻⁵ |
|---|---|---|---|---|
| Wt2 | 161 | 172 | 60.9 | 65.5 |
| Wt3 | 342 | 344 | 129.3 | 132.9 |
| Wt4 | 68 | 69 | 25.7 | 26.4 |
| KO1 | 2 | 3 | 0.76 | 1.14 |
| KO2 | 6 | 6 | 2.27 | 2.29 |
| KO3 | 3 | 3 | 1.13 | 1.15 |

Every library is within 11 calls, and five of six are within 2. SHMtrace used 497 reads per library here, the common usable depth after alignment, against 500 in Geneious. The residual differences are most likely due to the mapper (minimap2 against the Geneious mapper) and Geneious's homopolymer quality smoothing, which is not reproduced. Neither is tested.

Alternatives that were tried and rejected: dividing the error by three (KO 7–11 calls), using the mean quality of the non-variant bases or of the whole column (KO 163–178 calls), and a summed-quality convolution (KO 31–50 calls).

## Recommended mode against Geneious

| Library | Geneious | `recommended` (BH + leave-one-out KO null) | sites shared with `--mode geneious` |
|---|---|---|---|
| Wt2 | 161 | 72 | 63 of 72 |
| Wt3 | 342 | 297 | 295 of 297 |
| Wt4 | 68 | 68 | 68 of 68 |
| KO1 | 2 | 0 | 0 |
| KO2 | 6 | 3 | 3 of 3 |
| KO3 | 3 | 1 | 1 of 1 |

The recommended mode is almost a strict subset of the Geneious calls. It gives a lower background (0.5 against 1.5 × 10⁻⁵) and a higher WT/KO separation (110-fold against 49-fold). Nearly all of the disagreement is in Wt2.

## Site-level concordance against the paper's Geneious exports

Run on 2026-10-11 against the per-sample Geneious annotation exports in the data deposit (`deposit_Ver9`), with the export re-filtered to p ≤ 0.05, frequency < 50% and the analysis window. Sheets: `benchmark_results/jh4_sheet.csv`, `n114_sheet.csv`. Full tables: `benchmark_results/site_level_jh4/`, `site_level_n114/`.

**JH4 intron (497 reads, `--mode geneious`)**

| Library | SHMtrace | Geneious | Precision | Recall | Jaccard |
|---|---|---|---|---|---|
| Wt2 | 172 | 161 | 0.91 | 0.97 | 0.88 |
| Wt3 | 344 | 342 | 0.98 | 0.98 | 0.96 |
| Wt4 | 69 | 68 | 0.96 | 0.97 | 0.93 |
| KO1 | 3 | 2 | 0.67 | 1.00 | 0.67 |
| KO2 | 6 | 6 | 0.83 | 0.83 | 0.71 |
| KO3 | 3 | 3 | 1.00 | 1.00 | 1.00 |

Recommended mode recalls 0.38 (Wt2), 0.86 (Wt3) and 0.97 (Wt4) of the Geneious sites with precision 0.85–0.99.

**N114 V region (494 reads, 655-nt Geneious reference `seq_23353`, `--mode geneious`)**

| Library | SHMtrace | Geneious | Precision | Recall | Jaccard |
|---|---|---|---|---|---|
| C1 AID | 72 | 80 | 0.67 | 0.60 | 0.46 |
| D3 AID | 81 | 71 | 0.64 | 0.73 | 0.52 |
| D5 AID | 150 | 147 | 0.71 | 0.72 | 0.56 |
| C3 empty | 2 | 2 | – | – | 0.00 |
| C4 empty | 5 | 7 | 0.20 | 0.14 | 0.09 |
| D1 empty | 3 | 10 | 1.00 | 0.30 | 0.30 |

Counts and AID/empty separation match, but individual sites agree less than for JH4. **Cause, checked 2026-10-11: the deposited N114 sub500 FASTQs are not the reads Geneious analysed.**

- The Geneious assembly BAMs (`Submission Package - JIM/6_Raw Data/N114_Vregion_ONT/*sub500 assembled to seq_23353.bam`) share only 80 of 500 read names with the deposited C1 sub500 FASTQ and 151 of 500 with D5. That is what two independent random 500-read draws from the same library would share. For JH4 the deposited Wt3 sub500 FASTQ and the Geneious BAM share all 500 reads.
- The C1 Geneious BAM pileup reproduces the exported read counts (86 of 92 variant rows exact), so the BAM is the assembly the export came from.
- Running SHMtrace `--mode geneious` on the reads extracted from the Geneious BAMs gives the JH4-level agreement:

| Library | SHMtrace | Geneious | Shared | Precision | Recall | Jaccard |
|---|---|---|---|---|---|---|
| C1 AID (Geneious's reads) | 80 | 80 | 74 | 0.93 | 0.93 | 0.86 |
| D5 AID (Geneious's reads) | 148 | 147 | 138 | 0.93 | 0.94 | 0.88 |

So SHMtrace's Geneious emulation holds on N114 too. The deposit itself has a reproducibility gap: re-running Geneious on the deposited N114 sub500 FASTQs will not reproduce the deposited N114 annotation tables. Either re-deposit the reads Geneious used (extractable from the BAMs) or regenerate the exports from the deposited FASTQs.

## Still open

- Replace the N114 sub500 FASTQs in the data deposit with the reads in the Geneious BAMs (or regenerate the exports), then re-run `benchmark_results/n114_sheet.csv`.
- The remaining 5–10% site disagreement on both amplicons is single-read calls near Q40, where one read's quality decides the call.
