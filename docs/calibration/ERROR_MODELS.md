# Error models: strand-aware and quality-calibrated backgrounds

*11 October 2026. Data: the paper deposit (J꜀4 intron, 3 WT and 3 Aicda-/-; N114 V region, 3 AID and 4 empty vector, using the reads in Geneious's own BAMs).*

## Why new models

A nanopore miscall and an AID mutation look the same in one read. The v0.1 default mode learns a per-position error rate from mutation-free controls and asks whether a sample has more mismatches than that. Two things it ignored:

1. **Base quality.** In the J꜀4 Aicda-/- controls, mismatches at Q40 or above sit almost entirely at five positions: the germline differences at 349, 351, 433 and 436 (about 85% of reads) and a polymorphism at 543 (16%, in KO and WT alike). Everywhere else a Q40+ mismatch occurs about once per 100,000 bases. Low-quality mismatches are a hundred times more common. A test that pools every mismatch therefore wastes the information that the sequencer itself gives about each base.
2. **Strand.** Some ONT errors are strand-specific (they depend on the k-mer that passed the pore). A site miscalled on the minus strand in controls should not hide a mutation seen on both strands, and a mutation seen on one strand only deserves suspicion.

## The models

| Mode | What it tests |
|---|---|
| `qcal` | For a candidate with k supporting reads, sorted by base quality q₁ ≥ q₂ ≥ … , and for each j ≤ k: is it surprising to see j miscalls at quality ≥ q_j? The expected count is Σ over strands of depth × r(q_j), where r is the control rate of miscalls at that quality or better for that substitution and strand, blended with the site's own control count. The call's p-value is the smallest over j, multiplied by k. |
| `strand` | Each strand's alt count against a beta-binomial with that strand's control rate (lightly shrunk toward the pooled rate for that substitution class, never below it) and overdispersion ρ estimated between control libraries. Strand p-values are combined by Fisher's method (`strand_combine="both"` takes the larger instead). |
| `pbinom` | Exact Poisson-binomial over every read's error probability from base quality (LoFreq's test). Needs no controls. |
| `+filter` | Binomial test of the call's plus/minus split against the column's strand mix (needs ≥3 reads, α = 0.01). |

### Calibrating `qcal`

Quality scores overstate certainty: taken at face value, `qcal` made 29.5 false calls per J꜀4 sample in the spike-in test. So the control rates are multiplied by a factor c chosen by leave-one-control-out: each control is called against a null fitted on the others, and c is the smallest value on the grid 1, 2, 3, 5, 10, 20, 50 at which no more than 1% of tested positions reach p < 0.01. J꜀4 gave c = 10. The factor is saved with the panel and reported in `summary.json` (`qscale`, with the full trace).

## Spike-in validation

`shmtrace-validate` plants known substitutions into reads from a control that is *not* in the null, choosing sites as AID would (WRC/GYW weighted 4×, WA/TW 2×, G:C transitions favoured) and writing each one into 1, 2, 3, 5, 10 or 20 reads with the read's own base quality kept. The caller then runs as normal. 60 sites × 3 draws.

**J꜀4** (spiked into Aicda-/- KO1, null from KO2 and KO3):

| Mode | 1 read | 2 | 3 | 5 | 10 | 20 | Recall | False calls per sample | FDR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| geneious | 0.80 | 0.87 | 0.90 | 1.00 | 1.00 | 1.00 | 0.93 | 2.7 | 0.05 |
| v0.1 | 0.00 | 0.07 | 0.50 | 0.83 | 1.00 | 1.00 | 0.57 | 3.0 | 0.08 |
| **qcal** | 0.03 | 0.43 | 0.77 | 0.97 | 1.00 | 1.00 | 0.70 | **0.0** | **0.00** |
| strand | 0.00 | 0.00 | 0.17 | 0.77 | 1.00 | 1.00 | 0.49 | 2.0 | 0.06 |

**N114** (spiked into empty-vector C4, null from C3_1, C3_2 and D1):

| Mode | 1 read | 2 | 3 | 5 | 10 | 20 | Recall | False calls per sample | FDR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| geneious | 0.50 | 0.83 | 1.00 | 0.97 | 1.00 | 1.00 | 0.88 | 6.0 | 0.10 |
| v0.1 | 0.07 | 0.37 | 0.67 | 0.97 | 1.00 | 1.00 | 0.68 | 1.0 | 0.02 |
| **qcal** | 0.20 | 0.90 | 1.00 | 1.00 | 1.00 | 1.00 | 0.85 | **0.3** | **0.01** |
| strand | 0.00 | 0.13 | 0.40 | 0.90 | 1.00 | 1.00 | 0.57 | 0.0 | 0.00 |

The `+filter` variants gave identical numbers: planted mutations are spread over both strands, and the background calls left after testing were not strand-biased.

What this shows:
- `qcal` dominates the v0.1 default on both amplicons: more recall at every support level and fewer false calls.
- Geneious mode finds more single-read mutations, but those come with 3–6 false calls per 500-read sample. In a real Aicda-/- library that is most of what Geneious reports.
- `strand` is the most conservative. Splitting the reads in two halves the power at each strand, so it only pays off where strand-specific error is a real problem, which it is not in these runs.
- No mode can confirm a single-read mutation at 500 reads against this background without false calls. Single-read calls are the "unique mutations" that SHMTool and the manuscript count, so they need a control-matched background rate subtracted, not a per-site test.

Caveat: planted bases are idealised. A real mutation's read can carry its own alignment or quality artefacts, so real recall will be somewhat lower.

## Benchmark against Geneious (real libraries)

`shmtrace-benchmark` on the deposit at matched depth (J꜀4 498 reads, N114 492 reads, on the reads in Geneious's own BAMs). Each control is called against a null fitted from the *other* controls of its run, so control calls are an honest background estimate. Jaccard is site-level overlap with the Geneious export.

**J꜀4**

| Mode | Wt2 / Wt3 / Wt4 calls (Geneious 161 / 342 / 68) | Jaccard WT | KO calls (Geneious 2 / 6 / 3) | WT vs KO density ×10⁻⁵ | Fold |
|---|---|---|---|---|---:|
| geneious | 172 / 344 / 69 | 0.97–0.98 | 3 / 6 / 3 | 74.9 vs 1.52 | 49 |
| v0.1 | 72 / 297 / 68 | 0.38–0.97 | 0 / 3 / 1 | 56.2 vs 0.51 | 110 |
| **qcal** | **156 / 347 / 69** | **0.85–0.97** | **0 / 2 / 1** | **73.4 vs 0.38** | **192** |
| strand | 51 / 296 / 69 | 0.29–0.97 | 0 / 3 / 1 | 53.5 vs 0.51 | 105 |

**N114**

| Mode | AID density ×10⁻⁵ | Empty-vector density ×10⁻⁵ (calls in 4 controls) | Jaccard AID |
|---|---:|---|---|
| geneious | 28.9 | 1.57 (2, 1, 6, 10) | 0.85–0.97 |
| v0.1 | 11.6 | 0 | 0.24–0.51 |
| **qcal** | **24.2** | **0** | 0.41–0.75 |
| strand | 8.1 | 0 | 0.13–0.44 |

`+filter` changed at most one call per sample in either amplicon.

**Decision: `recommended` now means `qcal`.** It is the only model that keeps Geneious-level wild-type signal (it recovers the Wt2 calls the v0.1 default lost) while cutting the control background to a quarter of Geneious's. The v0.1 behaviour remains available as `--mode v0.1`. On N114 the D3 libraries keep 43–76% of Geneious's calls; the dropped ones are mostly single-read calls, which in the spike-in are where Geneious's false calls come from.

## Strand findings

- Deposit reads are about 50/50 strand on both amplicons.
- `fit_aligned` counts sites where one strand's control rate is far above the other's (`strand_specific_error_sites` in the null stats).
- Multi-read Geneious calls are not strand-biased. Single-read calls lean to the plus strand (54 of 76, p = 0.0008; see `CONFOUNDING.md`).

## Where the reads come from (Plasmidsaurus)

POD5/FAST5 files are not available, so the basecaller had to be inferred:
- Plasmidsaurus sequences amplicons on R10.4.1 flow cells with v14 chemistry and basecalls with "the latest" super-accurate model, without naming a version ([amplicon FAQ](https://plasmidsaurus.com/faq/amplicon), [results guide](https://app.plasmidsaurus.com/results_interpretation)).
- The FASTQs cap quality at Q50 with 41–59% of bases at Q40 or above, which fits Dorado sup models. Headers are renamed (`<n>_<RUN>_<sample>`), so the run ID and basecaller tags are gone.
- Reads below read-level Q10 are dropped before delivery. All of these runs used the Premium PCR service (ligation prep, all reads delivered, about 5,000 per sample), as confirmed by Abdulmateen.

Consequence: the error profile can shift between orders when Plasmidsaurus updates its basecaller. Background panels are tied to a run, same-run controls remain the default, and `qcal` re-measures how trustworthy the quality scores are for each control set rather than assuming it.

## Reproduce

```bash
shmtrace-validate --control KO1.fastq --null KO2.fastq --null KO3.fastq \
  --reference jh4-base.fasta --contaminant Contaminant_JH4_ref.fasta --preset jh4_intron \
  --profiles geneious,recommended,qcal,qcal+filter,strand,strand+filter --reps 3 --sites 60 --out val_jh4
shmtrace-benchmark --sheet jh4_sheet.csv --reference jh4-base.fasta --contaminant Contaminant_JH4_ref.fasta \
  --preset jh4_intron --modes geneious,recommended,qcal,qcal+filter,strand,strand+filter --out bench_jh4
```
