# Is the filtering selecting real SHM?

This check was run on 11 October 2026 on the six bundled J꜀4 libraries at 500 reads each, with phasing off and leave-one-out knockout nulls. To reproduce it, run `python docs/calibration/confounding_check.py` from the repository root.

The logic: AID deaminates C most often in WRC/GYW motifs, and Pol η introduces A:T mutations most often in WA/TW motifs. Sequencing error has no reason to prefer those motifs. So if the filter keeps real SHM, the kept calls should be enriched at these hotspots and the rejected candidates should look like knockout error.

## Hotspot targeting

| Set | n | G:C changes at WRC/GYW | A:T changes at WA/TW | Share at G:C |
|---|---|---|---|---|
| Expected if positions were random | | 0.20 | 0.50 | |
| Knockout, all candidates (error) | 727 | 0.16 | 0.54 | 0.57 |
| WT, kept by both modes | 426 | **0.57** | **0.71** | 0.35 |
| WT, kept by Geneious mode only | 159 | 0.45 | 0.66 | 0.28 |
| WT, kept as a single-read call (Geneious mode) | 76 | 0.52 | 0.57 | 0.30 |
| WT, rejected by Geneious mode | 738 | 0.12 | 0.58 | 0.53 |

- Kept wild-type calls are at AID hotspots 2.8 times as often as random positions. They are dominated by A:T changes, as expected for J꜀4 SHM.
- Rejected wild-type candidates look like knockout error. Their hotspot rate is at or below random and they lean towards G:C changes.
- **The recommended mode discards real SHM.** If the Geneious-only calls are a mix of SHM (hotspot rate 0.57) and error (0.16), their observed rate of 0.45 puts about 70% of them as real. Most of them are seen in only 1 or 2 reads.

## Recurrence at knockout error-prone positions

10% of all possible position–substitution slots are candidates in at least 2 of the 3 knockouts.

| | Wt2 | Wt3 | Wt4 |
|---|---|---|---|
| Kept calls at those slots | 7–8% | 9–10% | 13% |
| Rejected candidates at those slots | 29% | 28% | 26% |

Kept calls fall at error-prone positions no more often than chance. Rejected candidates are enriched there almost 3-fold.

## Quality

| Sample | Median read mean Q | Single-read candidates with variant base ≥ Q40 |
|---|---|---|
| Wt2 | 40.8 | 17% |
| Wt3 | 40.7 | 17% |
| Wt4 | 40.7 | 1% |
| KO1–3 | 40.8–40.9 | 1% |

Overall read quality is the same in wild type and knockout, so group differences are not explained by run or library quality. The quality filter works because real substitutions are called at normal base quality, while errors mostly sit at Q8–Q11.

## Confounds this cannot rule out

1. **No UMIs.** A call supported by many reads may come from an expanded clone, or from PCR copies of one molecule. Unique-call counts are robust to this; recurrence-weighted spectra are not.
2. **Strand bias (corrected 2026-10-11).** An earlier version said strand bias could not be tested. That was true only of the repository's example FASTQs (`examples_compare/`, all forward). The deposited JH4 sub500 FASTQs are mixed, 51–56% plus strand, so it can be tested (`strand_check.py`). In `--mode geneious`, none of the 517 multi-read calls across the six JH4 libraries is strand-biased (binomial p < 0.01 against the column's strand mix). Single-read calls are not: 54 of 76 are on the plus strand where about 52% is expected (p = 0.0008). So part of the single-read Q40+ calls is a strand-specific error, and a both-strand requirement or strand-aware null would remove it.
3. **The thresholds were chosen on these same six mice.** The background estimate is optimistic until it is checked on a new run.
4. **One run, one locus.** Wild type and knockout share a flow cell, which is good for the comparison. The error rates may still not transfer to another run or locus.
5. **These files were already screened for the contaminant haplotype.** On raw reads, a haplotype block is the biggest artefact. In the July calibration it accounted for about a 7-fold drop in background.
