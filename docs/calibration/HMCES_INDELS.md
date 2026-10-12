# Indels against a known biological effect

*12 October 2026. Data: [PRJNA820976](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA820976), "HMCES protects immunoglobulin genes specifically from deletions during somatic hypermutation". Illumina MiSeq 2×300, human Ramos B cells, IGH VDJ. All 30 FASTQ files checked against ENA's MD5.*

The spike-ins in [`INDELS.md`](INDELS.md) show the indel caller recovers indels that were planted into real reads. They cannot show it tracks a real biological difference. This set can: knocking out HMCES is published to raise SHM-associated deletions specifically, and the same study provides AID-knockout clones as the null.

It is also the first test of the indel caller on a platform it was not calibrated on. Everything in `INDELS.md` was measured on nanopore.

## The design

Fifteen Ramos clones: 6 wild type, 6 HMCES-knockout, 3 AID-knockout. Pairs merge to a 471-bp VDJ amplicon (97% merge). The reference is the majority consensus of the AID-knockout reads, which carry no SHM by construction. Settings: `--platform illumina`, `generic_amplicon`, window 26–446, 20,000 reads per clone from the first 100,000 pairs, indel calling on, the three AID-knockout clones as controls.

## Result

| Group | n | Indel density ×10⁻⁵ | Deletions per clone | Deletions as a share of indels | Substitutions per 10³ bp | Mutations at G:C |
|---|---:|---:|---:|---:|---:|---:|
| Wild type | 6 | 0.358 | 23.8 | 79% | 2.00 | 94% |
| **HMCES-KO** | 6 | **1.240** | **89.3** | **86%** | 5.63 | 91% |
| AID-KO (clean) | 2 | **0.000** | **0.0** | – | **0.00** | – |
| AID-KO (all 3) | 3 | 0.084 | 3.7 | 52% | 10.39 | 61% |

| Comparison | Welch on log values | Mann-Whitney |
|---|---:|---:|
| HMCES-KO vs wild type, indel density | **P = 0.000025** | P = 0.0022 |
| HMCES-KO vs wild type, deletions | **P = 0.000052** | P = 0.0022 |
| HMCES-KO vs wild type, insertions | P = 0.0054 | P = 0.0049 |
| HMCES-KO vs clean AID-KO, deletions | P = 0.000007 | – |
| Wild type vs clean AID-KO, deletions | P = 0.000111 | – |

**The published effect is reproduced, including its specificity.** HMCES-KO clones carry 3.5 times the indel density of wild type. Deletions rise 3.75-fold (23.8 to 89.3 per clone) against 2.3-fold for insertions (6.2 to 14.3), and deletions grow from 79% to 86% of all indels. Deletions are the stronger effect by two orders of magnitude in p-value, which is what "specifically from deletions" predicts.

**The AID-null controls behave.** Two of the three AID-knockout clones yield zero indels and zero substitution calls in 20,000 reads each. Both wild type and HMCES-KO separate from them at P ≤ 0.0002.

**The substitutions look like AID's.** 91–94% of calls in the mutating clones are at G:C, the expected Ramos SHM signature, and the AID-knockout clones contribute none.

Deletion sizes are large in both mutating groups (median 11–19 bp, maximum 123–129 bp) and span the amplicon, consistent with AID-initiated lesions resolved by end joining rather than with sequencing error. In the AID-knockout clones the few indels are short (median 3 bp, maximum 13 bp).

## One control is an outlier, and the spectrum shows it

AIDKO_3 is out of line with its two siblings: 33 substitution calls and 21 indels against zero and zero. Its substitutions are only 61% at G:C, against 91–94% in the genuinely mutating clones, so they do not look like AID's either. The library appears contaminated or mislabelled.

It is reported here both ways. Including it, the AID-knockout group mean substitution frequency (10.39 per 10³ bp) exceeds wild type (2.00), which is obviously wrong; excluding it, the controls are clean and every comparison strengthens. The HMCES-versus-wild-type conclusion does not depend on it, since that comparison never uses the controls except to fit the null.

Two honest notes about the automatic checks:

- An earlier version of this analysis, run before the calibration bug in [`../../CHANGELOG.md`](../../CHANGELOG.md) was fixed, flagged AIDKO_3 through `clustered_calls` (13 calls in 10 bp, p = 4.5 × 10⁻⁸). With the corrected, more conservative null it makes fewer calls and no longer trips that threshold. **The positional-clustering flag is therefore not a complete safety net.**
- What does catch it in both versions is the spectrum: `pct_GC` of 61% where the real signal is 91–94%. Comparing a control against its siblings remains the most reliable check, and no automatic test replaces looking at the control group.

## Limits

- **One locus, one cell line, one run.** The effect size here should not be read as a general indel rate.
- **No deduplication.** As in [`ILLUMINA.md`](ILLUMINA.md), these amplicons are resequenced far beyond one read per molecule, so distinct-event counts are inflated relative to independent molecules. The comparison between groups is still valid because every group is treated identically at the same depth.
- **Indel counts are per analysed clone at 20,000 reads**, not per molecule, and are not comparable with the nanopore numbers in `INDELS.md`, which are at 497 reads.

## Reproduce

```bash
python /Users/abdulmateenaderinto/Nanoshm/public_data/run_prjna820976.py 20000
```

The download, MD5 check, subsampling and reference-building steps sit with the data, outside the repository.
