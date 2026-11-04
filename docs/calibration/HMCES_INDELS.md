# Indels against a known biological effect

*12 October 2026; rerun 2 November 2026 (see "Rerun with the November 2026 phasing"). Data: [PRJNA820976](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA820976), "HMCES protects immunoglobulin genes specifically from deletions during somatic hypermutation". Illumina MiSeq 2×300, human Ramos B cells, IGH VDJ. All 30 FASTQ files checked against ENA's MD5.*

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

## Rerun with the November 2026 phasing

*2 November 2026.* Every Ramos clone, the AID-knockout ones included, carries a second IGH allele: a block of 25–30 linked differences in close to half its reads. Phasing drops a block's reads only while that keeps it under `phase_max_drop` (45% of the library). Until this change, a library where the allele sat above that refused to phase and the allele's sites were then called as substitutions: eight of the fifteen clones, AIDKO_3 among them. Phasing now leaves the sites of such a block out of calling, as a second allele's should be. The benchmark run was regenerated with that code and nothing else changed:

| | Before | After |
|---|---:|---:|
| Deletions per clone, Wt / HMCES-KO | 23.2 / 87.2 | 25.8 / 87.2 |
| HMCES-KO vs Wt, deletions (Welch on log) | 3.8-fold, P = 1.0 × 10⁻⁴ | 3.4-fold, P = 2.9 × 10⁻⁴ |
| Insertions per clone, Wt / HMCES-KO | 5.5 / 13.7 | 5.5 / 14.2 |
| Substitutions per 10³ bp, Wt / HMCES-KO / AID-KO (means) | 2.47 / 6.45 / 10.40 | 0.67 / 0.34 / 0.01 |
| AIDKO_3 as a control | 36 calls, 31.2 per 10³, 61% G:C | 8 calls, 0.04 per 10³, 75% G:C |
| AIDKO_3 analysed as a sample (`hmces_k3`) | 63 calls, 31.4 per 10³, 57% G:C | 35 calls, 0.26 per 10³, 54% G:C |

The "Before" column is the run the paper reported in October, which is not the table above (that one is from 12 October, an earlier version still). The deletion effect stands; the four Wt clones that phase gain a few deletions each because phasing now finds the allele as two to four blocks and drops a few dozen more reads, so a slightly different 20,000 are drawn.

**This corrects the section above.** Most of AIDKO_3's excess was the second allele, not contamination or a mislabel: its substitution frequency falls from 31.2 to 0.04 per 10³ once the allele is out. Some excess remains. As a sample it still makes 35 calls where its two siblings make 0 and 2, only 54% at G:C against 90–97% in the mutating clones, and the noisy-control rule still sets it aside as a decoy. So it is a control with a background its siblings lack, but a small one, and the group means in the first table that include it are no longer what a default run gives.

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
