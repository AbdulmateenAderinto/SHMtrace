# Group tests with the animal as the unit, and power

*12 October 2026. `shmtrace/replicates.py`, `shmtrace/power.py`; calibration `docs/calibration/replicates_calibration.py`, output in `replicates_calibration/`.*


**In the app:** the *Power* tab. Take the pilot from a group of the current analysis (choose the measure), or type the numbers in.

## Why

Reads and clones from one mouse are not independent. A clone's reads share its mutations, and mice differ from each other far more than reads within a mouse do. Testing every read of one group against every read of another (pseudoreplication) treats thousands of correlated observations as independent.

## The tests (`animal_tests` in every project result)

For each count-based measure, every non-control group is compared with the control group, one value per animal. The two measures are mutations per base, and unique mutations per base sequenced.

| Test | What it assumes |
|---|---|
| **Negative binomial**, p by parametric bootstrap | Counts are negative binomial around each group's rate, with one dispersion shared by the groups (a constant coefficient of variation between animals). The dispersion is fitted under the full model, as DESeq2 and edgeR do. The bootstrap replaces the chi-square approximation, which is unsafe with 2–10 animals per group. |
| **Permutation** on log10 rates | Nothing. Exact for small groups. With 3 vs 3 animals the smallest two-sided p is 0.1: six animals can't show more than that without a model. |
| **Welch** on log10 rates | Log-normal animals (the manuscript's test). |
| **Fold** | The ratio of the groups' fitted rates, and the ratio of geometric means with a 95% interval on the log scale. |
| *Pooled reads* | Reported only as a contrast. |

## Are the p-values honest?

Simulated experiments: animals log-normal around their group's rate with a between-animal CV of 1 (as in the J꜀4 wild types), negative binomial counts on top, and 500 reads × 526 bp per animal. Share of 300 experiments per row with p < 0.05:

| Animals per group | True fold | Welch on log10 | Negative binomial | Permutation | Pooled reads |
|---:|---:|---:|---:|---:|---:|
| 3 | 1 (no difference) | 0.037 | **0.050** | 0.000 | **0.913** |
| 3 | 2 | 0.093 | 0.120 | 0.000 | 0.960 |
| 3 | 4 | 0.207 | 0.297 | 0.000 | 0.990 |
| 5 | 1 (no difference) | 0.043 | **0.057** | 0.047 | **0.937** |
| 5 | 2 | 0.150 | 0.220 | 0.157 | 0.963 |
| 5 | 4 | 0.463 | 0.557 | 0.473 | 1.000 |

- **Negative binomial:** false positives at the nominal 5%, and the most power of the honest tests.
- **Welch on log10:** slightly conservative.
- **Permutation:** exact, but powerless with 3 vs 3.
- **Pooled reads:** "finds" a difference in 91–94% of experiments where there is none. That is what pseudoreplication does.

On the real J꜀4 libraries (Wt vs Aicda−/−, 3 vs 3):

| Measure | Fold [95% CI] | Negative binomial | Welch | Permutation | Pooled reads |
|---|---|---:|---:|---:|---:|
| Mutation frequency | 188 [10–4,250] | p = 0.022 | 0.009 | 0.1 (floor) | p ≈ 0 |
| Unique density | 138 [16–1,440] | p = 0.006 | 0.004 | 0.1 (floor) | p ≈ 0 |

## Power: how many mice?

`shmtrace-power` simulates whole experiments analysed with the animal-level test. It takes the rate, the variation between animals and the clumping within animals from a pilot project.

From the J꜀4 wild types (3 mice, 500 reads) as the pilot, power at α = 0.05 (Welch on log10, 1,000 simulated experiments per cell):

| Measure (between-animal CV) | Fold | 3 mice | 5 | 8 | 10 |
|---|---:|---:|---:|---:|---:|
| Mutation frequency, read-weighted (CV 2.0) | 2 | 0.06 | 0.12 | 0.16 | 0.21 |
| Unique density (CV 0.95) | 2 | 0.10 | 0.21 | 0.34 | 0.45 |
| Clone-weighted frequency (CV 1.0) | 2 | 0.09 | 0.19 | 0.33 | 0.42 |

These are at 500 reads per mouse; 250 or 1,000 reads change power by at most 0.05.

What this means for designing experiments:
- **Mice matter, reads don't.** Past about 250 reads per mouse, more reads barely help. The limit is how much mice differ from each other.
- **The read-weighted mutation frequency is the noisiest measure** when clones expand unevenly between mice (Wt2's dominant unmutated clone gives a 12-fold spread). Unique density or the clone-weighted frequency halve the between-mouse variation.
- **3 mice per group only detects large effects.** That's enough for AID-dependence (100-fold or more), but a 2-fold change between two AID-proficient groups, such as a knockout of an AID cofactor, needs more than 10 mice per group, or a less variable measure.
- **The CV estimate is itself uncertain.** It comes from 3 mice. A pilot with more animals would tighten every number here.

```bash
shmtrace-power --pilot project.json --group WT --metric density --fold 0.5 --animals 3,5,8,10,15
```
