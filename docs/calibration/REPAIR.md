# Repair routes from the spectrum

*12 October 2026. `shmtrace/repair.py`, `shmtrace-repair`.*

## What it reads

After AID turns a C into a U, the route that resolves it leaves its own mark:

| Route | Signature | Knockout that removes it |
|---|---|---|
| Replication over U | G:C transitions (C>T, G>A) | none; it is all that remains in Ung⁻/⁻Msh2⁻/⁻ |
| UNG excision, then translesion synthesis | G:C transversions; REV1 inserts dC, so mostly C>G and G>C | Ung⁻/⁻: G:C mutations become almost all transitions |
| MSH2/MSH6–EXO1, then Pol η | A:T mutations, more at WA than TW (Pol η errs on one strand) | Msh2⁻/⁻ or Polh⁻/⁻: A:T mutations fall sharply |

**Indices, per animal and per group:**

| Index | Definition |
|---|---|
| `pct_AT` | A:T mutations as a share of all |
| `pct_GC_Ts_of_GC`, `pct_GC_Tv_of_GC` | transitions and transversions as a share of G:C mutations |
| `pct_CG_of_GC_Tv` | C>G plus G>C as a share of G:C transversions |
| `A_over_T` | A mutations over T mutations |
| `WA_over_TW` | mutations at WA over those at TW |
| `WRC_over_GYW` | mutations at WRC over those at GYW |

Groups are compared animal by animal (Welch and permutation). A Fisher test on pooled mutations is shown only as the pseudoreplicated contrast.

**Decomposition:** a group's 12-type spectrum is fitted as a non-negative mixture of reference spectra, with weights summing to 1. The references are other groups in the same project, for example Ung⁻/⁻, Msh2⁻/⁻ or Ung⁻/⁻Msh2⁻/⁻ mice sequenced alongside. The 95% intervals come from resampling animals, and mutations within animals, in both the group and its references. On simulated mixtures the weights are recovered to within 0.05 (`tests/test_repair.py`).

```bash
shmtrace-repair project.json --case MKLN --control WT
shmtrace-repair project.json --case MKLN --reference ung_msh2=UngMsh2 --reference msh2=Msh2 --reference wt=WT
```


**In the app:** the *Repair* tab, after an analysis. Pick a group, the group to compare it with, and optionally reference groups for a decomposition.

## First use: MKLN against wild type (Sanger, J꜀4)

You expected more G:C transversions in the MKLN knockouts. Per mouse:

| Index | MKLN (4 mice) | Wt (3 mice) | Welch p | Permutation p |
|---|---:|---:|---:|---:|
| G:C transversions, % of G:C | 57.9 | 41.1 | 0.25 | 0.31 |
| A:T, % of all | 53.5 | 62.5 | 0.31 | 0.34 |
| C>G + G>C, % of G:C transversions | 50.2 | 58.7 | 0.56 | 0.63 |
| WA / TW | 2.26 | 1.18 | 0.45 | 0.89 |

**The direction matches the expectation, but it isn't significant with 4 vs 3 mice.** The mice disagree among themselves: K3 (88%) and K4 (63%) are high, while K1 (38%) and K2 (43%) are in the wild-type range (29–49%). K2 and K3 rest on only 23 and 19 mutations. Testing pooled mutations instead gives p = 0.12; the per-mouse spread is the honest uncertainty. Making this a finding would need more mice, or reference spectra from repair-deficient mice to decompose against.
