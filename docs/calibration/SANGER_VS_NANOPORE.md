# Do Sanger and nanopore agree in the same wild-type mice?

*12 October 2026. Sanger: the lab's J꜀4 clones from their `.ab1` traces, organised by mouse with `organize_jh4_ab1.py`, plus the GenBank records of clones whose traces are not in the trace folder. Nanopore: the six bundled libraries at 497 reads, in the recommended mode and, for comparison with the manuscript, the Geneious rule. Window 101–629 (529 bp) on both platforms. Figures: `figures/sanger_vs_nanopore/<mode>/`; script: `sanger_vs_nanopore.py`.*

## Answer

**Yes, on everything both platforms can measure.** They agree on:

- the substitution spectrum;
- where along the intron mutations fall;
- hotspot preference;
- which mouse is most mutated.

They also agree on something the manuscript could not see. **WT2 is dominated by one unmutated lineage, on both platforms**: 30 of its 34 traced Sanger clones carry no mutation, and one nearly unmutated clone holds 85% of its nanopore reads.

| Test | Recommended mode | Geneious rule | Manuscript (Ver13) |
|---|---|---|---|
| Unique mutations, Sanger / ONT | 415 / 548 | 415 / 585 | 515 / 571 |
| 12-type spectrum, Pearson r [95% CI] | **0.90 [0.66–0.97]** | **0.93 [0.76–0.98]** | 0.97 [0.91–0.99] |
| Slope (ONT on Sanger) | 0.98 | 1.02 | 1.11 |
| Mean / largest difference per type | 1.7 / 3.8 points | 1.5 / 3.0 points | ≤ 1.4 points |
| Platform × type interaction (RM ANOVA) | p = 0.18 | p = 0.20 | p = 0.875 |
| Smallest paired, Šídák-adjusted type difference | p = 0.49 (G>A) | p = 0.67 | p = 0.106 |
| G:C transitions / transversions, adjusted p | 0.19 / 1.00 | 0.17 / 0.95 | 0.80 / 0.98 |
| Position along the window, 25-bp bins, r | **0.89** (p = 6 × 10⁻⁸) | 0.89 | not tested |
| G:C mutations in WRC/GYW (window share 20%) | Sanger 58%, ONT 56% | Sanger 58%, ONT 55% | not tested |
| A:T mutations in WA/TW (window share 50%) | Sanger 80%, ONT 70% | Sanger 80%, ONT 69% | not tested |
| Mouse ranking by load | WT3 > WT4 > WT2 on both | same | not tested |

WT2's Sanger spectrum rests on 6 mutations, so its per-mouse proportions are noise. That is why the paired, per-mouse tests (panel B's error bars, the interaction) are looser than the pooled correlation.

## The Sanger data

The traces are in five sequencing orders (`berrumar_432625`, `434138`, `434141`, `434355`, `434545`; July 2024), 288 files. Each trace stores the order number, not the clone name, so a clone's identity is its file name and its plate well.

| What | Count |
|---|---:|
| Traces | 288 |
| Usable J꜀4 traces | 209 |
| J꜀4 paralogue or other contaminant | 64 |
| Failed reads (no or too short an alignment) | 15 |
| Wild-type and *Aicda*−/− GenBank records with no trace in these orders | 41 (mostly the `jhextra-w3` and `jhextra-w4` batches) |

Clone labels were reused between orders: "w4-2" in three orders is three different colonies (3, 31 and 0 mutations, none shared). Every trace is therefore one colony, identified by order and well.

The GenBank records are faithful to the traces. Inside the window, 87 of 95 differences are a low-quality (Q < 20) or N call that was corrected by hand. One record, `3-w2-3`, carries four C's in 9 bp (591–599) where the trace reads the reference at Q23–27. The trace wins, and that trace fails alignment anyway.

## Which clone belongs to which mouse

A trace goes to the mouse its label names, unless its sequence contradicts the label. Four do:

| Trace | Label | Assigned | Evidence |
|---|---|---|---|
| 434141 E8 | w2-1 | **WT3** | 18 mutations; WT2's other 34 traces carry at most 3, while 10 of WT3's 19 carry 13 or more. It sits where w3-1 belongs (G8 w3-3, H8 w3-4, A9 w3-5), and w2-1 is already at A7 on the same plate. Filed under W2 by hand, so the manuscript counted it as WT2. |
| 434141 F8 | w2-2 | **WT3** | 13 mutations; same plate logic (w2-2 is already at B7). |
| 434141 A11 | w4-9 | **MKLN K2** | Carries the same 7-mutation lineage as k2-1 (E11, same plate, where the K2 series starts right after the WT4 column) and k2-24 (order 434355). Not in any WT4 nanopore read. |
| 432625 E6 | a1-1 | **unassigned** | 37 AID-type mutations (64% of G:C changes at WRC/GYW), impossible for an *Aicda*−/− mouse. Not the paralogue: 50 differences from J꜀4, 93 from the paralogue. It shares no lineage with any clone, so its mouse is unknown. It had no GenBank record, so the manuscript did not use it. |

The GenBank records `88-w2-28` and `96-w2-36` (filed as WT2, 22 and 24 mutations) are **renamed copies of K4 traces** 434355 H11 `88-k4-28` and H12 `96-k4-36`: identical sequences, same well-order numbers, same clone numbers. They are excluded.

Every other wild-type clone sits in the mouse its name gives. This holds even where the hand-made folders disagree: "wt3" clones filed under W2 match WT3's load (p = 0.25) and not WT2's (p = 0.0008), and a 25-mutation WT3 clone appears in both the W2 and the W3 folders. Shared lineages (five or more identical mutations) never cross mice apart from the cases above.

The organised set (`~/Nanoshm/sanger_ab1_by_mouse`, with `manifest.csv` and `README.md`):

| Mouse | Traces | GenBank-only clones | Usable J꜀4 clones | Unique mutations |
|---|---:|---:|---:|---:|
| WT2 | 34 | 1 | 35 | 6 |
| WT3 | 21 (2 relabelled) | 20 | 39 | 303 |
| WT4 | 13 | 17 | 25 | 106 |
| *Aicda*−/− A1 | 5 | 1 | 6 | 1 |
| *Aicda*−/− A2 | 6 | 1 | 7 | 1 |
| *Aicda*−/− A3 | 0 | 1 | 1 | 0 |

## AID dependence on each platform

| | Wt geometric mean [95% bootstrap CI] | *Aicda*−/− | Fold | P (Welch, log) |
|---|---|---|---:|---:|
| ONT, recommended | 56.7 [26–127] × 10⁻⁵ | 0.55 | 103 | 0.0017 |
| ONT, Geneious rule | 60.7 [26–131] | 1.43 | 42 | 0.0060 |
| **Sanger** | 337 [32–1469] | ≤ 55 | ≥ 6 | **0.27** |
| *Manuscript: ONT / Sanger* | *59 / 828* | *1.2 / 20* | *47 / 41* | *0.0037 / 0.0045* |

**With the clones assigned correctly, Sanger alone no longer separates wild type from *Aicda*−/−.** WT2's own clones carry 6 mutations in 35 clones (32 × 10⁻⁵), the knockouts' level. The manuscript's Sanger WT2 was inflated by WT3 and MKLN clones filed under W2, and its knockout denominators by JH6 clones.

Nanopore separates the groups clearly, because 497 reads also catch WT2's minority mutated lineages (146 unique mutations). This is the manuscript's own argument for depth, and the corrected data make it more strongly than the published figure did.

## WT2: the same expanded clone on both platforms

| | Sanger | Nanopore |
|---|---|---|
| Sequences with no mutation | 30 of 34 traced clones (88%) | the largest clone holds 85% of junction reads, 7% of them mutated ([`CLONES.md`](CLONES.md)) |
| Mutations per 10³ bp, per sequence | 0.3 | 1.7 |
| Per clone (nanopore junction) | — | 12.0 |

Two independent PCRs and two platforms see one almost unmutated lineage dominating this mouse's germinal-centre B cells. Sanger's 34 clones, almost all from that lineage, miss the rest; nanopore's 497 reads do not.

## What changed from the manuscript's Sanger numbers

`figures/sanger_vs_nanopore/figS_sanger_inclusion` (from `sanger_inclusion_figure.py`) draws this: what each mouse's Sanger value was built from in the manuscript and now, each mouse's density before and after, and the group comparison.

The manuscript's values come from Geneious assemblies built folder by folder:

| Manuscript assembly | Clones | By clone name |
|---|---:|---|
| W2 | 35 | 20 WT2, **14 WT3**, 1 w1 |
| W3 | 33 | **11 WT2**, 20 WT3, 2 w1 |
| W4 | 24 | 24 WT4 |
| AID 1 / 2 / 3 | 10 / 9 / 9 | 4 A1 + **6 JH6**; 6 A2 + **3 JH6**; 1 A1, 1 A2, 1 A3 + **6 JH6** |

Its denominators (37/36/38) also count 19 paralogue sequences as unmutated J꜀4 clones. Its W2 assembly additionally held the two WT3 clones at 434141 E8/F8 and the two renamed K4 clones. Figures 2B and 3 and Supplementary S2 should be regenerated from the organised set. The spectrum and positional conclusions stand. The Sanger knockout comparison does not: it should be reported as underpowered.

## A bug this comparison found in nanoshm

A Sanger sample of fewer than four clones cannot apply the germline rule to itself, so nanoshm switched the rule off. A3's single clone then reported the strain's five germline differences from the reference as five mutations. Now a sample with too few clones uses the germline changes measured across the project's other clones (`project.sanger_germline`; test in `tests/test_platforms.py`).

## Reproduce

```bash
python docs/calibration/organize_jh4_ab1.py "~/Documents/WT JH4 SANGER AB1" ~/Nanoshm/sanger_ab1_by_mouse ~/Nanoshm/sanger_raw
python docs/calibration/sanger_vs_nanopore.py --ont-mode recommended
python docs/calibration/sanger_vs_nanopore.py --ont-mode geneious
```

The same folders run in the app. Drop a mouse's folder (or choose it) as one sample, set the platform to Sanger, and the `.ab1` traces and the GenBank-only FASTA in it are read together.

## Limits

- **Some traces are missing.** 41 wild-type and knockout clones exist only as GenBank records, read at a flat Q30. Their traces (the `jhextra-w3`/`jhextra-w4` orders and a few others) would let their quality be used too.
- **The Sanger reads start at J꜀4**, so they carry no V(D)J junction, the strongest identity check. A misfiled clone with few mutations and no relatives in its true mouse cannot be caught; mutation load still catches a misfile between mice as different as WT2 and WT3.
- **Three mice per group.** The agreement tests are descriptive, as in the manuscript, and the repeated-measures ANOVA assumes sphericity.
