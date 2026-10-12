# Do Sanger and nanopore agree in the same wild-type mice?

*12 October 2026. `docs/calibration/sanger_vs_nanopore.py`. Sanger: the lab's J꜀4 clones, assigned to mice by clone name ([`SANGER.md`](SANGER.md)), with three misfiled clones removed (below), in nanoshm's Sanger mode. Nanopore: the six bundled libraries at 497 reads, in the recommended mode and, for comparison with the manuscript, the Geneious rule. Window 101–629 (529 bp) on both platforms. Figures: `figures/sanger_vs_nanopore/<mode>/`.*

## Answer

**Yes, on everything the two platforms can both measure.** They agree on:

- the substitution spectrum;
- where along the intron mutations fall;
- hotspot preference;
- which mouse is most mutated;
- that mutation depends on AID.

They don't agree on absolute mutation density. Neither did the manuscript's numbers, and they shouldn't: a Sanger clone and a nanopore read are sampled differently.

| Test | Recommended mode | Geneious rule | Manuscript (Ver13) |
|---|---|---|---|
| Unique mutations, Sanger / ONT | 395 / 548 | 395 / 585 | 515 / 571 |
| 12-type spectrum, Pearson r [95% CI] | **0.92 [0.74–0.98]** | **0.95 [0.81–0.98]** | 0.97 [0.91–0.99] |
| Slope (ONT on Sanger) | 1.04 | 1.07 | 1.11 |
| Mean / largest difference per type | 1.5 / 3.0 points | 1.4 / 2.8 points | ≤ 1.4 points |
| Platform × type interaction (RM ANOVA) | p = 0.50 | p = 0.60 | p = 0.875 |
| Smallest paired, Šídák-adjusted type difference | p = 0.15 (G>T) | p = 0.24 (G>T) | p = 0.106 |
| G:C transitions / transversions, adjusted p | 0.06 / 0.95 | 0.13 / 1.00 | 0.80 / 0.98 |
| Position along the window, 25-bp bins, r | **0.89** (p < 10⁻⁷) | 0.88 | not tested |
| G:C mutations in WRC/GYW (window share 20%) | Sanger 57%, ONT 56% | Sanger 57%, ONT 55% | not tested |
| A:T mutations in WA/TW (window share 50%) | Sanger 72%, ONT 70% | Sanger 72%, ONT 69% | not tested |
| Mouse ranking by load | WT3 > WT4 > WT2 on both | same | not tested |

## Figures

- **`fig1_density_wt_vs_ko`** (the manuscript's Figure 2A/B): wild type against *Aicda*−/− on each platform.
- **`fig2_spectrum_concordance`** (Figure 3A–E): pooled From→To matrices, 12-type proportions with paired tests, the Pearson scatter, G:C transitions and transversions, and mutations per sequence.
- **`fig3_agreement_beyond_spectrum`** (new): mouse-against-mouse spectrum correlation, position along the window, hotspot preference against the window's own motif content, and mutation load mouse by mouse, with the clone-weighted nanopore value.
- **`figS_per_mouse_matrices`** (Supplementary S2/S3): one matrix per mouse and platform, shaded within wild type and within knockout.

Each figure is written as PNG, PDF and SVG.

## AID dependence on each platform

| | Wt geometric mean [95% bootstrap CI] | *Aicda*−/− | Fold | P (Welch, log) |
|---|---|---|---:|---:|
| ONT, recommended | 56.7 [26–127] × 10⁻⁵ | 0.55 | 103 | 0.0017 |
| ONT, Geneious rule | 60.7 [26–131] | 1.43 | 42 | 0.0060 |
| Sanger | 759 [348–1612] | ≤ 70 | ≥ 11 | 0.024 |
| *Manuscript: ONT / Sanger* | *59 / 828* | *1.2 / 20* | *47 / 41* | *0.0037 / 0.0045* |

The nanopore numbers reproduce the manuscript's.

The Sanger knockout value is an upper bound. Two of the three knockout mice have no Sanger mutation at all, so they are plotted at one mutation's worth. With 4 clones (A1) and 1 clone (A3), that floor is 47 and 189 × 10⁻⁵. Sanger cannot see a background this low with 10 clones. This is the manuscript's own point about depth, and it shows more sharply once the clones are counted correctly (next section).

## What changed from the manuscript's Sanger numbers, and why

The manuscript's Sanger values come from Geneious assemblies built folder by folder (`W2/W3/W4 Sanger assembled to Jh4 -529 bp`, `AID 1–3 assembled to …`). Their BAMs show what each assembly holds:

| Manuscript assembly | Clones in BAM | By clone name |
|---|---:|---|
| W2 | 35 | 20 WT2, **14 WT3**, 1 w1 |
| W3 | 33 | **11 WT2**, 20 WT3, 2 w1 |
| W4 | 24 | 24 WT4 |
| AID 1 | 10 | 4 A1, **6 JH6 amplicon** |
| AID 2 | 9 | 6 A2, **3 JH6 amplicon** |
| AID 3 | 9 | 1 A1, 1 A2, 1 A3, **6 JH6 amplicon** |

The manuscript's denominators (37/36/38 = 111) are these 92 clones plus the 19 paralogue sequences. Its legend calls those an Agouti-background haplotype. Either way they are not the C57BL/6 J꜀4 target, and they carry no SHM even in wild type.

This has three consequences for the manuscript:

1. **"W2" and "W3" are each a mixture of two mice.** Pooled wild-type figures are barely affected, because the same clones are in the pool. The mouse-matched tests, though (paired t, the repeated-measures interaction, per-mouse matrices S2), pair a nanopore mouse with a Sanger mixture of two mice. Assigned by name, the paired analysis above reaches the same conclusion, so the claim stands; the numbers in the manuscript's Figure 3 and Table should be replaced.
2. **The *Aicda*−/− Sanger denominators are 54% JH6 clones** (15 of 28), a different amplicon. That made the knockout density look about 2–3 times lower than its own J꜀4 clones support, and so made the Sanger fold change look larger (41× against ≥ 11× here).
3. **The paralogue sequences were counted as unmutated J꜀4 clones,** which lowers every wild-type Sanger density. Removing them, counting each mouse's own clones and removing the three misfiled MKLN clones puts the wild-type geometric mean at 759 against 828.

## Were clones mixed between mice, or only misfiled?

The clones were collated by hand, and the folders mix mice. Names and folders disagree for 25 wild-type clones. Two kinds of evidence decide which is right, and neither uses a label.

**Mutation load.** WT3 is far more mutated than WT2 on both platforms, so a clone's load says which mouse it can belong to:

| Clones | n | Mutations per 10³ bp | Like | Unlike |
|---|---:|---:|---|---|
| named WT3, filed in the W2 folder | 11 | 26.4 | named WT3 in its own folder (28.5; p = 0.25) | WT2 (7.5; p = 0.0008) |
| named WT2, filed in the W3 folder | 11 | 1.6 | named WT2 in its own folder (p = 0.44) | WT3 (p = 0.0001) |

So for those, **the name is right and only the folder was wrong.**

**Shared lineages.** Clones sharing five or more identical mutations descend from one B cell, so they must come from one mouse. Nine such lineages exist among all 318 J꜀4 clones. Seven stay within the mouse their names give, including a 25-mutation WT3 clone found in both the W2 folder (9 Aug) and the W3 folder (29 Jul), which confirms that the WT3-named clones in the W2 folder are WT3's. Two cross mice:

| Lineage | Shared mutations | Clones | Verdict |
|---|---:|---|---|
| A | 8 | MKLN K4: k4-7, k4-22 (two records), k4-11, sequenced 29 Jul – 1 Aug 2024; **"WT2": w2-28, w2-36**, 9 Aug | K4 clones misfiled as WT2 |
| B | 7, identical | MKLN K2: k2-24, k2-1; **"WT4": w4-9**, same day as k2-1 | a K2 clone misfiled as WT4 |

Neither lineage is in the wild-type mouse it was filed under: no WT2 nanopore read carries two of lineage A's eight mutations, and no WT4 read carries two of lineage B's seven. Each lineage otherwise appears only in MKLN clones sequenced on several dates. **These three are real mix-ups, not naming slips**: a tube or colony from an MKLN mouse ended up among the wild types. They are removed from the wild-type sets; `~/Nanoshm/sanger_jh4_corrected/CORRECTIONS.txt` records the change, and the original files are untouched.

Removing them changes WT2 the most. Its Sanger load halves (7.5 → 3.7 × 10⁻³) because the two clones carried 22 and 24 mutations, and its unique mutations fall from 67 to 37. The spectrum agreement does not change (r = 0.92).

Limits: the Sanger reads start at J꜀4, so they carry no V(D)J junction, the strongest clone barcode. A clone that is a misfiled copy of a lightly mutated clone from another mouse, sharing no lineage, cannot be caught this way. Mutation load still catches a misfile between mice as different as WT3 and WT2.

## A bug this comparison found in nanoshm

Assigned correctly, A3 is a single clone, and it came out with 5 mutations: 945 × 10⁻⁵, higher than one wild-type mouse. These five changes are germline. They are carried by every Sanger clone of every mouse and by about 95% of all nanopore reads, knockouts included. They are the five positions the germline rule removes everywhere else.

A Sanger sample of fewer than 4 clones cannot apply the germline rule to itself, so nanoshm switched the rule off for it. **Now a sample with too few clones uses the germline changes measured across the project's other clones** (`project.sanger_germline`; regression test in `tests/test_platforms.py`). A3 now has no mutations.

## What does not agree, and why that is expected

- **Absolute density.** Sanger counts unique mutations over 18–34 clones and nanopore over 497 reads. More sequences resample the same mutations, so the per-sequence density falls with depth. Compare within a platform, not between them.
- **Mutation load per sequence** (Figure 3E, Figure 3D here). Nanopore per read is lower than Sanger per clone in every mouse. In WT2 it is 1.7 against 3.7 × 10⁻³, because one almost unmutated clone holds 85% of WT2's nanopore reads ([`CLONES.md`](CLONES.md)). Counted per clone, nanopore WT2 is 12.0. Separate PCRs, and a mouse whose response is dominated by a few lineages, make WT2 the least stable mouse on either platform. Separate PCRs can sample a lineage that dominates one reaction and not the other.
- **Per-mouse spectra against each other** (Figure 3A here). Each Sanger mouse matches its own nanopore spectrum no better than another mouse's (mean r 0.57 against 0.71). The spectrum is AID's signature, not the mouse's. WT2's and WT4's Sanger spectra rest on 37 and 68 mutations, too few for a 12-type profile per mouse; WT3's, with 290, correlates at r = 0.93 with its own nanopore.

## Reproduce

```bash
python docs/calibration/sanger_vs_nanopore.py --ont-mode recommended
python docs/calibration/sanger_vs_nanopore.py --ont-mode geneious
```

The Sanger FASTAs are read from `~/Nanoshm/sanger_jh4_corrected` (`--sanger-dir`): one file per mouse, built from the GenBank clone exports by clone name, with the three misfiled clones removed. The repeated-measures ANOVA assumes sphericity. With three mice, nothing here is powered to show equivalence; agreement is descriptive, as in the manuscript.
