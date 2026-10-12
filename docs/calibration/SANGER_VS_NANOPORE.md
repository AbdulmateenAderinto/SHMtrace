# Do Sanger and nanopore agree in the same wild-type mice?

*12 October 2026. Sanger: the lab's J꜀4 clones from their `.ab1` traces, organised by mouse with `organize_jh4_ab1.py`, every clone now read from its trace. Nanopore: the six bundled libraries at 497 reads, in the recommended mode and, for comparison with the manuscript, the Geneious rule. Window 101–629 (529 bp) on both platforms. Figures: `figures/sanger_vs_nanopore/<mode>/`; script: `sanger_vs_nanopore.py`.*

## Answer

**Yes, on everything both platforms can measure.** They agree on:

- the substitution spectrum;
- where along the intron mutations fall;
- hotspot preference;
- which mouse is most mutated.

They also agree on something the manuscript could not see. **WT2 is dominated by one unmutated lineage, on both platforms**: 30 of its 36 Sanger clones carry no mutation, and one nearly unmutated clone holds 85% of its nanopore reads.

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

WT2's Sanger spectrum rests on 6 mutations, three of them A>T, so its per-mouse proportions are noise (50% A>T). A mean of the three mice would let that one mouse pull Sanger's A>T bar to 24%. Panels B and D therefore show each platform's pooled proportion over all its unique mutations, with a 95% Wilson interval, and every mouse as a point sized by its number of mutations. Pooled, Sanger A>T is 11.3% against nanopore's 8.9%. The paired per-mouse tests and the interaction are reported as before, and are looser than the pooled correlation for the same reason.

## The Sanger data

The traces are in eight sequencing orders, 447 files. The earliest, `406809` (11 October 2023), holds the first J꜀4 clone of each mouse (`8-w2-jh4`, `10-w4-jh4`, `15-a1-jh4`, …). Two orders were run on 9 January 2024: `414310` (the `jhextra-w3` and `jhextra-w4` clones) and `414309` (mostly clones labelled "jh6", plus K1 and w1 "extra" clones). The other five were run in July 2024: `432625`, `434138`, `434141`, `434355` and `434545`. Each trace stores the order number, not the clone name, so a clone's identity is its file name and its plate well. Geneious cut the GenBank record names at 23 characters (`23-jhextra-w3-2-user_ad` is the trace `23-jhextra-w3-2-user_added`), which the organiser allows for.

The January traces confirm their GenBank records. All 29 that have one carry exactly the same mutations as their record, and one (`jhextra-w4-10`, 5 mutations) is a clone the GenBank set never had.

| What | Count |
|---|---:|
| Traces | 447 |
| Usable J꜀4 traces, wild type, *Aicda*−/− and MKLN | 278 |
| Noisy or failed reads: no stretch of good quality, or no or too short an alignment | 71 |
| J꜀4 paralogue (Agouti-type) or other contaminant | 82 |
| Wild-type and *Aicda*−/− GenBank records without a trace | 0 |

Clone labels were reused between orders: "w4-2" in three orders is three different colonies (3, 31 and 0 mutations, none shared). Every trace is therefore one colony, identified by order and well.

The GenBank records are faithful to the traces. Inside the window, 87 of 95 differences are a low-quality (Q < 20) or N call that was corrected by hand. One record, `3-w2-3`, carries four C's in 9 bp (591–599) where the trace reads the reference at Q23–27. The trace wins, and that trace fails alignment anyway.

## Which clone belongs to which mouse

**The rule.** A clone belongs to the mouse the original GenBank organisation gives it, unless the sequence clearly shows otherwise. "Clearly" is calibrated against clones from *other* mice. MKLN clones, compared with the wild-type nanopore reads, share at most 7 identical mutations with a single read by coincidence: the J꜀4 intron has strong hotspots, and WT3's reads carry about 11 mutations each. Somatic mutations are not inherited, so inbred mice do not share them; the germline differences they do share are removed before any of this. A clone therefore moves only if a single read of another mouse's nanopore library shares **8 or more** of its mutations, or a whole batch it belongs to is shown to be another mouse's.

**The default is the GenBank folder** the student filed each clone in (W2, W3, W4, A1, A2, A3); a trace with no GenBank record falls back on its label. **Only a trace nanoshm cannot read is excluded**: no stretch of good quality (a noisy chromatogram), or no or too short an alignment. A clean but short trace stays, and the density counts the bases each clone actually read.

**Clones labelled "jh6" are J꜀4-intron clones.** They align across the whole window with only the strain's five germline differences, exactly like the a1/a2/a3 clones; an earlier version of this analysis excluded them by label as "another amplicon", which was wrong. The tag after "jh6" names the mouse (`jh6-a3-g9` is an A3 clone). The manuscript's single *Aicda*−/− 3 mutation is in `81-jh6-a3-g9`.

**What the rule gives.** The GenBank folders and the clone names disagree for 30 wild-type clones; the names were written on the plates, the folders were made by hand later.

- Seven clones filed under W2 clear the bar for WT3: five labelled w3/wt3, plus two labelled w2 (E8/F8, below). Each shares 8–13 mutations with one WT3 read and at most 3 with any WT2 read. **No clone clears it for WT2.**
- Nine more clones filed under W2 but labelled w3/wt3 carry 1–26 mutations without a decisive single read. They are WT3's as a group. All 22 clones both named and filed as WT2 carry 0 or 1 mutations; these nine differ from them beyond chance (p = 1.5 × 10⁻⁶) and not from WT3's clones (p = 0.28). They also sit in the same plate blocks as the five proven WT3 clones (434138 B11–F11, 414310 D4–H4).
- Twelve clones filed under W3 but labelled w2 are near-unmutated (eleven carry no mutation). Only 1 of WT3's 23 clones is unmutated, so they cannot be WT3's; they are WT2's, as their names say.

So the GenBank default plus clear evidence gives the same wild-type grouping as the clone names, except for the traces below. The knockouts follow their folders: A3's folder holds `15-a1-jh4`, `16-a2-jh4` and `17-a3-jh4`, and they stay in A3.

**Where your Geneious WT2 value came from.** In the manuscript's W2 assembly, the 189 variants passing p ≤ 0.05 and frequency < 50% are carried by: clones named WT3, 124; the two renamed K4 records, 24; E8/F8, 23; combinations of those, 16; WT2's own 16 clones, 1–2. In the W3 assembly, the 11 clones named WT2 carry one variant of their own. WT2's own clones were nearly unmutated in the original analysis too; the other mice's clones in its folder hid it.

**The Agouti-type clones do not change this.** Aligned to the Agouti reference with its own background removed, the 19 from wild-type mice (1 WT2, 4 WT3, 14 WT4) carry no mutations, and none of them is in the GenBank folders. Including them would only add unmutated clones. Applied blindly, the folder default would put a WT3-like population into WT2 and raise its Sanger frequency to 6.5 × 10⁻³, three times its own nanopore value.

**WT2's low count is real.** WT2 is the least mutated mouse on both platforms, and its nanopore library is dominated by one nearly unmutated clone, so 6 unique mutations in 34 Sanger clones is what that mouse looks like.

Five traces contradict their labels:

| Trace | Label | Assigned | Evidence |
|---|---|---|---|
| 434141 E8 | w2-1 | **WT3** | Its lineage is in WT3's nanopore library and not WT2's: single WT3 reads carry 13 of its 18 mutations (12 reads share three or more), and no WT2 read shares more than two. It also sits where w3-1 belongs on the plate (G8 w3-3, H8 w3-4, A9 w3-5), and w2-1 is already at A7. Filed under W2 by hand, so the manuscript counted it as WT2. |
| 434141 F8 | w2-2 | **WT3** | Same: WT3 reads carry 10 of its 13 mutations (35 reads share three or more); no WT2 read shares more than two. w2-2 is already at B7 on the plate. |
| 434141 A11 | w4-9 | **MKLN K2** | Carries the same 7-mutation lineage as k2-1 (E11, same plate, where the K2 series starts right after the WT4 column) and k2-24 (order 434355). Not in any WT4 nanopore read. |
| 414309 F7 | jh6-a1-e6 | **removed** | 5 confident AID-type mutations (Q26–56; four at WRC/GYW or WA/TW, mostly A:T transitions) in an *Aicda*−/− sample, with no GenBank record: a contaminant, like a1-1 below on a smaller scale. |
| 432625 E6 | a1-1 | **unassigned** | 37 AID-type mutations (64% of G:C changes at WRC/GYW), impossible for an *Aicda*−/− mouse. Not the paralogue: 50 differences from J꜀4, 93 from the paralogue. It shares no lineage with any clone, so its mouse is unknown. It had no GenBank record, so the manuscript did not use it. |

The GenBank records `88-w2-28` and `96-w2-36` (filed as WT2, 22 and 24 mutations) are **renamed copies of K4 traces** 434355 H11 `88-k4-28` and H12 `96-k4-36`: identical sequences, same well-order numbers, same clone numbers. They are excluded.

A caution from the same test: five K4 clones carry exactly the same five changes and nothing else, one expanded K4 lineage, and 4 of WT2's 500 nanopore reads (0.8%) carry all five. No other mouse's reads do. A germline variant would appear across K4's clones and in about half of a heterozygous mouse's reads, so this is most likely low-level carry-over between samples. It does not affect the Sanger grouping, and at 4 reads it sits at the edge of what nanopore calls in WT2.

The organised set (`~/Nanoshm/sanger_ab1_by_mouse`, with `manifest.csv` and `README.md`):

| Mouse | Traces | GenBank-only clones | Usable J꜀4 clones | Unique mutations |
|---|---:|---:|---:|---:|
| WT2 | 36 (13 moved from the W3 folder) | 0 | 36 | 6 |
| WT3 | 41 (16 moved from the W2 folder) | 0 | 41 | 303 |
| WT4 | 26 | 0 | 26 | 106 |
| *Aicda*−/− A1 | 11 | 0 | 11 | 1 |
| *Aicda*−/− A2 | 10 | 0 | 10 | 1 |
| *Aicda*−/− A3 | 10 | 0 | 10 | 2 |

## Is the 529-bp window defensible?

Yes, and the reason is biological as much as technical (`figS_sanger_window`, `sanger_window_figure.py`). The J꜀4 coding segment runs from reference 53 to 104 (`…CTCCTCAG`), and the intron starts at 105 with the splice donor `GTAAG`. Positions 101–629 are therefore the J꜀4 intron, the region the assay measures because it is not translated and so not shaped by antigen selection, plus four bases of coding sequence. Its 3' end sits where the forward Sanger reads begin to fade.

| Window | Clones spanning it end to end | Unique wild-type mutations | What the extra bases are |
|---|---:|---:|---|
| **101–629 (529 bp)** | most | 415 | the intron |
| 74–629 (556 bp) | the same | 425 | the last 27 bases of the J꜀4 coding segment |
| 101–656 (556 bp) | about 10% fewer | 429 | intron where the reads are ending; none reach 670 |

The background is not what separates them. Sanger error in the *Aicda*−/− clones is 3 mismatches in 7,375 bases inside the window and none in either flank.

- **Against 101–656**, 529 is cleaner in the sense that matters for a frequency: F = mutations ÷ (sequences × length) assumes every sequence covers every base, and nine more clones fall short of 656 than of 629.
- **Against 74–629**, read quality is the same, so "cleaner" is the wrong word. The defensible reason is that 74–100 is coding sequence, which in a productively rearranged allele is part of the antibody (FR4) and can be selected; wild-type mutation density there is less than half the intron's (5.6 against 13.9 per 10³ bp).

Either way the choice hardly moves the numbers: each mouse's count changes by at most 3% between the three windows. Some clean clones do not span every base, so mutations ÷ (n × 529) slightly understates the frequency. The densities here divide by the bases each clone actually read, which is the same as n × 529 for full-length reads.

## AID dependence on each platform

| | Wt geometric mean [95% bootstrap CI] | *Aicda*−/− | Fold | P (Welch, log) |
|---|---|---|---:|---:|
| ONT, recommended | 56.7 [26–127] × 10⁻⁵ | 0.55 | 103 | 0.0017 |
| ONT, Geneious rule | 60.7 [26–131] | 1.43 | 42 | 0.0060 |
| **Sanger** | 339 [32–1447] | 23 [17–38] | 15 | **0.15** |
| *Manuscript: ONT / Sanger* | *59 / 828* | *1.2 / 20* | *47 / 41* | *0.0037 / 0.0045* |

**With the clones assigned correctly, Sanger alone separates wild type from *Aicda*−/− 15-fold but not significantly.** All three knockouts now carry one or two mutations in 10–11 clones (17–38 × 10⁻⁵). WT2's own clones carry 6 mutations in 36 clones (32 × 10⁻⁵), close to that level. The manuscript's Sanger WT2 was inflated by WT3 and MKLN clones filed under W2, and its knockout denominators by JH6 clones.

Nanopore separates the groups clearly, because 497 reads also catch WT2's minority mutated lineages (146 unique mutations). This is the manuscript's own argument for depth, and the corrected data make it more strongly than the published figure did.

## WT2: the same expanded clone on both platforms

| | Sanger | Nanopore |
|---|---|---|
| Sequences with no mutation | 30 of 36 clones (83%) | the largest clone holds 85% of junction reads, 7% of them mutated ([`CLONES.md`](CLONES.md)) |
| Mutations per 10³ bp, per sequence | 0.4 | 1.7 |
| Per clone (nanopore junction) | — | 12.0 |

Two independent PCRs and two platforms see one almost unmutated lineage dominating this mouse's germinal-centre B cells. Sanger's 36 clones, almost all from that lineage, miss the rest; nanopore's 497 reads do not.

## Can the nanopore data say which clones are WT2's?

Yes, for any clone that carries enough mutations, because a clone's mutations travel together in its relatives' reads.

**Each mouse's nanopore lineages, searched for among all Sanger clones** (lineages = reads grouped by V(D)J junction, consensus of 4 or more mutations; a hit = 4 or more shared):

| Nanopore library | Mutated lineages | Found in a Sanger clone | Where those clones are |
|---|---:|---:|---|
| WT3 | 33 | 15 | all assigned WT3, including E8/F8, `wt3-6` and `wt3-2` from the W2 folder; and `a1-1`, the contaminant in the A1 sample, which carries all 37 mutations of a 21-read WT3 lineage |
| WT2 | 12 | 2 | one is the K4 five-mutation set (carry-over); one shares 4 of 16 with a WT3 clone, within coincidence |
| WT4 | 3 | 0 | — |

None of WT2's mutated lineages is in a clone assigned to another mouse, so WT2's mutated clones were not filed as WT3's.

**Why WT2's Sanger count is low.** Resampling WT2's nanopore reads as if Sanger-cloning them separates two things:

- *How many sequences are mutated agrees.* 5 of 36 Sanger clones carry a mutation (14%), against 17% of nanopore reads (P = 0.40). On both platforms about 85% of WT2's sequences are one unmutated lineage, seen independently in separate PCRs.
- *How mutated they are does not.* The five carry 1, 1, 1, 1 and 3 mutations. Five reads drawn from WT2's mutated nanopore reads would carry about 23 unique mutations (95% range 7–51). With only about five non-dominant clones, Sanger sampled lightly mutated ones and missed the heavily mutated minority lineages. That is unlucky (36 random reads give a median of 28; P = 0.04 for 6 or fewer), not impossible.

So WT2's low Sanger value comes from a mouse dominated by one unmutated clone, sampled with too few clones to reach its minority lineages. It does not come from misfiling. A per-clone measure, which nanopore's junctions allow (WT2 12.0 × 10⁻³ per clone against 1.7 per read), is the fair way to compare this mouse with the others.

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
python docs/calibration/sanger_inclusion_figure.py && python docs/calibration/sanger_window_figure.py
python docs/calibration/build_sanger_report.py      # this page as HTML (needs: pip install markdown)
```

The same folders run in the app. Drop a mouse's folder (or choose it) as one sample, set the platform to Sanger, and the `.ab1` traces and the GenBank-only FASTA in it are read together.

## Limits

- **The Sanger reads start at J꜀4**, so they carry no V(D)J junction, the strongest identity check. A misfiled clone with few mutations and no relatives in its true mouse cannot be caught; mutation load still catches a misfile between mice as different as WT2 and WT3.
- **Three mice per group.** The agreement tests are descriptive, as in the manuscript, and the repeated-measures ANOVA assumes sphericity.
