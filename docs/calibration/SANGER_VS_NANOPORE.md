# Do Sanger and nanopore agree in the same wild-type mice?

*12 October 2026. Sanger: the lab's J꜀4 clones from their `.ab1` traces, one folder per mouse, every clone read from its trace. Nanopore: the six bundled libraries at 497 reads, in the recommended mode and under the Geneious rule. Window 101–629 (529 bp) on both platforms.*

## Answer

**Yes, on everything both platforms can measure.** They agree on:

- the substitution spectrum;
- where along the intron mutations fall;
- hotspot preference;
- which mouse is most mutated.

They also agree on something a mutation count alone does not show. **WT2 is dominated by one unmutated lineage, on both platforms**: 30 of its 36 Sanger clones carry no mutation, and one nearly unmutated clone holds 85% of its nanopore reads.

| Test | Recommended mode | Geneious rule |
|---|---|---|
| Unique mutations, Sanger / ONT | 415 / 548 | 415 / 585 |
| 12-type spectrum, Pearson r [95% CI] | **0.90 [0.66–0.97]** | **0.93 [0.76–0.98]** |
| Slope (ONT on Sanger) | 0.98 | 1.02 |
| Mean / largest difference per type | 1.7 / 3.8 points | 1.5 / 3.0 points |
| Platform × type interaction (RM ANOVA) | p = 0.18 | p = 0.20 |
| Smallest paired, Šídák-adjusted type difference | p = 0.49 (G>A) | p = 0.67 |
| G:C transitions / transversions, adjusted p | 0.19 / 1.00 | 0.17 / 0.95 |
| Position along the window, 25-bp bins, r | **0.89** (p = 6 × 10⁻⁸) | 0.89 |
| G:C mutations in WRC/GYW (window share 20%) | Sanger 58%, ONT 56% | Sanger 58%, ONT 55% |
| A:T mutations in WA/TW (window share 50%) | Sanger 80%, ONT 70% | Sanger 80%, ONT 69% |
| Mouse ranking by load | WT3 > WT4 > WT2 on both | same |

WT2's Sanger spectrum rests on 6 mutations, three of them A>T, so its per-mouse proportions are noise (50% A>T). A mean of the three mice would let that one mouse pull Sanger's A>T bar to 24%. Panels B and D therefore show each platform's pooled proportion over all its unique mutations, with a 95% Wilson interval, and every mouse as a point sized by its number of mutations. Pooled, Sanger A>T is 11.3% against nanopore's 8.9%. The paired per-mouse tests and the interaction are reported as before, and are looser than the pooled correlation for the same reason.

## The Sanger data

447 traces in eight sequencing orders: 278 usable J꜀4 clones from wild-type, *Aicda*−/− and MKLN mice, 71 noisy or failed reads, and 82 from the J꜀4 paralogue (Agouti-type) or another contaminant. Clone labels were reused between orders ("w4-2" in three orders is three different colonies), so every trace is one colony, identified by its order and well.

**Assigning clones to mice.** A clone stays with the mouse its label gives unless its sequence clearly says otherwise. "Clearly" is calibrated against clones from other mice: MKLN clones share at most 7 identical mutations with any single wild-type nanopore read by coincidence, so a clone moves only if a single read of another mouse's library shares 8 or more of its mutations. Somatic mutations are not inherited, and the strain's shared germline differences are removed first. Only a trace SHMtrace cannot read is excluded (no stretch of good quality, or no or too short an alignment); a clean but short trace stays, and the density counts the bases each clone read. A clone with AID-type mutations in an *Aicda*−/− sample is set aside as a contaminant.

Clones labelled "jh6" are J꜀4-intron clones: they align across the whole window with only the strain's five germline differences, like the other knockout clones, and the tag after "jh6" names the mouse.

| Mouse | Usable J꜀4 clones | Unique mutations |
|---|---:|---:|
| WT2 | 36 | 6 |
| WT3 | 41 | 303 |
| WT4 | 26 | 106 |
| *Aicda*−/− A1 | 11 | 1 |
| *Aicda*−/− A2 | 10 | 1 |
| *Aicda*−/− A3 | 10 | 2 |

**WT2's low count is real.** WT2 is the least mutated mouse on both platforms, and its nanopore library is dominated by one nearly unmutated clone, so 6 unique mutations in 36 Sanger clones is what that mouse looks like (below).

## Is the 529-bp window defensible?

Yes, and the reason is biological as much as technical. The J꜀4 coding segment runs from reference 53 to 104 (`…CTCCTCAG`), and the intron starts at 105 with the splice donor `GTAAG`. Positions 101–629 are therefore the J꜀4 intron, the region the assay measures because it is not translated and so not shaped by antigen selection, plus four bases of coding sequence. Its 3' end sits where the forward Sanger reads begin to fade.

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

**Sanger alone separates wild type from *Aicda*−/− 15-fold but not significantly.** All three knockouts now carry one or two mutations in 10–11 clones (17–38 × 10⁻⁵). WT2's own clones carry 6 mutations in 36 clones (32 × 10⁻⁵), close to that level.

Nanopore separates the groups clearly, because 497 reads also catch WT2's minority mutated lineages (146 unique mutations). This is the argument for depth.

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
| WT3 | 33 | 15 | all assigned WT3, apart from one AID-type contaminant in an *Aicda*−/− sample, which carries all 37 mutations of a 21-read WT3 lineage |
| WT2 | 12 | 2 | one is a five-mutation MKLN K4 lineage in 4 WT2 reads, most likely carry-over between samples; one shares 4 of 16 with a WT3 clone, within coincidence |
| WT4 | 3 | 0 | — |

None of WT2's mutated lineages is in a clone assigned to another mouse, so WT2's mutated clones were not filed as WT3's.

**Why WT2's Sanger count is low.** Resampling WT2's nanopore reads as if Sanger-cloning them separates two things:

- *How many sequences are mutated agrees.* 5 of 36 Sanger clones carry a mutation (14%), against 17% of nanopore reads (P = 0.40). On both platforms about 85% of WT2's sequences are one unmutated lineage, seen independently in separate PCRs.
- *How mutated they are does not.* The five carry 1, 1, 1, 1 and 3 mutations. Five reads drawn from WT2's mutated nanopore reads would carry about 23 unique mutations (95% range 7–51). With only about five non-dominant clones, Sanger sampled lightly mutated ones and missed the heavily mutated minority lineages. That is unlucky (36 random reads give a median of 28; P = 0.04 for 6 or fewer), not impossible.

So WT2's low Sanger value comes from a mouse dominated by one unmutated clone, sampled with too few clones to reach its minority lineages. A per-clone measure, which nanopore's junctions allow (WT2 12.0 × 10⁻³ per clone against 1.7 per read), is the fair way to compare this mouse with the others.

## A bug this comparison found in SHMtrace

A Sanger sample of fewer than four clones cannot apply the germline rule to itself, so SHMtrace switched the rule off. A3's single clone then reported the strain's five germline differences from the reference as five mutations. Now a sample with too few clones uses the germline changes measured across the project's other clones (`project.sanger_germline`; test in `tests/test_platforms.py`).

## Running it

Each mouse's folder runs in the app. Drop it (or choose it) as one sample, set the platform to Sanger, and the `.ab1` traces and the GenBank-only FASTA in it are read together.

## Limits

- **The Sanger reads start at J꜀4**, so they carry no V(D)J junction, the strongest identity check. A misfiled clone with few mutations and no relatives in its true mouse cannot be caught; mutation load still catches a misfile between mice as different as WT2 and WT3.
- **Three mice per group.** The agreement tests are descriptive, and the repeated-measures ANOVA assumes sphericity.
