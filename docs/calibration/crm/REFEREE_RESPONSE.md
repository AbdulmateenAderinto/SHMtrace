# Response to the external referee review of 2026-10-18

Source: `SHMtrace_external_referee_2026-10-18/REFEREE_REVIEW.md`.

Every claim below was checked against the actual source or output before acting. Where the referee was
right, the fix is applied and named. Where a claim could not be checked from what is in the repository,
it is listed under "needs the author" rather than being quietly acted on or quietly dropped.

The manuscript rebuilds: summary 145 words, main text 6,993 (cap 7,000), 52 references.

---

## 1. Confirmed and fixed

### Errors of fact in the text (several were ours, not the referee's misreading)

| # | What was wrong | What it says now | Verified against |
|---|---|---|---|
| 1 | Methods said the TCR reference was "verified against the project's own J_H4 reference". It is not: `tcr_prepare.py` does no such check, and J_H4 could not verify a TCR construct in any case. | States that the reconstructed length matches the SAM header's own declaration and that no read disagreed at any covered base, and that identity to the physical construct was **not** independently verified because the deposit does not publish that sequence. | `tcr_prepare.py` |
| 2 | "the Ramos clones of Figure 7F" | Figure 7D. 7F is the mOrange reporter. | Figure 7 panel map |
| 3 | "The AID-null clones can measure is the error itself" (ungrammatical) | "The AID-null clones measure the observed mismatch background" | — |
| 4 | S16D legend said a junction "aligns equally well at more than one position"; Methods describe a 5-base tolerance. | "has another placement within five matched bases of its best one" | `crm_methods.py` ambiguity criterion |
| 5 | S16A: "reads per unique junction differ 1.1–23× between samples" conflates a ratio range with a fold difference. | "range from 1.1 to 23 across samples" | `switch_pacbio.py` output |
| 6 | Table 1 omitted all six newly added datasets; Summary and eTOC said "three platforms" although PacBio makes four. | Six rows added (PRJNA831666, GSE77788, SRP159221, SRP075898, figshare 12803501, Devilder); "four platforms" throughout. | Table 1, `crm_text.py` |
| 7 | **Figures S8D and S9B** axis read "% of mutations at G:C", but `gc()` divides by `sum(spectrum.values())` — the share of **all** mutations, which is what Results and Table S4 report. | "% of all mutations" on both panels. The legends were already correct. | `figS8_sanger_cohort2.py:57`, `figS9_sanger_mkln.py:40` |
| 8 | **Table S3** printed "≥ 141-fold" and "≥ 85-fold" while the table's own footnote disclaims that these are lower confidence bounds. The main text already used the plain folds. | Both "≥" removed; the footnote explains the one-event convention. | `build_crm.py` STATS, footnote |
| 9 | "lowering the read-weighted frequency by **up to -62%**" and "fell to as much as **-40% below** the truth" — double negatives, because `min(err)` is a signed shortfall. | `crm_numbers.py` now emits the **magnitude** of the shortfall (and prefixes "+" if a key is ever not a shortfall), so the prose reads "by up to 62%". | `crm_numbers.py:182–198`, verified against `out/mask_sim_jh4/lineage_sim.json` (rw: −40.2 at 60%, −61.9 at 80%; cw: −26.3, −33.8) |
| 10 | Figure S6 panel and legend turned a **Spearman** coefficient into "% of variation explained". A rank statistic squared is not a share of variance. | Both now report the correlations: r = 0.12–0.20 over three pairs in J_H4, 0.29–0.45 in Sμ, described as weak rank agreement. | `figS6_background.py`, values printed by the rebuild |
| 11 | S13 legend: "checks … never applied" alongside "an accidental copy is counted once" — the second is false, because the QC pass changes nothing. | Flagged-not-removed is now separated from the authored sample-sheet exclusion, which is the only thing that drops a clone. | `project.py:805–810` ("Nothing was changed"), `qc.py` `sheet_excluded` |
| 12 | Methods invoked DESeq2 and edgeR as validation of the one-endpoint dispersion fit. | Analogy removed; the model was already fully specified, and it now says no cross-feature moderation applies with one endpoint. (References 53 → 52.) | `crm_methods.py` |
| 13 | "None of these is a curated population" — ambiguous antecedent; the TCR libraries **do** carry read-level frequencies, which S18D plots. | Scoped to the Devilder records, with the TCR frequencies stated and framed as endpoints of selection. | `crm_methods.py`, Figure S18D |
| 14 | "ClairS rejected 11 of its 12 candidates, **including neither** of the two carried by 20% or more" — says the opposite of what was meant. | "…; of the two carried by 20% or more of reads, a level at which it recalled nearly all clean planted mutations, neither passed." | `crm_numbers.py:486` |
| 16 | **The mOrange separation folds were wrong.** Results gave 105-, 8,102- and 2,033-fold. Two compounding faults: the denominator was floored at a hard-coded `1e-6` (so the constant, not the data, set the fold — `1e-5` would have given Ugi 810, `1e-7` 81,022), and `mutation_frequency_x1e3` is rounded to four decimals, so dead-AID libraries with real mutations stored as exactly `0.0` and were then floored. | Corrected to **20-, 1,002- and 103-fold**. UNG was out 20× from the rounding alone. Frequencies are now computed from counts; only a library that genuinely saw nothing is floored, at one event on its own denominator (the Figure 4F / Table S3 convention). Results states that only the Ugi arm rests on a floored denominator. | `out/public/morange/project.json`; `crm_numbers.py` |
| 21 | **Figure S18B drew no guide positions**, and the deposit does not publish them (referee minor 11). | Obtained the paper's supplementary. Supplementary Figure 9 draws each sgRNA on the ABL exon 6 sequence with its PAM marked but **does not number them**, so per-guide identity is still unavailable. The protospacers were instead located computationally: each PAM was found in the amplicon consensus rebuilt from the reads, and the 20 bases 5' of it taken as the protospacer. S18B now shades those nine + strand protospacers, unlabelled, as the source figure presents them. **All five libraries' highest peaks fall inside one** (#4 at 49 in 46–65; #6 at 95, #7 at 98 and #15 at 94 in 83–102; #9 at 154 in 146–165), which is the positional relationship the referee asked to see. | Ma et al. Supplementary Figure 9; `out/tam/footprint.json` |
| 20 | **Figure S18C explained the zero at three guides as "guides target outside the amplicons recovered here".** Ma et al.'s Supplementary Figure 12 says all four pools target the same 100-bp window, positions 201–300 of ABL exon 6, and that only seven guides fit inside it, so the ten-guide pool adds three that flank it. The guides are not aimed elsewhere. | Legend corrected to say so: the zero at three guides is too little signal to clear the five-fold threshold. The earlier wording offered a technical excuse the source paper's own design contradicts. | Ma et al. Supplementary Figure 12 |
| 19 | **Table S3 kept the spectrum's Pearson 95% CI** after Results dropped it, under an explicit "Pearson r [95% CI]" label. | Dropped there too, so the two agree. The column now reads "Pearson r (descriptive); slope". A Fisher interval over 12 compositional categories pooled from three mice describes the types, not the animals. | `build_crm.py` STATS |
| 18 | **526 bp vs 529 bp** (referee's numerical table). The group simulation's "500 reads × 526 bp" sat a few lines from the 529-bp analysis window with nothing distinguishing them, and the referee asked whether 526 was stale. | **Intentional, and now stated.** 526 is the mean number of bases a read covers, not the window: reads do not all span every position. Methods now gives the denominator as 262,913 base observations and says so; the J_H4 libraries average 525 covered bases per read (range 516–528), so the simulation sits inside the observed spread. Caveat left for the author: 262,913 matches no single library (nearest KO3 at 263,019) nor the exact mean (262,393), and its derivation is not recorded. | `docs/calibration/replicates_calibration.py`, one level above `crm/`, which is why an earlier search missed it |
| 17 | Results said "Calibrated against the dead-AID libraries of the same experiments", which the referee read as each arm calibrating itself. | "Calibrated against **that arm's** dead-AID libraries" — only the empty arm's three are controls, which is what the code does. Settles major issue 2's first half. | `project.json` control flags |
| 15 | Resource availability claimed files "**have been** deposited at [data deposit DOI]" with the identifier unresolved. | "will be deposited … and will be publicly available as of the date of publication". `lab_details.py` now prints a reminder to switch the tense once the identifiers resolve. | build checklist |

### Overclaims bounded

- **Introduction.** "Neither approach accounts for the errors that a particular library makes, so sequencing
  errors … **are** wrongly counted" → "In the configurations evaluated here, neither measures the errors a
  particular library makes, so mismatches that recur at one position … **can** be counted as SHM mutations."
- **Discussion.** "a per-base error model **cannot know** an SHM assay's library-specific background" →
  "a per-base error model fitted without matched controls **had no way to reach** this assay's library-specific
  background."
- **Figure 2E lead-in.** "SHM rarely reaches the allele fractions where these callers recovered planted
  mutations" → "The calls SHMtrace **retained** rarely reach such fractions", since the distribution is
  conditional on detection and on the <50% rule.
- **Spectrum correlation.** The Fisher interval is gone from Results: "the spectra pooled across the three mice
  were similar over the 12 substitution types (descriptive Pearson r = 0.92; slope = 0.97)". A Fisher interval
  over compositional categories does not carry a population-confidence reading.
- **Sanger provenance.** Results now gives the analysis n: "the 184 of their 447 archived chromatograms that
  the sample sheet files to these mice, 134 of which yielded clones". The full flow stays in the S7 legend.

### Figures redrawn

- **Figure 1C.** "FDR ≤ 0.05" → "empirical FDR target 0.05" in both places. The germline box now reads
  "window first, then alleles in ≥ 50% of reads set aside and reported, not dropped, before any test; indels
  optional (opt-in: known ancestor)". The ordering claim was checked: `variants.py` applies the window, the
  end trim and the frequency rule and only then computes a p-value, so the diagram now states the implemented
  order. (The hypothesis family is positions × 3 regardless, so a set-aside site does not change it.)
- **Figure 2D.** Red curve renamed "decoy-estimated false calls". The log x-axis was running to 10 although a
  p-value threshold cannot exceed 1; now capped at 1.
- **Figure S17A.** Panel letter moved to the upper-left of the schematic (it sat below everything).
- **Figure S17C.** The panel now states one wild-type library per treatment against three AID-null, pooled as a
  descriptive contrast; the legend no longer calls treatment a "nuisance factor". Counts confirmed from
  `out/public_smu/public_smu.json`: WT 1 per treatment, AID-null 3 per treatment.
- **Figure S18B.** The title claimed "The position follows the guide" — the deposit names the sgRNAs but does
  not publish their protospacer coordinates, so the panel cannot show that relationship. Retitled "Each
  library's excess peaks at its own position", and the caption says guides are not drawn and why.
- **Figure S18C.** Retitled "Localised sites in four guide-pool libraries"; the joining line is gone (it read
  as a dose-response across four unreplicated libraries); the annotation says descriptive contrast.
- **Figure S18D.** The pooled mean bar is replaced by per-founder groups (h5, h12, h7, h8) with the baseline
  labelled as h5's, so the h5 descendants' different relationship to the baseline is visible.

### Checked and found to be correct as written

- **mOrange "positions 75–693" against a 679-bp amplicon.** Not an error. The project's analysis window is
  [31, 649] in the 679-bp reference rebuilt here, which begins at position 45 of the study's numbering — a
  constant offset of 44 (75−31 = 693−649 = 44). The mapping is now stated in Methods, as the referee asked.
  Verified in `out/public/morange/project.json`.
- **The S7 flow (447 / 184 / 182 / 134) and the arithmetic the referee reconciled.** Unchanged; only its
  visibility in Results changed.

---

## 2. A floor-dependent significance claim at OT1 and OT2, found in the figure audit and now fixed

Not raised by the referee. Results said the MLLT1/MLLT3 double knockout sat "below the AID-competent
libraries (p = 0.034 and 2.6 × 10⁻⁵)". Both came from a Welch test on `log10(mutation_frequency_x1e3)`
with frequencies hand-floored at `1e-4`, and the floor decided them:

| floor | OT1 p | OT2 p |
|---|---|---|
| 1e-3 | 0.0085 | 5.6e-05 |
| **1e-4 (used)** | **0.034** | **2.6e-05** |
| 1e-5 | 0.070 | 1.5e-05 |
| 1e-6 | 0.106 | 9.8e-06 |

OT1 crossed the significance threshold on the constant alone, and OT2 was degenerate: every dKO library saw
nothing, so after flooring that group had zero variance and SciPy warned of catastrophic cancellation.
The comment directly above that code already rejected folds for this reason, then used a floored log for a test.

**Fixed.** Both p values are gone from Results, which now reads "neither separable from the controls
(p = 0.45 and 0.56) and both far below the AID-competent counts above". The counts are unambiguous and are
stated in the preceding sentence: 1, 1 and 0 calls at OT1 and 0, 0 and 0 at OT2, against 28, 23 and 27 at OT2
and comparable numbers at OT1. The floored computation was deleted from `crm_numbers.py` so it cannot be
reused, with a comment recording why.

**The pipeline's own log test was checked and kept.** `_welch_log` censors a zero to a third of the smallest
observed value — data-derived, documented as pulling towards no difference, and stable: OT1's reported
p = 0.035 comes out at 0.0347, and 0.030 if a half is used instead. That is a different thing from a constant
chosen independently of the data, so the remaining p values stand.

## 3. Figure and table audit, 18 October

Every main figure, every table and every scripted supplementary figure was re-derived from its source files.
Two defects were found and fixed (the mOrange folds, section 1 item 16; the OT1/OT2 p values, section 2) and
one stale count (Table 1's ten loci / seven laboratories). Everything else reproduced exactly:

| | checked | result |
|---|---|---|
| Figures 1-7 | every annotated number | exact |
| Table S1 | 12 libraries | read-fate categories sum to `delivered`; cross-validates Figure 4G |
| Table S2 | 275 real-library + 52 spike-in rows | exact |
| Tables S3-S5 | computed rows | exact |
| Figures S2, S3, S5-S10, S15-S18 | annotated numbers | exact |

Worth recording from that pass:

- **Figure S10 cross-validates a Methods claim.** Its sequences sum to 1,795 for SHMtrace against 1,796
  published, the single Msh2 record the Methods describes, and its mutations to 2,804 of the 2,809 the paper
  counts.
- **Figure S15** reproduces in all 16 conditions, including the 1/25 all-null discovery rate and the 100%
  set-aside under one noisy control.
- **Figure S16** confirms the known ambiguity defect: 11 of 12 libraries report exactly 0.0% and the twelfth
  0.2%, against 19-47% on direct re-alignment. Still deliberately unpatched; see major issue 5.
- **Figure S5's** one above-100% value is legitimate: it is a 500-read draw against that library's *mean* at
  1,000 reads, and the axis says so.

The numbers in the text cannot drift from their sources, because `crm_numbers.py` computes them from the JSON
at build time. The risk is in the derivation, not the transcription, which is where both defects were. A scan
of every script for the two patterns involved - a constant floor, and a falsy real value silently replaced -
leaves only two hits, both cosmetic and feeding no reported number (a log-axis start in `fig2_call.py`, a
confidence-band clip in `figS2_quality.py`).

**STAR Protocols**, audited last, shares `crm_numbers` with the main manuscript, so its 28 placeholders
cannot drift and it hard-codes no result values. Its prose had drifted anyway, because it was last built on
16 October and not re-read after the main paper was corrected. One sentence contradicted the manuscript
outright: STAR said controls without calls are placed at their one-event density "so fold changes are lower
bounds", while Table S3's footnote says such folds "are not statistical lower confidence bounds". Fixing it
exposed a second, internal to STAR: Expected Outcomes told users to look for "at least 141-fold" eighteen
lines earlier. Both corrected, and the deposit statement matched to the main paper's future tense.

## 4. Needs the author (could not be settled from the repository)

Nothing outstanding in this section. The three items that sat here have all been settled: the 526 bp
simulation length (section 1, item 18), Table S3's Pearson interval (item 19), and the Figure S18B guide
positions (item 21). What remains for the author is section 5.

---

## 5. The nine major issues

Six are done, one needs an authorial decision, two remain.

| # | issue | status |
|---|---|---|
| 1 | Benchmark vs production route | **Done.** `crm/route_check.py` runs both over all seven datasets. J_H4 and N114 agree exactly (12/12); the four public sets differ, because the leave-one-control-out margin is fitted over a different set. Every separation conclusion survives, and *Pim1* weakens from 1.8 to 0.8, strengthening the limit-of-detection claim. In Methods. |
| 2 | mOrange control fitting | **Done.** Only the empty arm's three dead-AID libraries are controls, and Exp3 is set aside as noisy; Results now says "that arm's dead-AID libraries". Settled with the fold correction (section 1, item 16). |
| 3 | Mode matrix | **Done.** Methods now sets out all four control counts: none, one, two, three or more, with what each does to the background, the margin, the correction, indels and the noisy-control rule, and notes that the output's `multiple_testing` field records which was used. The TAM reporter arm is named as the one-control case whose calibration inverted. |
| 4 | Clone-weighted estimand | **Needs the author.** The referee proposes "lineage-consensus burden". Mechanical to apply but it is the paper's vocabulary, across Results, Discussion, Methods, two figures and a table. Not chosen unilaterally. |
| 5 | CSR ambiguity | **Done.** `crm/switch_ambiguity_sensitivity.py` drops every ambiguous junction, 2–21% of each sample, and re-runs the comparisons. Thirteen of fifteen hold; two move in opposite directions and are flagged provisional in Methods and the S16 legend. |
| 6 | TCR and Devilder variant recovery | **Done for Devilder, impossible for the TCR.** All 150 published mutations across the 42 selected clusters are recovered, none missed and none added. No equivalent list is published for the TCR libraries, whose deposit carries alignments only; Methods says so. |
| 7 | Effect sizes and multiplicity | **Done.** The fifteen S16C comparisons are unadjusted and reported with both group means; Benjamini–Hochberg across them changes nothing, nine significant before and after. Figure S18A is a single pre-specified contrast, so no adjustment arises. In Methods. |
| 8 | Table 1 and the inventory | **Done.** Fourteen AID-targeted loci from thirteen laboratories, recounted on the original convention and split by class (section 1, item 20). |
| 9 | STAR Methods completeness | **Partly.** Equations, versions and the route reconciliation are in. What remains is the deposit identifiers and the bench details, which are the 14 items in `NEXT_STEPS.md` and need the author. |

