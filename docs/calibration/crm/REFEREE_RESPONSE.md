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

## 2. Needs the author (could not be settled from the repository)

1. **526 bp vs 529 bp (Methods, group-test simulation).** The script that produced "pooled reads 0.91–0.94;
   negative binomial 0.050–0.057; Welch 0.037–0.043" is not in `crm/`, `src/shmtrace/` or `tests/`. Either the
   526 is a stale length and should be 529, or the simulation deliberately used a trimmed window — the text
   should say which. **Not changed.**
2. **Guide, protospacer and PAM positions for Figure S18B** (referee minor 11). These are not in the SRP075898
   deposit as downloaded; they would come from the paper's supplementary material, which would mean a new
   download. **Not fetched.**
3. **Table S3's Pearson 95% CI for the spectrum.** The interval is gone from Results but remains in Table S3
   under an explicit "Pearson r [95% CI]" label. Decide whether to drop it there too.

---

## 3. Still open: the nine major issues

These need analysis or an authorial decision, and none has been acted on. In the referee's own priority order:

1. **Benchmark vs production route (major 1).** Which route produced each Figure 4 library and each headline
   Summary value. Needs a route-by-route reconciliation table.
2. **mOrange control fitting (major 2, question 5).** Which dead-AID library IDs fitted each repair arm, and
   whether Results matches what Methods describes.
3. **A mode matrix (major 3).** Explicit operating modes for zero, one, two and three controls, and what the
   product does in the TAM configuration where calibration inverted. Includes whether exported labels are still
   FDR labels in that mode.
4. **The clone-weighted estimand (major 4).** The referee proposes "lineage-consensus burden", on the grounds
   that the endpoint changes mutation **inclusion** as well as lineage **weights**. This is a rename plus a
   definition, applied consistently across Results, Discussion and Methods.
5. **CSR ambiguity (major 5).** The module reports 0% ambiguity on all 12 PacBio libraries while direct
   re-alignment measures 19–47%. This is already documented as a confirmed defect in `DATASET_LEDGER.md` and
   deliberately not patched; the referee wants its effect on the displayed endpoints tested.
6. **TCR variant recovery (major 6, question 9).** A table of which published TCR and Devilder variants are
   recovered and which are not.
7. **Effect sizes for the new endpoints (major 7).** Whether the S18A motif and S16C p values are adjusted
   across contrasts, and effect differences reported beside them.
8. **Table 1 and the inventory (major 8).** Partly done — the six rows and "four platforms" are in. The
   locus and laboratory totals still need a documented recount.
9. **STAR Methods completeness (major 9).** Equations, software versions, container digests, and the deposit
   statements. Overlaps the 14 submission blockers already tracked in `NEXT_STEPS.md`.

The referee's own summary is that none of these requires new experiments or new datasets — only a smaller set
of clearly defined claims supported by the analyses already run.
