# Response to the external review of 16 October 2026

Review: `~/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/nanoshm/output/reviews/SHMtrace_manuscript_review_2026-10-16.md`
Recommendation: major revision before submission.

Status key: **open** / **investigating** / **fixed** / **won't fix (with reason)**.

---

## Major issue 1 — OT1/OT2: two different analyses (fixed, c5ad350)

**Verified. The reviewer is right, and the cause is identified.**

The manuscript reads OT1/OT2 from two different pipelines:

| Consumer | Source file | Produced by | What it is |
|---|---|---|---|
| Results text (`crm_numbers.py:449-475`) | `~/Nanoshm/public_data/PRJNA1481264/offtarget_ot{1,2}_2000.json` | `run_ircm_offtargets.py` → `project.run_project` | The full method: decoy-based FDR threshold selection **and** noisy-control QC |
| Figure 4 (`figs/fig4_controls.py:30`) and Table S2 (`make_tables.py`) | `out/benchmark_real_ot{1,2}.json` | `benchmark_real.py` | A **caller-comparison harness**: plain leave-one-control-out null, no threshold selection, no control QC, so SHMtrace, Geneious and LoFreq see identical reads |

Both use the same window, depth 2,000, `recommended` mode and `qcal` error model. They differ in procedure, so they legitimately give different numbers. Neither is "wrong"; they answer different questions. The error is presenting them as one analysis.

Counts, confirming the reviewer's C01-C06:

| | OT1 WT | OT1 dKO | OT1 AID-null | OT2 WT | OT2 dKO | OT2 AID-null |
|---|---|---|---|---|---|---|
| Text (`run_project`) | 79/56/57 | 1/1/0 | 0/2/1 | 28/23/27 | 0/0/0 | 0/0/1 |
| Fig. 4 / Table S2 (`benchmark_real`) | 60/46/52 | 2/1/3 | 0/0/0 | 35/30/21 | 1/0/0 | 0/0/0 |

Why the AID-null controls differ: for OT1, `run_project` emitted the note *"Control KO2 is far noisier than the other controls and was left out of the decoys that set the false-discovery threshold (it still contributes to the error model)."* Excluding KO2 from the decoys loosens the threshold, which both raises the wild-type counts (79/56/57 vs 60/46/52) and lets KO2 keep 2 calls of its own. `benchmark_real.py` has no such exclusion, so its threshold is stricter and every control lands on 0. This ties issue 1 directly to issue 4.

**Second, independent defect — Figure 4F pools the perturbation with wild type. Confirmed.**

`separation()` (`figs/fig4_controls.py:42`) defines the numerator as *every non-control library*. For OT1/OT2 that is wild type **plus** the MLLT1/MLLT3 double knockout. Recomputed from `benchmark_real_ot*.json`:

| | OT1 | OT2 |
|---|---:|---:|
| Wild type / AID-null | 52.78 | 28.18 |
| Double knockout / AID-null | 1.83 | 1.00 |
| **(WT + dKO) / AID-null — what Fig. 4F plots** | **9.83** | **5.32** |

This exactly reproduces the reviewer's table and the plotted 9.8 and 5.3. Pooling hides the very comparison the paragraph exists to make: the double knockout sits at the AID-null level (1.83 and 1.00), which is the biological result.

Scope of this bug: **OT1 and OT2 only.** Every other dataset in Figure 4 has a single AID-competent case group, so its separation is unaffected:
`jh4` WT2-4; `n114` C1/D3/D5; `pim1`, `pim1_deep`, `su_lps` WT1-3; `su_ch12` WT1-4.

**Also note (feeds issue 7):** the two files use different denominators. `offtarget_*.json` reports `frequency_bases` as a flat nominal 2,000 reads x window width (1,936,000 for OT1 in every library); `benchmark_real_*.json` reports actual callable bases per library (1.86-1.91M). The manuscript must say which denominator each displayed number uses.

**Also note:** the text's 1,179-fold and 7,139-fold come from the geometric-mean ratio with zero-call controls floored at `1e-6` (`crm_numbers.py:458`). Reproduced: OT1 = 1,181. A fold computed against a floor is an artefact of the floor, not a measurement, and it is why the text and figure differ by two orders of magnitude.

### Resolved (user decision, 16 Oct): run_project is authoritative for biology; no floored folds; KO2 stays excluded

Done in c5ad350:
- Figure 4F gives OT1 and OT2 one row per AID-competent group ("OT1 Wt", "OT1 dKO"). SHMtrace now reads 53 and 1.8 at OT1, 28 and 1.0 at OT2. Panel G, whose counts belong to the dataset rather than the case group, names each dataset once instead of repeating its controls on both rows.
- The Results quote no fold. They give per-library counts, group mean densities with SDs, and the Welch tests on log density and on frequency. OT1 wild type 3.32 against 0.05 per 10^5 bases (p = 6.2 x 10^-4); OT2 1.39 against 0.02 (p = 4.3 x 10^-4). The double knockout is 0.03 and 0.00, separable from neither control group (p = 0.45 and 0.56) and below the AID-competent libraries (p = 0.034 and 2.6 x 10^-5).
- The paragraph states that control KO2 was left out of the decoys at OT1, so that threshold rests on two controls.
- The Figure 4 legend defines both OT groups and says panels A-I come from the caller-comparison run while the text's OT counts come from the full pipeline.

Still open from this issue: Table S2 is still labelled as covering every caller and dataset (C29), and its OT rows come from the caller-comparison run; the table needs a column or note naming the run behind each row.

---

## Status of the remaining major issues

| # | Issue | Status |
|---|---|---|
| 1 | OT1/OT2 analysis provenance | **fixed** (c5ad350); Table S2 labelling still open (C29) |
| 2 | Platform-wide calling claim contradicts Sanger Methods | **fixed** (see below) |
| 3 | Empirical FDR presented as stronger control than established | **wording fixed**; empirical validation of the whole procedure still open |
| 4 | Noisy-control policy ("sets aside" vs "still contributes") | **partly fixed** (c5ad350): wording reconciled, three-control minimum and provenance stated. Still open: side-by-side all-controls vs exclusion results, and validation under heterogeneous nulls |
| 5 | "Reproduced unaided" overstates the Sanger evidence | **fixed** (e54003f) |
| 6 | SWHEL host identity and change of statistical test | **fixed** (e54003f) |
| 7 | Estimands undefined (three different frequency metrics) | **fixed** (74aed94); Figure 3F truth denominator still open |
| 8 | >=50% mask can discard genuine SHM in expanded clones | **fixed**: tested on mutated dominant clones; removed sites now exported |
| 9 | Cross-platform concordance supports similarity, not equivalence | **fixed** (e54003f) |
| 10 | Independent units and benchmark uncertainty | **fixed** (see below) |
| 11 | Biological QC identifies suspicious calls, does not prove origin | open |
| 12 | Reproducibility archive and Methods incomplete | open |

Inconsistency table C01-C31 and the figure-by-figure plan are tracked in the review itself; items resolved here are C01-C07 (cause identified) and C30 (Figure 1 conflation, pending issue 2).

Explicitly **not** contradictions, per the reviewer: the low Mkln1 knockout spectrum correlation (already disclosed, see `CATALOGUE.md`), 40/0.18 vs ~229-fold in Table S3 (rounding), 526 vs 529 simulated bases (excluded bases). These need explanation in the text, not new numbers.

---

## Decisions taken

1. **OT authoritative source: `run_project`** for the biology; `benchmark_real` for the caller comparison only, with Figure 4 labelled so the two roles are separate. (16 Oct)
2. **No folds against floored zeros.** Groups are compared by density with a test. (16 Oct)
3. **KO2 stays excluded** at OT1, and the exclusion is now stated in the Results, with the rule's minimum and provenance in the Methods. (16 Oct)

Note on provenance, checked rather than assumed: the rule was introduced on 13 Oct 2026 for the Ramos clones (a public dataset), applied unchanged to all eleven datasets, and was already fixed when the OT libraries were run on 15 Oct. It is therefore derived from data, not pre-registered, and the Methods now say so.

---

## Major issue 2 — the calling claim covered Sanger, which it does not (fixed)

**Verified.** The Summary, eTOC, graphical abstract and Figure 1 all said, or implied, that every candidate is tested against matched AID-null error. The Sanger Methods say the opposite and are the accurate description: every substitution at base quality >= 20 is counted, and the control background is reported beside the calls — "Nothing is removed by this". The SWHEL re-analysis has no AID-null control at all.

Fixed by making the two routes visible everywhere the claim appears:

- **Figure 1C** is now a fork rather than one chain. A teal route for reads (control error model -> test and FDR) and an orange route for Sanger (count the mismatches -> background reported) rejoin at the shared germline filter and window (step 3) and automatic checks (step 4). Panel A notes it describes the read route; panel B says steps 2 and 3 are shared, while trace trimming replaces step 1 and depth matching cannot apply.
- **Figure 1's title and legend** now name both routes and state that the AID-null background is reported beside Sanger calls and never subtracted.
- **Summary**: "tests each candidate mutation in a read ... A Sanger clone is one molecule read once, so its substitutions are counted and the same controls give the background they carry." 150 words.
- **eTOC and graphical abstract**: same split; the GA's scatter is now labelled "pooled substitution spectrum, nanopore vs Sanger, r = 0.92".
- **Results**: the heading is now "SHMtrace calibrates reads against the AID-null controls", and a new paragraph states that no per-call test is possible for a Sanger clone, what is done instead, and that the Sanger results test trace reading, counting, bookkeeping and group comparisons rather than the control-calibrated test.

Follow-up in Figures S12-S14 and the text (fixed):
- **S12**: cards retitled "From the Sanger clones here" / "From the nanopore reads here"; the legend says the split describes the two workflows evaluated, not the platforms (Sanger clones spanning the junction can be grouped into lineages, UMI consensus gives per-molecule values, read span is assay-dependent). "exactly" removed; "known false-discovery rate" -> "estimated"; trees cross-reference S10 -> S11 (C12); the nanopore fold now says controls without calls count at one event.
- **Singletons (C13, C14)**: SHMtrace individually calls 30% of planted single-read mutations on J_H4 and 83% in the polyclonal simulations, so "only estimated" (S12) and "in aggregate rather than individually" (Discussion, Limitations) were wrong. All now say single-read mutations *not called* are estimated in aggregate.
- **Length gates (C16)**: S13/S14 said reads "longer than the amplicon" are dropped; the presets drop reads over 2,000 bp (J_H4, a 1,230-bp amplicon) and 900 bp (N114, 655 bp). Both figures now say the limit is set per locus (Methods give the values).
- **>= 50% rule (part of issue 8)**: S13, S14 and the S13 legend now call it a heuristic that sets changes aside "as likely germline", and the legend notes it can remove a real mutation shared by a dominant clone. The fuller test of issue 8 is still open.
- **S14 error-rate equation**: showed only the per-position ratio. The code (`null_model.qrate`) uses the pooled rate g as a floor and shrinks each position's own count toward it with m = 500 pseudo-observations: q_x = kappa * max(g, (k_x + m g)/(n_x + m)). The equation and stage text now match.

## Open decisions (user)

None outstanding.


---

## Major issue 8 — the >=50% rule can remove real mutations (fixed)

**Verified, and the cost measured.** The rule sets aside any change carried by at least half the reads covering
it. The existing demonstration (Figure 3E) used an *unmutated* dominant clone, which the rule cannot touch, so it
could not show the failure the reviewer describes.

Three new scenarios put the trunk mutations in the dominant clone itself (`lineage_sim.SCENARIOS`,
`expanded_mutated_40/60/80`), three libraries each, every J_H4 AID-null control planted in turn:

| Dominant clone | Planted mutations reaching the threshold | Recovered | Read-weighted frequency | Clone-weighted |
|---|---:|---:|---|---|
| 40% of reads | 0 | – | within 5% of truth | within 14% |
| 60% of reads | 0-4 per library | none | to 40% below truth | to 26% below |
| 80% of reads | 0-4 per library | none | to 62% below truth | to 34% below |

Nothing that crosses the threshold survives, because the rule runs before any test. The loss is bounded by the
dominant clone's own trunk rather than by library size, and falls hardest on the read-weighted frequency, which
those mutations dominate; the clone-weighted frequency, which counts each clone once, is hurt less.

`lineage_sim.score` now reports `planted_at_mask` and `planted_at_mask_called`, so the cost is measured on every
future run rather than argued.

**Export (the reviewer's second request).** The removed sites were recorded internally but never surfaced.
`project.public` now gives each sample `high_fraction_sites_removed` and `high_fraction_sites` (position,
alternate base, reads, depth, frequency). This also separates the two cases the reviewer wanted distinguished: a
true germline difference sits near 100% (the five J_H4 lab-haplotype sites are at 99.7%), while a dominant
lineage's trunk sits near that clone's read share.

Wording: the rule is called a heuristic rather than a test wherever it appears (Results, Limitations, Figures S13
and S14), and the Figure S13 legend notes it can remove a real mutation shared by a dominant clone.

Still open from this issue: SHMtrace does not yet *recover* such mutations, only report them. A clone-aware
exemption (a change confined to one junction-defined lineage) would be the fix, and is not attempted here.

---

## Remaining numbered inconsistencies (fixed)

- **C23 "six loci".** Checked against git history: when the phrase was written (ed190f8, 13 Oct), Table 1 had exactly
  six SHM target loci (J_H4, N114, Smu, Pim1, CH12 VDJ, Ramos IGH VDJ) from five sources. OT1, OT2, mOrange and SWHEL
  V_H10 were added later, and two laboratories with them, so the count went stale. The same rule now gives ten loci
  from seven laboratories on three platforms (Summary, eTOC, Introduction). Table 1 states the rule and gains the
  SWHEL row it lacked.
- **C25.** Figure S10D marks the total-mutation row (a) and footnotes that the two columns use different tests
  (one-way ANOVA in the paper, Welch on per-animal log frequency here); the legend says they agree in direction, not
  inference.
- **C30.** Figure 1A no longer attributes pseudoreplication to read-level callers; panel A is about error
  calibration only, and group comparison stays in panel D. Both FDR mentions in Figure 1 now say "estimated".
- **C11.** `make_review_package.py` rebuilds the package from the current files: PDFs converted fresh, every image
  copied from crm/figures (the files the supplement embeds), Table S4 rendered from the supplement page itself.

All 31 numbered items are now addressed.

---

## Major issue 7 — the three measures were not defined, and one was misdescribed (fixed)

**Verified, including an error the reviewer did not see.** A new STAR Methods section, "What each mutation measure
counts", gives each measure explicitly with a worked example (wild-type J_H4 library Wt2, 500 reads, 529-bp window):

| Measure | Definition | Wt2 |
|---|---|---|
| Unique-variant density | D / B x 10^5; B is base observations (covered positions x mean depth), not reference positions | 221 / 263,865 = 83.8 x 10^-5 |
| Read-weighted frequency | sum of each read's called sites / (reads x window) x 10^3 | 747 / 264,500 = 2.82 x 10^-3 |
| Clone-weighted frequency | mutations in >= half a clone's reads, over clones x window, x 10^3 | 207 / (24 x 529) = 16.3 x 10^-3 |

The section also states that "unique" means distinct within a library and not of independent origin, and that
unique-variant density falls with depth at a constant mutation rate, so it is compared only at matched depth.

**The error.** Four places said the Sanger frequency divides by "the bases each clone read at base quality >= 20".
It does not: `frequency_bases` is clones x 529, a fixed length per clone, verified against `out/sanger_mkln.json`
(every mouse's bases equal clones x 529 exactly). SHMtrace therefore uses the same *kind* of denominator as the
published studies; what differs is the length (529 against 556 and 577) and which clones are counted (paralog
amplicons and QC failures excluded). The Table S4 footnote, the Results, and Figure S12 now say that, which also
makes the "counted as published" row mean what it says: the same mutations over the published denominator.

Figure 5A is retitled "Within-platform genotype contrasts", since its two densities are not on a common scale.

Still open from this issue: Figure 3F does not yet state its truth denominator at each depth.


---

## Major issue 10 — units, uncertainty and the statistics specification (fixed)

**Unit of analysis.** Table 1 gains a column naming the independent replicate for every dataset: mouse, subclone,
host animal, culture, clone, molecule, experiment, replicate or simulated library. Three footnotes record what the
records do not establish, rather than implying more independence than exists:

- the three N114 subclones per vector came from one transduced population, so they are independent subclones but
  not independent transductions;
- the J_H4 Sanger row is the same six animals as the nanopore row, so the two are a platform comparison within
  animals, not independent experiments;
- PRJNA918596's "four replicates" are the deposit's own numbering, and what they correspond to experimentally is
  not stated in its records.

**The permutation test cannot reach p < 0.05 at n = 3.** Verified in code and in the saved simulations: with three
libraries per group there are C(6,3) = 20 allocations, so the smallest two-sided p is exactly 0.10. In all 300
simulated experiments at n = 3, at every effect size including 4-fold, the smallest p observed was 0.100 and the
share below 0.05 was 0. Figure 3G's zero bars for that test were therefore arithmetic, not a failure to detect,
and the legend and Methods now say so. The Welch test on log rates is named as the primary test, chosen in
advance, with the other two as sensitivity analyses.

**Monte Carlo uncertainty.** Each Figure 3G cell is 300 experiments, so a share carries a 95% Wilson interval no
wider than ±5.6 points; the legend states this. Under the null at n = 3 the Welch test rejected in 3.7%
(95% CI 2.1–6.4) and the negative-binomial test in 5.0% (3.1–8.1), against 91.3% (87.6–94.0) for pooling reads.

**The negative-binomial test is now specified**: mean n·r, variance n·r + φ(n·r)^2 with one φ shared by both
groups and fitted by maximum likelihood under the full model; likelihood-ratio statistic against a common rate;
null from 499 parametric bootstrap samples under that common rate with the fitted φ, because the asymptotic χ²
does not hold at these group sizes; a group with no calls keeps rate zero rather than being floored.

**All three measures are now reported where relevant**: Table S1 gives unique-variant density, read-weighted
frequency and clone-weighted frequency for every library, with a footnote saying they are not interchangeable.

Still open from this issue: Figure 7F does not yet show experiment-level points with their pairing.
