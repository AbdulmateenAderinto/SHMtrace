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
| 3 | Empirical FDR presented as stronger control than established | **fixed**: wording, the Benjamini-Hochberg reason corrected, and the whole procedure tested on fresh simulations (Figure S15) |
| 4 | Noisy-control policy ("sets aside" vs "still contributes") | **fixed**: wording reconciled, minimum and provenance stated, and all-controls vs exclusion run side by side on identical libraries under heterogeneous nulls (Figure S15) |
| 5 | "Reproduced unaided" overstates the Sanger evidence | **fixed** (e54003f) |
| 6 | SWHEL host identity and change of statistical test | **fixed** (e54003f) |
| 7 | Estimands undefined (three different frequency metrics) | **fixed** (74aed94); Figure 3F truth denominator still open |
| 8 | >=50% mask can discard genuine SHM in expanded clones | **fixed**: tested on mutated dominant clones; removed sites now exported |
| 9 | Cross-platform concordance supports similarity, not equivalence | **fixed** (e54003f) |
| 10 | Independent units and benchmark uncertainty | **fixed** (see below) |
| 11 | Biological QC identifies suspicious calls, does not prove origin | open |
| 12 | Reproducibility archive and Methods incomplete | **fixed** apart from 11 laboratory details only the authors can supply |

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

Closed in the figure-by-figure pass: Figure 3F's legend now states the truth denominator (every planted
mutation present in the analysed reads, over the estimator's own denominator of analysed positions x mean
depth) and why unique density falls with depth.


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

Closed in the figure-by-figure pass: Figure 7F now shows every experiment, SHMtrace's share and the
every-mismatch share joined across arms within each experiment, with each dead-AID library as a cross.


---

## Major issue 12 — the archive and the missing Methods (fixed, except what only the authors have)

**The archive was further along than the review package could show.** The data deposit (Ver19) already carries a
file manifest with SHA-256 checksums, sample manifests, FASTQ and GenBank validation, read accounting, and a
folder per analysis. The reviewer saw only the manuscript, so this was invisible.

What was genuinely missing, and is now added:

- **`ARCHIVE_ENVIRONMENT.txt`**: interpreter, package versions, external tools with their reported versions, and
  the commit the manuscript was built from.
- **`ARCHIVE.md`**: one row per analysis giving the command, its inputs and its outputs, for all 18 analyses
  behind the figures, tables and quoted numbers. `make_archive_manifest.py` regenerates both and checks each
  output exists, so a stale row reads "missing" instead of passing silently. All 18 are present.
- **Resource availability** now points at both, and states which alignment path each analysis uses (mappy
  bindings for SHMtrace, the standalone minimap2 executable for the LoFreq comparison), settling minor 15.
- **The Methods gaps the reviewer listed**: the error-rate floor and shrinkage (pooled rate as a floor, each
  position's own count shrunk toward it with m = 500), the probability floor (half a pseudo-count for a class
  the controls never showed, so no p value is exactly zero), the indel caller (class rate as floor, site rate
  shrunk toward it, margin fitted as for substitutions, Benjamini-Hochberg, quality >= 30, three reads minimum),
  UMI consensus (3 edits per tag end, directional merging, groups of three or more) and the single-read
  quality-eligibility correction. The bootstrap was specified under issue 10.
- **Every build now lists the remaining submission blockers** by name instead of leaving a placeholder to be
  noticed by eye.

**Still outstanding, and only the authors can supply them (11):** GL7 catalogue number and RRID; proteinase K
catalogue number; Monarch kit catalogue number; the Aicda-/- MGI RRID; FACS buffer composition; lysis buffer,
volume and incubation; lysate volume per PCR; the SRA BioProject; the data deposit DOI; the SHMtrace Zenodo DOI;
and the companion article citation once accepted. These are not guessed: catalogue numbers and RRIDs would be
fabricated if invented, so they are left as named placeholders that the build reports.

The AI-use declaration now says the assistance covered drafting the analysis and figure code, language editing
and figure layout, which is what it covered.

---

## The figure-by-figure plan

Worked through in the review's own order. Figures are built with the project `.venv` (matplotlib 3.11); the
system interpreter's 3.10 cannot render Figure 2's underlined sequence context.

### Figure 2
- **Q >= 20 labelled.** Panels A and B counted alt reads at any base quality while the test counts only bases at
  quality >= 20. Each count now gives both (165: 20 reads, 19 at Q >= 20; KO2's and KO3's single reads there are
  below Q20 and never entered the test). The axis says "any quality" and the Q20 line in C is labelled.
- **Panel E retitled** "Allele fractions of SHMtrace calls".
- **The shaded band.** It ran 10-60% and read "somatic callers work here", which claimed a general operating
  range. It now covers only the 5-20% range Figure 3H actually tested, with the measured recall printed beside
  it (DeepSomatic 50/77/93%, ClairS 0/27/97%), and the legend calls it a measured range for those settings on
  J_H4, not a limit of either caller. The main text's "allele fractions at which these callers work" is reworded
  to match at equal length.

### Figure 3
- **Uncertainty.** A, B and F carry the range over the three held-out controls; D the range over the three
  simulated libraries. A and B print the pooled FDR for SHMtrace and the Geneious rule.
- **Panel E.** Each library is drawn as a stem from its own truth to its own estimate, because the truth differs
  between libraries (one expanded-clone library's truth is 12.1 against a mean of 15.5), which a single mean bar
  concealed.
- **Panel H clipping.** Widened; title, label and axis no longer clip.
- **Denominators.** The legend gives 90 planted mutations per point in A and B, defines the low-load scenario and
  the support bins in D, states F's truth denominator (every planted mutation in the analysed reads, over the
  estimator's own denominator) and why unique density falls with depth for N114.

### Figure 4
Panel G already names each dataset once and the OT groups were split earlier (issue 2). The legend now cites
Figures S4-S6 and Tables S1 and S2.

### Figure 5
Panel A was retitled "Within-platform genotype contrasts" earlier (issue 5). The legend now cites Figures
S7-S10 and Tables S4 and S5.

### Figure 6
- **Nested node counts.** Panel D is titled a read-level reconstruction, states that 48 of the clone's 69 reads
  fit one path and 21 do not, and says node counts are nested, so they do not sum.
- **Reads that fit no tree.** Panel C now spans all 500 reads per mouse, so the reads whose junction could not be
  read (21-45%) are visible instead of being divided away.

### Figure 7
- **Panel subtitles.** Each panel carries its validation question in italics beneath the title.
- **Experiment-level points in F.** It had plotted three arm means and shaded the dead-AID control as the range
  of those means, hiding a 43-66% spread across its nine libraries. It now shows every experiment: SHMtrace's
  share and the every-mismatch share joined across arms within each experiment, and each dead-AID library as a
  cross. The legend says SHMtrace's shares rest on 4-40 called C:G mutations per library.
- **Panel E.** Molecule counts per bin, and the aggregate reduction as a point with a 95% bootstrap interval
  over molecules: 0.17% (0.10-0.28), 57-fold (35-99). The interval is now in the text too.
- **Every bar, interval and symbol defined**, with the sidedness of both tests in A and the one-event censoring
  the Welch test uses (minor 8).

### Figure S1
It was shared with the JIM submission, whose figure must not change, so the article now builds its own copy
(`scicanvas/figS1_workflow.js`). In that copy the Sanger route says a clone is one molecule, counted and not
tested, over clones x 529 bp; the >= 50% rule is named a heuristic; "lineage trees" reads "read-level trees".
S1 now agrees with Figure 1C instead of contradicting it.

### Figure S2
95% Clopper-Pearson bands on every curve, bins with no mismatch drawn at their upper bound instead of dropped,
and a second panel giving the control bases per bin (above Q40 they thin out, which is where the curves wander).
Run matching is stated from the data: all six J_H4 libraries from one order (DHJGKB), all six N114 from another
(75K3B3); PRJNA1481264's records do not say which libraries shared a flow cell, so nothing stronger is claimed.

### Figure S3
- **Target power, model and the required n.** `power_s3.py` re-simulates from the pilot out to 60 animals:
  unique density reaches 80% power at about 24 per group. The legend states alpha, sidedness, the test, the
  count model and the pilot.
- **CV uncertainty.** The CV comes from 3 mice, so its 95% interval has 2 degrees of freedom and an effectively
  unbounded upper end. The band shows the lower bound (0.44) and twice the estimate (1.9) as a stated
  pessimistic case, and the legend says the CV, not the simulation, limits this calculation.
- The Discussion said "more than ten animals per group"; it now says about twenty.

### Figure S4
The three components are JIM figures and cannot be renumbered, so the legend identifies each part by position
and says what its own letters mean. **The reported "density axis labelled frequency" could not be reproduced**:
the only quantitative axis in the three components reads "Mutation density (x10^-5 mut/bp)" and the matrices
carry no density axis. The legend now states the measure explicitly, so the point is covered either way.

### Figure S5
- **Repeated downsampling.** `depth_draws.py` repeats the whole analysis five times at each depth with a
  different draw, controls redrawn and recalibrated each time. The concern is borne out: at 30 reads J_H4 Wt2
  gives 7 to 60 calls across five draws, against 212 to 244 at 500.
- **Zero handling.** An AID-null library with no calls has a density of zero, which has no logarithm. Those
  draws are drawn on a labelled "0" row rather than floored silently, and they are common: J_H4 KO1 and KO3
  called nothing at any depth in any draw.
- **93%/74%.** Already phrased as "of the calls each made at its largest depth" in the text, not as saturation;
  the figure and legend now say the same, explicitly against the largest depth tested. Both numbers now come
  from the repeated draws rather than one: J_H4 is unchanged at 93%, N114 moves from 74% to 70%.
- The two stacked components are replaced by one continuously lettered figure (A-E).

### Figure S6
The title claimed the background "is position-specific" while the panel's own number said otherwise: between
controls Spearman r is 0.12-0.20 at J_H4 and 0.29-0.45 at Smu, so one control explains 1-4% and 9-20% of
another's position-level variation. The panel now says that in full, and the legend draws the right conclusion:
recurrent artifacts are a tendency of the chemistry and sequence context, not a fixed per-position fingerprint,
which is why the model pools over substitution type, strand and base quality and shrinks each position toward
the pooled rate. Nothing in the method changes; the figure now argues for what the method does. Each profile
panel also gains the control-only inset the review asked for.

### Figure S7
The legend now carries the measured path from traces to clones: 184 of the 447 chromatograms are filed to these
six mice, 2 are excluded in the sample sheet with a stated reason, and of the 182 analysed 48 yielded no clone
(22 contaminant or paralogue, 12 with no stretch of >= 20 bases at Q >= 20, 11 unmapped, 3 too short), leaving
134 clones, 103 wild-type and 31 Aicda-/-. The counts come from the pipeline's own per-trace fates. The legend
also says these drops are properties of the trace, applied before any mutation is counted, and that the
bookkeeping checks never remove a clone.

### Figures S8-S10
**Table S5** asks all three studies the same five questions: same clones, same endpoint, same denominator, same
test, and where the published values came from. Several answers are no. The Maea and Mkln1 denominators differ
from the papers'; for SW_HEL the denominator and the test both differ, so that claim agrees in direction only
and the table says so. The table also separates original published values (Maea Source Data) from values
recomputed from printed tables, and says where Figures S8F and S9D use the one-tailed test the original studies
used while Table S4 reports two-tailed throughout.

### Figure S11
The ancestor annotation sat below the first branch and fell outside the axes in a single-row panel, so Wt4's was
clipped; each node's read count was drawn inside its own circle, which grows with the read share, so the largest
nodes were the least readable; mutation labels were 5.6 pt, below the 6 pt floor for print. The note now sits
above the branch with room reserved, counts sit above their circles, labels are 7.2-7.5 pt, the panel is larger,
and the figure is also written as SVG.

### Figures S12-S14 and the graphical abstract
Done earlier under issues 2, 6 and 11.

### Cross-references
Figure S6 and Tables S1, S3 and S5 were never cited. All 14 supplementary figures and all 5 tables are now
cited, from the legends, which are not word-capped; the main text holds at 6,999 words.


---

## Major issues 3 and 4 — the validation halves (fixed, Figure S15)

Both issues ask for the same thing from two directions, so one experiment answers both. The margin, the
threshold and the noisy-control rule are all chosen from the controls, so they have to be tested together,
and on data none of them has seen.

**Design** (`crm/fdr_validation.py`). A case library is built from one *Aicda*-/- J_H4 library's reads, with
a known set of planted mutations or none at all. The controls are two further *Aicda*-/- libraries, each
split in two, so no control shares a molecule with the case while two controls come from each sequencing
library, as two controls from one run do. Every library is then analysed by `run_project` itself -- the same
entry point as every result in the paper -- so the threshold and the rule are selected afresh each time.
This is the distinction the review drew: `benchmark_real.py` fits a leave-one-out null for a caller
comparison, while this exercises the complete procedure.

Four conditions, from the review's own list: the controls as delivered; one control carrying extra errors at
six recurrent positions in 6% of its reads; the same artifact in the case and in no control; and the same
artifact in the case and in every control. Control number (2 or 3), depth (250, 500, 1,000) and mutation
burden (none, or about 400 planted) are varied, and **every condition runs with the rule on and off on
identical libraries** -- the side-by-side issue 4 asks for. 110 runs in the main grid, 80 more for the
all-null cells, 30 for the estimator arm.

### What it found

**The measures the review asked for, which the paper did not report.** All-null libraries, 500 reads,
3 controls:

| controls | rule | libraries with any call | control set aside | mean false calls |
|---|---|---|---|---|
| as delivered | on | 1/25 (95% upper 20%) | 3/25 | 0.04 |
| as delivered | off | 0/25 (95% upper 14%) | 0/25 | 0.00 |
| one made noisy | on | 0/25 (95% upper 14%) | 25/25 | 0.00 |
| one made noisy | off | 0/25 (95% upper 14%) | 0/25 | 0.00 |
| artifact in case only | either | 5/5 | 0/5 | 6.00 |
| artifact shared | either | 4/5 | 0/5 | 1.20 |

The first two columns make the review's point exactly: a mean of 1.2 false calls looks small, but 4 of 5 of
those libraries discovered something, and in a library with nothing to find one call is a false-discovery
proportion of 1. Both are now reported.

**Three findings.**

1. **An artifact the controls do not carry is called every time**, under either setting, in every run. The
   controls cannot measure what they never show, and no threshold repairs it. Shared with the controls, the
   same artifact is largely removed (6.00 to 1.20 false calls). This is the limitation the Discussion already
   stated; it is now quantified rather than asserted.
2. **With real mutations present the per-library FDP stays near target** -- 3.8% with the controls as
   delivered, 2.9% with one noisy, 2.9% shared -- and rises above 5% only when the case carries an artifact
   the controls lack (5.9%). Six false calls among no true calls is an FDP of 1; the same six among four
   hundred is under 2%, which is why the per-library proportion and the pooled ratio are both given.
3. **The rule's fold comparison degenerates when the peers are clean.** Five times a median of zero is zero,
   so a control with at least five held-out candidates is set aside whenever its peers have none, however
   mild its excess. It fired on undoctored controls in 3 of 25 all-null runs for this reason. The Methods now
   say so, and the rule names every library it sets aside.

**The rule, side by side.** Where one control carried the artifact it was set aside in all 25 runs, and no
all-null library called anything with the rule on or off: in the condition the rule was written for it
neither gained nor cost measurably at this depth. With planted mutations it recovered about 3 true calls of
400 (393.0 vs 389.6) for +0.1 percentage points of FDP. The one all-null false discovery in 50 runs came
with the rule on, in a run where it had set aside a control -- the trade the review anticipated, visible but
not separable from chance at this number of runs. With two controls it never fires, as the minimum requires.

**Decoy-tail sensitivity** (`crm/fdr_estimator.py`). SHMtrace averages the held-out controls at each
threshold; requiring the threshold to satisfy the noisiest control is the conservative alternative already
in the code. Measured rather than argued: the conservative choice lowers mean FDP from 3.8% to 2.5% and
costs about 14 of 403 true calls (~3%). Where the case carries an artifact the controls lack it helps little
(5.9% to 5.1%), because no summary of the decoys can represent an artifact the decoys do not contain. The
choice is a modest trade, not what decides the failure mode.

**Benjamini-Hochberg.** The manuscript said the p values are "discrete and conservative, so rather than rely
on their distribution, as Benjamini-Hochberg does". That is not the objection, as the review says: BH
controls the FDR for valid discrete p values, conservatively. The reason is candidate selection -- a
position is tested only where a read already differs, so the tested set is chosen by the same data that
supply the evidence, and a candidate's null p value is conditioned on that selection rather than
superuniform. The control-based threshold is not exposed to this, because the decoys are selected and tested
by the identical rule. Corrected in both the Results and the Methods.

**Still not claimed.** These are simulations on one locus, built from three *Aicda*-/- libraries of one
sequencing run; the artifact conditions rest on 5 runs each, with correspondingly wide intervals. An
independent sequencing run would test the thresholds further, and the Limitations say so.


---

## Minor issues 1-15

| # | Issue | Status |
|---|---|---|
| 1 | Correspondence address truncated | **fixed** earlier: the `.ca` address is escaped so the asterisk no longer eats it |
| 2 | Incomplete sentence in the Discussion | **fixed** earlier, with the reviewer's own replacement |
| 3 | AID-null controls vs biological pseudoreplication | **already separate**: the Introduction says the controls address the first problem and animal-level design the second |
| 4 | Unqualified claim that SHM is measured by cloning | **fixed**: scoped to genomic reporter loci, with repertoire and UMI-consensus sequencing acknowledged |
| 5 | "the weakest link ... is bookkeeping" | **fixed** earlier |
| 6 | "uracils can only be copied" | **fixed** earlier: now "are expected to be copied as C→T or G→A" |
| 7 | Absence of effect inferred from non-significance | **fixed**: effects and intervals now in Table S3 |
| 8 | Figures 3/7 legends: define every mark and test | **fixed**: 7D bars defined, 7F band replaced by per-experiment points, all four group tests in 3G named two-sided |
| 9 | "clone" ambiguous between bacterial and B cell | **fixed**: both defined in STAR Methods |
| 10 | Cell lines: authentication, mycoplasma, transduction vs subclone | **partly fixed**: transduction and subcloning already described (Table 1 footnote a); authentication and mycoplasma are now named blockers |
| 11 | Sex per mouse | **blocker**: the Methods give age, background and "both sexes"; per-mouse sex is now a named blocker |
| 12 | Benchmark framing | **fixed**: germline callers already framed as expected negative-range comparators; the somatic callers now carry the same chemistry caveat as LoFreq |
| 13 | Runtime: hardware, scope, human effort | **fixed**: hardware recorded from the machine, scope stated, curation already excluded |
| 14 | Sensitive mode and repair indices | **fixed**: sensitive mode now has a table row showing its cost; "repair-pathway indices" was an over-claim and is corrected |
| 15 | minimap2 vs mappy versions | **fixed** under issue 12 |

### The ones worth reading

**4.** The Motivation claimed somatic hypermutation "is still measured mostly by cloning PCR products and
sequencing colonies one at a time" as a statement about the field. It is now a statement about genomic
reporter loci such as the J_H4 intron, and says that repertoire sequencing and UMI consensus scale far
better but are aimed at expressed, antigen-selected V regions rather than a short non-coding amplicon from
genomic DNA.

**7.** Two claims leaned on significance. The *Ldha* result said "no detectable difference"; Table S3 now
gives 1.23-fold [95% CI 0.46-3.30], an interval spanning a two-fold decrease to a three-fold increase, so
the comparison establishes much less than "no difference". The repair-context sentence said the
experiment-level tests "did not reach significance despite these differences", and that was simply wrong:
the transition share at C:G falls from 95% under UNG inhibition to 4% under UNG overexpression, and paired
within experiment that swing is 92 points (95% CI 74-110, p = 0.002). The manuscript was understating its
own result.

**9.** The field uses "clone" for two different things and so did the paper. A *Sanger clone* is one
bacterial colony, a unit of sequencing that says nothing about which B cell the molecule came from; a
*B cell clone* is a lineage defined by a shared V(D)J junction, which is what the clone-weighted frequency
and the mutation trees use. The Methods now say this, and note that the Sanger window (positions 101-629)
does not reach the junction at position 53, so lineages cannot be identified in those clones at all --
checked against the preset and the Sanger output, not assumed.

**14.** "Repair-pathway indices" appears nowhere in the code. What SHMtrace reports is the share of calls
at G:C, the share that are transitions, and the shares in WRC/GYW and WA/TW motifs: spectrum and sequence
context, not a mechanistic decomposition. The Results now say "motif shares", and the Methods state plainly
that the products of uracil excision, mismatch repair and replication over uracil overlap, so a spectrum
alone does not separate them; where this paper discusses repair the pathway is set by genotype or vector,
not inferred. The sensitive mode's claims now have a table row that shows the cost as well as the gain:
on J_H4, single-read recall rises from 30% to 72%, and false calls per library from 2.4 to 11.9.

**10 and 11** need the authors. Three new named blockers: per-mouse sex, how the N114 hybridoma was
authenticated, and its mycoplasma status -- each phrased so that "not done" is an acceptable answer, since
the review asks for an explicit statement either way. The build now lists 14 blockers.
