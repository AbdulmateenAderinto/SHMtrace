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
| 3 | Empirical FDR presented as stronger control than established | open |
| 4 | Noisy-control policy ("sets aside" vs "still contributes") | **partly fixed** (c5ad350): wording reconciled, three-control minimum and provenance stated. Still open: side-by-side all-controls vs exclusion results, and validation under heterogeneous nulls |
| 5 | "Reproduced unaided" overstates the Sanger evidence | open |
| 6 | SWHEL host identity and change of statistical test | open |
| 7 | Estimands undefined (three different frequency metrics) | open — confirmed by the OT denominator finding above |
| 8 | >=50% mask can discard genuine SHM in expanded clones | open |
| 9 | Cross-platform concordance supports similarity, not equivalence | open |
| 10 | Independent units and benchmark uncertainty | open |
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

Remaining from this issue: Figures S12-S13 still carry the exclusive platform claims the reviewer flagged (C15 is fixed; S12's "only estimated" singleton claim, the S10 -> S11 cross-reference and the length-gate contradiction C16 are not).

## Open decisions (user)

None outstanding.
