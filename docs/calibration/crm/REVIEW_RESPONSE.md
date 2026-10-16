# Response to the external review of 16 October 2026

Review: `~/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/nanoshm/output/reviews/SHMtrace_manuscript_review_2026-10-16.md`
Recommendation: major revision before submission.

Status key: **open** / **investigating** / **fixed** / **won't fix (with reason)**.

---

## Major issue 1 — OT1/OT2: two different analyses (investigating)

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

### Decision required before fixing
Which source is authoritative for the OT paragraph? See "Open decisions" below.

---

## Status of the remaining major issues

| # | Issue | Status |
|---|---|---|
| 1 | OT1/OT2 analysis provenance | investigating — cause identified, decision needed |
| 2 | Platform-wide calling claim contradicts Sanger Methods | open |
| 3 | Empirical FDR presented as stronger control than established | open |
| 4 | Noisy-control policy ("sets aside" vs "still contributes") | open — linked to issue 1 |
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

## Open decisions (user)

1. **OT authoritative source.** Recommendation: use `run_project` (`offtarget_*.json`) for the biology, because it is the method the manuscript describes; keep `benchmark_real` for the caller comparison only, and relabel Figure 4's panels so the two roles are visibly separate. The alternative is to re-run the OT libraries through `benchmark_real` only, which would mean the headline OT numbers no longer come from the full method.
2. **OT fold convention.** Recommendation: stop quoting a fold against floored zeros. Report wild type vs AID-null and double knockout vs AID-null as separate three-group comparisons with counts and the test, and state the one-event convention wherever a fold is unavoidable.
3. **Noisy-control exclusion (issue 4).** Whether KO2 of OT1 should be excluded at all, and whether the rule was fixed before the public data were analysed. This determines both the OT numbers and the honesty of the FDR claim.
