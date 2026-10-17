# Response to your 17 October review

Thank you for the review. It found real defects, and several of them were in work done in the previous
revision pass. This package shows what was done with each point, so you can check it rather than take it
on trust. Every number below is recomputed by the packaging script and recorded in `evidence/`.

Repository state: `{head}` (working tree clean). Your patch was made against `83a908c`.

## What happened to your patch

**Applied unmodified** as commit `3317fd9`, so its provenance stays separate from what followed. Before
applying, each confirmed defect was checked against the code. Two of them also got an independent test
beyond your regression tests:

- **Shared alleles.** Your harness fix was applied to a scratch copy and two shared-artifact seeds were run
  before touching the repository; both gave 0 calls, matching your result.
- **Benchmark window.** Your note that only J_H4 and N114 were affected was checked by computing the old and
  new windows for all seven benchmark datasets. For the other five the end trim sets the window and the two
  old shifts cancel exactly, as you said.

Everything else in the patch was read before applying: the `build.js` profile change, the tests, the
figure and table builders, `validation_summary.py`. Nothing does more than its stated purpose.

## Outputs regenerated here, not copied from your package

`crm/out/` is not under version control, so adopting your fixes meant either importing your output files or
regenerating them. They were regenerated with the patched code:

| Output | Result against yours |
|---|---|
| FDR validation, `fdr_validation` + `fdr_null_reps` + `fdr_estimator` | {fdr_same}/{fdr_total} runs identical in every result field |
| J_H4 and N114 benchmark, `benchmark_real.py` | {bench_same}/{bench_total} rows identical in calls and base counts |
| Table S2 workbook | {s2_summary} |

{fdr_order_note}

The pre-review outputs are kept in `crm/out/pre_review_2026-10-17/` for comparison.

## Two places this differs from your patch

1. **Figure 5A title.** You shortened "Within-platform genotype contrasts" to "Genotype contrasts" because it
   collided with panel B. The collision was real. But "within-platform" was the wording the first external
   review asked for (its issue 5), because panel A compares genotypes within each platform and should not
   be read as a cross-platform comparison of densities. The title is now **"Within-platform contrasts"**,
   which clears panel B with room to spare. See `figures/Figure5.png`.
2. **Summary.** Your rewrite is kept; it is more accurate on the two routes and on what the re-analyses
   recovered. It had dropped two things, both restored: the motivating clause ("because single reads err at
   rates close to the mutation rate") and "from the same run", which is the design point the control
   calibration rests on. The Summary is {summary_words} words (limit 150).

The exact source changes beyond your patch are in `patches/changes_beyond_your_patch.patch`; apart from the
two files above they are the tracker and handover corrections below.

## REVIEW_RESPONSE.md and HANDOVER.md

You treated these as claims, and several were contradicted by the implementation. Your patch did not touch
them, so after applying it they still stated false things as fact. Each is now corrected in place and marked,
rather than deleted, so the record shows what was believed and when:

- **Repair code.** "Repair-pathway indices appear nowhere in the code" is corrected. `repair.py` exists; the
  real reason the paper does not report it is that the project workflow does not use it.
- **Q20.** Q ≥ 20 is no longer described as the read-route test cutoff.
- **Benchmark harness.** It is no longer described as selecting no threshold.
- **Validation results.** These are replaced with the regenerated ones.

C29 is closed by your Table S2 About sheet. The handover's closing advice previously cited the false
repair-code claim as an example of a correct check; it now records what went wrong. Both documents are
included as `05_REVIEW_RESPONSE.md` and `06_HANDOVER.md`.

## Where things stand

- Main text **{main_words}** words (cap 7,000); Summary {summary_words}; eTOC {etoc_words}; highlights {highlights} characters.
- **{tests_passed} tests pass**, including your seven; `ruff check src` is clean.
- All supplementary figures and tables are cited; the figure files the documents embed all exist.
- **{blockers} submission blockers** remain in the shared registry. As you noted, not all of them appear in this
  manuscript: some are only in the companion protocol. The build still prints the registry total rather
  than the manuscript-specific count; that has not been changed.

## What we would like you to check

1. **The two deviations above.** Do you agree, or would you handle the Figure 5 collision differently?
2. **The AI-use declaration**, which now names Anthropic Claude and OpenAI Codex. It is accurate now that your
   changes are adopted, but it is the authors' statement, and the journal's current requirements should be
   confirmed before submission.
3. **The rendered PDFs.** Your layout fixes (US Letter page height, S4 and S11 split across two captioned
   pages, repeated table headers) were applied and rebuilt here. Please confirm they render as you intended.
4. **Your view on an independent validation run.** It is the main open suggestion from your review: the
   stress tests resample three existing AID-null libraries from one run, so they test the procedure
   conditionally, not its transfer to a new preparation. Is it needed before submission, or is it better
   stated as a limitation?

## Not done, and said so

- **Mean-depth rounding.** The precision issue you deferred (mean depth rounded before a denominator is
  rebuilt from it) is still deferred, for the reason you gave: fixing it means refreshing every affected
  metric consistently.
- **Optional suggestions.** Provenance columns in exported tables, an opt-in control and downsampling
  sensitivity report, and a related-junction confusability test are recorded as open, not started.
- **Duplicate bullet.** Your review's "Figure 2 and the read-weighted-frequency definition" bullet appears
  twice (section 10). It does not affect anything; it is only noted so the copy you keep is clean.

## Contents

| Path | What it is |
|---|---|
| `00_README.md` | this note |
| `01_Manuscript.pdf`, `.docx`, `.txt` | the article |
| `02_Supplement_Document_S1.pdf`, `.docx`, `.txt` | Document S1 |
| `03_Highlights_eTOC_Motivation.pdf`, `.docx` | front matter |
| `04_Table_S2_benchmark.xlsx` | Table S2 |
| `05_REVIEW_RESPONSE.md`, `06_HANDOVER.md` | the corrected tracker and handover |
| `figures/` | every standalone figure image; `.scifig` sources for Figure 1 and the graphical abstract |
| `patches/` | source diffs since your snapshot, and beyond your patch |
| `data/` | the regenerated validation and benchmark outputs |
| `evidence/` | the comparison against your outputs, the build log, the test result |
