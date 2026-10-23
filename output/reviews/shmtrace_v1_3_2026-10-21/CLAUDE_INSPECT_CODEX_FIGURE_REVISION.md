# Claude handoff: inspect the Codex SHMtrace figure revision

Date: 22 October 2026

## Requested task

Independently inspect the figure revisions Codex made. Assess scientific fidelity, clarity, visual quality at publication width, reproducibility, and suitability for the SHMtrace Cell Reports Methods manuscript. Do not assume the revisions are improvements merely because checks pass. Identify regressions, omissions, misleading presentations, and recommendations that were only partially implemented.

Start with an inspection report. Do not overwrite original figures or the canonical manuscript. If correcting issues, work in a separately named copy and preserve the Codex version for comparison.

## Locations

Repository:
`/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2`

Review directory:
`/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21`

Revision directory:
`/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21/codex_figures_v2`

Within the revision directory:

| Resource | Purpose |
|---|---|
| `SHMtrace_Codex_Figures_v2.pdf` | 44-page review pack: cover plus 43 figure/companion pages, placed at 174-mm figure width |
| `index.html` | Local before/after gallery with individual export links |
| `README.md` | Scope, boundaries, dependencies and rebuild instructions |
| `FIGURE_INTEGRATION.md` | Revised captions, provisional destinations and manuscript integration notes |
| `figure_manifest.json` | Ordered list of included assets and brief change descriptions |
| `revised_legends.json` | Machine-readable revised legends |
| `figures/` | Editable vector SVG, vector PDF, PNG and TIFF exports |
| `sources/revised/` | Active copied builders and layout code |
| `sources/originals/` | Snapshots of pre-revision source files |
| `sources/*_data.json` | Frozen numerical data extracted from existing SciCanvas builders |
| `sources/*_comparison_table.csv` | Comparison tables separated from S14–S16 |
| `inputs/` | Copies of the specific saved analysis resources used by the revision |
| `review/before_*.png` | Baseline images from the Codex 1 manuscript |
| `manuscript/original_operations.json` | Baseline manuscript construction operations; not an updated manuscript |
| `manuscript/original_legends.md` | Original full legends, retained for experimental detail and comparison |
| `validation.json`, `CHECKS.md`, `print_check.json` | Recorded verification results |
| `original_hashes.json`, `input_hashes.json` | Hash manifests for protected originals and recorded original data inputs |
| `input_snapshot_map.json`, `input_snapshot_hashes.json` | Mapping and hashes for the copied input resources |
| `scratch/` | Intermediate summaries and build logs; not authoritative source data |

Portable archive, beside the revision directory:
`SHMtrace_Codex_Figures_v2_sources_and_exports.zip`

The archive was created before this handoff and therefore does not contain this handoff file. It includes sources, copied inputs, gallery baseline images, PDF/SVG/PNG exports and documentation. TIFF exports remain in the local `figures/` folder and are omitted from the archive.

Earlier critique and requested direction:
`CLAUDE_FIGURE_REVISION_HANDOFF.md` in the review directory.

Baseline manuscript:
`codex_version/SHMtrace_CellReportsMethods_complete_v1.3_Codex_1.docx` in the review directory.

Canonical project context:
- `docs/calibration/crm/SESSION_HANDOVER_2026-10-21.md`
- `docs/calibration/crm/HANDOVER.md`
- `docs/calibration/crm/figs/`
- `docs/calibration/crm/scicanvas/`
- `docs/calibration/crm/manuscript/`

Some historical handover statements disagree about numbering and word counts. Establish the current state from the actual baseline manuscript, operations and figure manifest rather than combining incompatible historical notes.

## What Codex changed

### Main figures

1. **Figure 1:** Replaced the rounded-card workflow with an open schematic emphasizing matched inputs, control-informed calibration and outputs. Added explicit separation of read-library variant calling, Sanger clone counting and PacBio switch-junction analysis. The pileups are labelled schematic.
2. **Figure 2:** Enlarged the paired candidate examples, moved the quality comparison into a wide middle panel, removed the somatic-caller range shading/commentary from the allele-fraction histogram, and retained the empirical threshold illustration.
3. **Figure 3:** Retained calling-performance panels A–D. Moved downstream simulations to provisional S25. Preserved actual zero false calls on a symmetric-log axis, and labelled the pooled spike-in error summary as observed FDP.
4. **Figure 4:** Rearranged narrow panels into wider small multiples, shortened method labels with a key, added explicit zero-call strips and heatmap colour scales, and moved read accounting to provisional S26. The original statistical tests and one-event convention for ratios/tests remain unchanged.
5. **Figure 5:** Replaced the dominant bars in C with paired mouse observations, deemphasized mean bars in E, labelled motif availability, and preserved the positional profile and separately marked four-base terminal interval. Absolute platform densities are qualified as different measurement units.
6. **Figure 6:** Lightened individual positional profiles, separated clone composition from read/clone weighting, and clarified reconstruction, nested read counts and the interpretation of reads that do not fit a path.
7. **Figure 7:** Emphasized individual indel observations and split E into calibrated active-AID calls and original Q > 20 mismatch summaries. Marker shape identifies experiment, and trajectories connect the same experiment across repair conditions.

### Supplementary figures

- **S1:** Experimental design separated from software processing.
- **S2:** More readable quality/support plots and explicit uncertainty definitions.
- **S3:** Three representative outputs replace six small instructional cards.
- **S4:** Actual zero FDP retained on a symmetric-log scale; the scale and default target are explained.
- **S5–S7:** Clearer fallback/exclusion conditions, axes, zero-call explanations and planted-event count descriptions. These are smaller presentation changes, not complete redesigns.
- **S8:** Estimation/power and UMI results split across separate pages under the same figure number and original A–D panel letters.
- **S9:** Count matrices replace count-plus-percentage matrices. Library matrices share a count scale within each plate. Zero-event spectra are marked undefined; tiny totals are qualified. The pooled N114 control group was checked against its actual `pMX-Empty` label.
- **S10:** Clarified that recovery uses each library's maximum observed discovery across all depths/draws, not only its 1,000-read result.
- **S11:** Separate control tracks replace overlays and insets, with a common control scale within each locus. Individual-read maps remain.
- **S12:** Absolute frequencies lead; ratios are explicitly secondary. The one-event zero convention is disclosed.
- **S13:** Per-mouse matrices simplified to counts. Existing saved coverage curves were retained as vectors, with denominator clarification and print-width stroke normalization; raw traces were not reprocessed.
- **S14–S16:** Larger reanalysis plots, separately readable numerical comparison tables, and editable CSVs. Numerical plot data come from the existing injected SciCanvas data, including their existing precision.
- **S17:** Nested-window comparison explicitly qualified. Three selected trees now have dedicated pages: WT2 clone 0, WT3 clone 0 and WT3 clone 14. These are not one example per mouse. Fit annotations are scoped to mutation-bearing reads.
- **S18–S19:** Local interpretation/layout clarifications; much of the existing plot structure remains.
- **S20:** Each sample's junction composition is displayed separately instead of a group-average stack.
- **S21:** Dense platform illustrations replaced by a compact measurement comparison table.
- **S22:** Redundant workflow replaced by a comparison of distinct analysis routes.
- **S23:** Six formula cards reduced to three readable statistical stages using saved example data.
- **S24:** Distinguishes any-discovery probability, mean false calls and per-library FDP; interval definitions are stated.
- **S25:** Original Figure 3E–G become A–C; original H is separated into false positives (D) and power (E).
- **S26:** Original Figure 4I read accounting, preserved as a dedicated companion figure.

## Important implementation details

The active build has two stages. `sources/revised/build_data_figures.py` loads the copied figure builders, intercepts their save calls, and adjusts/splits their Matplotlib figures. `sources/revised/custom_figures.py` then supplies replacement diagrams, count matrices, external reanalysis layouts and separate background tracks. Some second-stage figures intentionally replace intermediate first-stage exports.

Inspect both the copied builders and these orchestration files. Looking only at `fig3_calibration.py`, for example, will not reveal the final split/layout.

A read/write guard in `build_data_figures.py` resolves recorded reads to copied inputs, redirects derived JSON outputs into `scratch/`, and rejects other attempted writes outside the revision folder. Some constants and optional-resource checks still assume the existing repository layout; this is not a fully independent software distribution.

Additional scripts:
- `sources/revised/finish_coverage.py`: retained S13 vector curves, denominator annotation and scaled stroke weights.
- `sources/package_figures.py`: PDF pack, gallery, captions, TIFF exports and integrity report.
- `sources/snapshot_inputs.py`: copies the resources actually read.
- `sources/validate_and_archive.py`: print sizes, PDF text bounds, hash checks, numerical cache comparisons and archive creation.

`setup_copy.py` and `revise_sources.py` record initial preparation. Do not rerun them over the finished revision: they initialize/reset copied sources. For normal rebuilds follow the README and run the active builders in order.

Use the project's `.venv/bin/python` for the data-figure builders (Matplotlib 3.11 and SciPy). The packaging scripts use bundled Python with pypdf, reportlab, Pillow and pdfplumber. If rendering with LibreOffice becomes necessary, use the bundled runtime rather than the user's desktop installation; no Word rebuild was performed in this figure task.

## Verification recorded by Codex

The current `validation.json` reports:

- 43 figure/companion pages and a 44-page review PDF.
- 196 protected original files unchanged.
- 70 recorded original data inputs unchanged.
- All 72 copied resources match their source hashes.
- No detected text outside individual figure bounds.
- No detected physical SVG print-minimum violations, including scaled S13 artwork.
- Eleven derived summaries match their existing original caches exactly:
  - `figS2_numbers.json`
  - `figS15_estimator.json`
  - `immcantation_summary.json`
  - `figS7_umi.json`
  - `controls_design.json`
  - `figS15_numbers.json`
  - `fig2_numbers.json`
  - `titration_summary.json`
  - `fig4_tests.json`
  - `fig7_morange.json`
  - `fdr_calibration.json`

Codex rendered the PDF pack and visually inspected the pages, with repeated inspection of corrected panels. These checks establish file preservation, selected numerical parity and presentation checks; they do not prove every plotted element or interpretation is correct. A text-bounds pass does not detect all overlaps or scientifically misleading visual encodings. No new method-validation benchmark was run.

## What is not finished or should be challenged

1. **The Word manuscript is unchanged.** Revised figures have not been integrated into a rebuilt article. Main-text panel references, captions, supplementary ordering and any resulting word count must be updated together if these figures are adopted.
2. **S25/S26 are provisional appended identifiers.** Final supplementary numbering should follow first citation. S14–S16 companion tables do not yet have final submission table numbers.
3. **The revised legend file is an integration aid.** Compare every caption against the actual artwork and full original legend. Check especially removed shading, changed counting displays, separated pages, old panel letters, sample sizes, denominators and the survival of experimental detail.
4. **Some recommendations were only partially implemented.** S5–S7, S10, S12 and S18–S19 retain substantial existing structure. Decide whether the remaining density and labelling are acceptable.
5. **Simplified diagrams may now be too text-heavy or sparse.** Assess Figure 1 and S1/S21–S23 for explanatory value, not merely absence of rounded rectangles. Check whether retaining all these supplements still duplicates material.
6. **Do not equate count-preservation with presentation fidelity.** Verify group inclusion, matrix orientation, pooling, per-library versus per-animal units, pairing identities, geometric-mean floors and the scope of error bars/intervals.
7. **Check S17 fit terminology independently.** Confirm the denominator and selected clone identity from saved tree data/source logic; do not infer that every clone read is represented by a plotted mutation-bearing node.
8. **Review S13's retained vector artwork.** Its styling is not identical to the newly drawn figures. Check whether the denominator clarification is sufficient and whether all physical line/text sizes are appropriate.
9. **No additional inferential tests were introduced.** Correlations remain descriptive; no new agreement analysis, formal FDR proof, or new biological replication claim should be inferred.
10. **Inspect export consistency.** SVG, PDF and raster views should agree; evaluate at 174-mm width, not only a large screen. Figure captions in the review PDF are brief design notes, not complete submission legends.

## Suggested inspection sequence

1. Read the README, figure manifest and integration map.
2. View the entire PDF pack and use the before/after gallery for changes that seem questionable.
3. Inspect active orchestration and copied builder code for Figures 1–7 first.
4. Check saved data/group identities for S9, pairing in Figure 5C/7E, reconstruction denominators in Figure 6/S17, and per-sample composition in S20.
5. Audit legend/panel mapping and information relocation.
6. Reproduce selected figures from the copies, then rerun the verification script if needed.
7. Report issues by severity with exact asset/panel and source-code locations. Separate scientific errors from aesthetic preferences and optional improvements.

## Requested output from Claude

Produce an independent Markdown inspection report with:

- high-priority scientific or interpretive problems;
- visual/layout regressions and publication-size readability problems;
- source/rebuild or export inconsistencies;
- incomplete recommendations and duplicated material;
- concrete fixes, with the affected panels and code paths;
- a recommendation for which revised figures to adopt, revise further or revert;
- a precise list of manuscript integration changes still required.

Preserve existing results and originals. If proposing an analytical change rather than a presentation fix, identify it separately, explain why it is needed, and do not silently replace reported numbers.
