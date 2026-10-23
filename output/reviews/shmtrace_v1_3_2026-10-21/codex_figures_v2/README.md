# SHMtrace figures: Codex revision 2

Created as a separate revision of the figures reviewed with the author. Original manuscript, sources, exports and saved analysis files are preserved.

## Start here

- `SHMtrace_Codex_Figures_v2.pdf`: sequential review pack, figures at 174-mm print width.
- `index.html`: before/after gallery with individual downloads.
- `FIGURE_INTEGRATION.md`: panel moves, revised captions and manuscript integration instructions.
- `figures/`: vector PDF and SVG, 600-dpi PNG and TIFF exports.
- `sources/revised/`: editable builders and the active layout code.
- `sources/originals/`: copied pre-revision source files.
- `inputs/`: copies of the specific saved analysis resources used, not raw sequencing archives.
- `sources/*_data.json`: frozen data from existing SciCanvas builders.
- `sources/*_comparison_table.csv`: editable numerical comparison tables.
- `validation.json`: original and input file integrity checks.

## Main changes

1. Replaced the rounded-card workflow with an open layout, keeping the different Sanger/read-library/PacBio routes explicit.
2. Reorganized Figure 2 around the paired candidate examples, with a full-width quality comparison.
3. Reduced Figure 3 to calling performance. Preserved its inference simulations in new Figure S25; separated false-positive calibration from power.
4. Enlarged Figure 4's small multiples and added colour scales and zero-call strips. Preserved read attrition in new Figure S26.
5. Made Figure 5's mouse-level platform comparisons paired points; labelled motif availability explicitly.
6. Separated clone composition from weighting in Figure 6 and clarified tree interpretation.
7. Split Figure 7E into calibrated-call and raw-mismatch tracks with experiment identities.
8. Simplified workflow supplements, separated S8's UMI panels, standardized count matrices, enlarged selected trees, separated S14–S16 numerical tables, and displayed PacBio composition per sample.

The original Word manuscript has not been overwritten or silently renumbered. This is a figure revision package, not a submission-ready reassembled manuscript. S25/S26 are provisional appended identifiers; use the integration map when adopting the figures. Full original legend details remain archived in `manuscript/original_legends.md`.

## Scientific boundaries

- No sequencing data were regenerated and no caller settings, statistical test definitions or thresholds were changed.
- Descriptive plotting summaries are recomputed from saved outputs, as in the original builders.
- The source statistical comparisons and results are preserved. No new inferential tests or agreement analyses were added.
- The S13 coverage curves remain the existing saved vector curves, with a denominator clarification; raw traces were not reprocessed.
- S14–S16 reanalysis plots reuse the numerical data frozen in the original SciCanvas builders, including their existing precision.
- S17's three selected trees are WT2 clone 0, WT3 clone 0 and WT3 clone 14. Their fit annotations refer to mutation-bearing reads, not every read in a clone.
- Empirical decoy calibration is not represented as a formally proven general FDR guarantee.

## Rebuild

Use the repository's `.venv/bin/python` for Matplotlib 3.11 and SciPy. Run from the repository root:

```sh
.venv/bin/python -W ignore output/reviews/shmtrace_v1_3_2026-10-21/codex_figures_v2/sources/revised/build_data_figures.py
.venv/bin/python -W ignore output/reviews/shmtrace_v1_3_2026-10-21/codex_figures_v2/sources/revised/custom_figures.py
```

The first command runs copied builders with revised layouts. The second applies the custom diagrams, count matrices, reanalysis layouts and separate control tracks. Run both in this order: some custom figures intentionally replace the first command's intermediate exports.

The read guard resolves recorded inputs to their copied snapshots and redirects derived JSON writes to `scratch/`. It rejects other writes outside this revision folder. The code still expects this repository layout for locating modules and checking optional resources.

`finish_coverage.py` and `package_figures.py` use bundled Python with pypdf, reportlab and Pillow. For the coverage PNG, render its vector PDF with bundled `pdftoppm -r 600 -png -singlefile` before packaging. `package_figures.py` creates the PDF pack, gallery, captions, TIFFs and integrity report.

`setup_copy.py` and `revise_sources.py` record initial preparation; do not rerun them over a finished revision, since they initialize/reset copied builders. Normal iteration uses the active revised builders above.

The portable ZIP contains figures, sources, the gallery, inputs and documentation. High-resolution TIFFs remain in the local `figures/` folder and are omitted from the ZIP to keep it smaller. No raw sequencing reads or third-party mutability-model files are bundled.
