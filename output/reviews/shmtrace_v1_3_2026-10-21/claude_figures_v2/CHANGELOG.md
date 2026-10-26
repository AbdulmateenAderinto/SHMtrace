# SHMtrace v1.3: Claude figure revision 2 (22 Oct 2026)

Revision 2 starts from Claude revision 1 (`../claude_figures_v1/`, see its CHANGELOG for the main-figure redesigns). It adds the corrections Codex made that revision 1 lacked, plus the fixes from `../CLAUDE_INDEPENDENT_REVIEW_OF_CODEX_WORK.md`.

Nothing outside this folder was changed. The canonical sources, the Codex 1 manuscript, `codex_figures_v2` and Claude revision 1 are all intact; `../claude_figures_v1/original_hashes.txt` was re-checked. No analysis result, default or statistical test was changed.

**Deliverables:**
- Manuscript: `SHMtrace_CellReportsMethods_complete_v1.3_Codex1_ClaudeFigures2.docx`
- Rendered PDF: `render/` (90 pp.)
- Figures: `figures/`
- Sources: `src/`
- Manuscript build: `manuscript/revise_operations.py` → `author_docx.py`

**Rebuild:**
1. Run the changed `src/` scripts with the project `.venv`.
2. Run `python manuscript/revise_operations.py`.
3. Run `python manuscript/author_docx.py`.

## 1. Codex corrections adopted

**Already present in revision 1.** Codex's manuscript-level corrections came in through the Codex 1 base:
- citation keys (minimap2, *Mkln1*)
- the 95% recall denominator
- S21 trace accounting and cross-references
- the Fig 2 → 3C pointer
- "two loci in three assay contexts"
- "two datasets" in the engineered-diversification Methods
- the Karst → S8 destination
- the two-to-four-control FDP range
- reference metadata
- the Table S2/S4 title fixes

**New in revision 2.** Codex's v2 figure corrections that revision 1 had not touched. The code was ported from `codex_figures_v2/sources/revised/custom_figures.py` and `build_data_figures.py` into `src/figS_matrices_tracks.py` and `src/figS20_switch_pacbio.py`, using the same data sources and counting:

| Figure | Correction |
|---|---|
| S9 page 1 | N114 density on a linear axis, so the zero-call subclone sits at 0 instead of at a one-event value on a log axis. Count matrices with a shared scale; the one empty-vector call is not read as a spectrum. |
| S9 page 2 | One count matrix per library on one scale. Zero-call libraries are marked "spectrum undefined" and 1–2-call libraries "descriptive only". |
| S11 | Controls on separate tracks sharing one scale per locus, replacing overlays and insets. The read maps are unchanged. |
| S13 | Per-mouse Sanger counts of distinct Q≥20 substitutions on one scale. The AID-null totals of 1–2 are marked descriptive only. |
| S17 page 1 | Counts, plus a nested-window scatter stating that the comparison is descriptive, not independent validation. |
| S20B | One composition bar per mouse instead of a group mean. |
| Fig 3 | "Observed FDP" label on the spike-in panels. |

**Not adopted:**
- **S17 tree selection.** Codex's choice is a regression; see review H1. The original per-mouse trees (Wt2:0, Wt3:0, Wt4:9) are kept.
- **Codex's Fig 1, S1, S3, S21–S23.** They lose content; the originals are kept.
- **Fig 4 and Fig 7A/B zero strips.** These still place zero on a log axis.

## 2. Review fixes (H2–H8 and others)

The FDP terms below are not new numbers; they are the same pooled realised proportions under the right name.

| Item | Change |
|---|---|
| H2: Sanger/nanopore units | Fig 5A labels and legend now say both are distinct-variant densities: distinct Q≥20 changes per clone base versus calibrated calls per read base. Results ¶46 is reworded to match. This was verified: Wt2 135 unique / 18,964 bases = 711.9 per 10^5. |
| H3 | Discussion: "the same AID-null contrast" → "separated the genotypes in the same direction". |
| H4 | Discussion: the 2-fold comparator becomes "2- to 32-fold, depending on germline and masking". |
| H5 and margin meaning | Restored in STAR Methods (not counted): the margin is an empirical factor, and because the same controls choose the margin and supply the decoys, the estimate is conditional, not a guarantee. |
| H6 | Fig 4 legend: "with three libraries per group, ns does not show a caller fails". |
| H7 | Discussion: "approximately twenty" → "about 24, and more if three mice underestimate the variation", matching the S8 legend and `power_s3.json` (0.73 power at 20, 0.91 at 30). |
| H8: FDR vs FDP | Realised values are now "observed false-discovery proportion (FDP)" in the Summary, Results ¶26/28/33, the Discussion and the Fig 3 legend and labels. "FDR" is kept for targets. |
| r = 1.0000 | Discussion → "r > 0.9999" (the Results and S16 were already fixed). |
| S10 wording | Results ¶39: "of discoveries at the largest tested depth" → "of the most each library yielded at any depth", matching the S10 denominator. |

## 3. Checks

- **Word count:** 6,981 of 7,000 (19 under). It was counted with `manuscript/word_count.py`, which uses the project's `word_count()` and reproduces Codex 1's 6,823 exactly. The Summary is 150 of 150.
- **Print check:** `print_check` passes on every new or changed figure: Fig 3, Fig 5, S9 (three pages), S11, S13, S17 page 1 and S20.
- **Rendering:** every corrected sentence was confirmed in the rendered PDF, and the changed pages were inspected.
- **Not run:**
  - `build_crm.py`'s legend check and `renumber_supp.py`, because this revision is replayed from recorded operations.
  - The analytical point in review §6 is untouched: the `v0.1` and `pbinom` ablations give identical numbers and need checking before they are reported, and the naive control-subtraction baseline has not been run.

## 4. Ablation baselines (added 22 Oct)

**Source runs.** `crm/spikein.py` gained two labelled profiles and was run on all four loci as `crm/spikein.py <locus> --profiles recommended,ctlnull_bh,qualonly_pbinom --no-lofreq --tag _ablation`, writing `crm/out/spikein_*_ablation.json`. These runs use the same planted libraries as Figure 3: the recommended mode reproduced the saved spike-ins field for field on all four loci.

**Table S1** (`Table_S1_benchmark.xlsx`, built by `src/make_tables.py`, a copy of `crm/make_tables.py`; the canonical `crm/out` workbook is untouched):
- six new rows on the Spike-ins sheet: "per-position control null, BH" and "read qualities only, no controls" for J_H4, N114 V_H, Sμ and *Pim1*;
- two new notes on the About sheet;
- the "Observed FDR" column renamed "Observed FDP".

**Results ¶28.** One new sentence. Without controls, a read-quality test made 5.7–46.6 false calls per library. A per-position control test made 1.3–2.8 but recalled 27–64%, against SHMtrace's 63–90%.

**STAR Methods.** The two ablations, the command, and the reproduction check.

**SI index.** The Table S1 entry now mentions the ablations.

**Trims to fit the cap.** Each sentence removed is still covered elsewhere:
- Figure 3 legend: the S8 pointer (the Results cite S8);
- Figure 4 legend: the S1E pointer (the Results cite S1E);
- Figure 7 legend: the "grouped by question" sentence (the figure carries the headers).

**Word count:** 6,995 of 7,000 (5 under).

## 5. Author comments on the v2 draft (22 Oct)

The commented draft is kept as `SHMtrace_CellReportsMethods_complete_v1.3_Codex1_ClaudeFigures2_author-comments.docx`. It contained 9 comments and no text edits (diffed against a fresh build). No data or statistics changed.

| Comment | Change |
|---|---|
| Fig 4: unsure about the change in orientation | A–E returned to the original orientation: callers on the x axis, density on a log y axis. The zero-call strip now sits below an axis break, and group sizes stay in the titles. |
| Fig 4: H has some minor clashing | The three series are offset slightly along x for display, so markers and selection rings no longer overlap. The legend has its own space above the data, and the reference label is shortened to "1%: if calibrated". |
| Fig 4: OT1 should be identified as an AID off-target in CH12 | Rows now read "CH12 OT1/OT2". The legend reads "AID off-target loci 1/2 in CH12 cells". |
| Fig 5: bar graphs for C; "AID-null" rather than "null" | C is back to bars (mean ± SD) with mouse-shaped points. "AID-null" is used in the A labels, the J_H4 title and the legend. |
| S1: why split? | No change; this is a question for the author (see chat). |
| S9 (all three pages): add the legend back | The colour-scale legend is restored on every page, with group headings (Wt / *Aicda*^−/−; pMX-*Aicda* / empty vector). S13 and S17 page 1 got the same, for consistency. |
| S17: C lacks r² | C now shows "Pearson r = 0.989, r² = 0.977"; the legend matches. |

**Wrapping and edges.** The Fig 6D subtitle now runs to two lines. An edge scan of all 37 manuscript images found text touching the border in Figs 3, 6 and 7 and in S1E and S6; their margins were widened. One edge contact remains, in S14: a SciCanvas export from the Codex base whose table is complete but meets the right edge; it was left unchanged.

**Checks.** `print_check` passes on every changed figure. Word count: 6,989 of 7,000.

## 6. SciCanvas figures (22 Oct)

The builders are copied into `scicanvas/`; the canonical builders are untouched. Each figure was built with the project's `crm/scicanvas/build.js` into `figures/`, producing an editable `.scifig`, PDF, SVG and 1,000-dpi TIFF.

**Approach.** These were targeted revisions, not redesigns. I judged Codex's sparse redesigns of Figure 1, S1, S21 and S22 a regression, and the author did not want S1 split up.

| Figure | Change | Why |
|---|---|---|
| 1 | A is labelled "(schematic)". The comparator is named as "LoFreq, Geneious rule". "empirical FDR target 0.05" becomes "empirical decoy estimate, target 5%". C step 2 is "Candidates tested", with "set an empirical threshold: 5%". The Sanger step becomes "Count Q ≥ 20 changes". The PacBio input reads "HiFi reads, junctions only". | Critique points 2, 3, 5 and 6. PacBio had no route in C or D, so it is now marked as junction-only. |
| S1, page 2 | New panel F: target vs paralog identity (moved from S23B), with a one-line reading guide. | It keeps unique data that S23 dropped. The page is now 128 mm tall. |
| S3 | The "Also:" file lists and the footer are removed. The six panels and their "Watch for" notes stay. Panel B now says the "estimated false-call share [is] targeted at 5%", replacing "held under 5%". | The critique asked to move file names to the documentation. The reading guidance is useful in a methods paper. |
| S14 | The results table is pulled 20 px in from the right edge. | It touched the image border. |
| S21 | The three prose cards are replaced by a table (unit, span, calling rule, only-here, Wt/AID-null, limits, agreement). It is built on Codex's corrected builder (184/134 trace accounting, pooled-category r, STAR Methods pointer). | The critique asked for a compact table. "count different things" is replaced: both are distinct-variant densities (review H2). |
| S23 | Read preparation and competitive alignment are removed; four statistical stages (A–D) remain. The title is now "How SHMtrace decides what is a call". The page is 126 mm tall. | The critique asked to focus the figure on background, test and threshold. The read fates were already in S1E, and the identity scatter moved to S1F. |

**Not changed:**
- **S22.** The critique said to keep it only if it serves a distinct purpose, but the author had asked for PacBio to be added to it (comment 47), so removing it is the author's call.
- **S15 and S16.** No issues remaining. S16 was already fixed in revision 1.

**Legends:**
- **Figure 1:** "empirical false-discovery threshold" becomes "an empirical, decoy-based threshold".
- **S1:** gains (F).
- **S3:** the file-list sentence is updated.
- **S21:** the bottom section describes the table.
- **S23:** rewritten for A–D.

**Checks:**
- Every SciCanvas build passes the Cell Press check (text 6–8 pt, strokes 0.5–1.5 pt). S21 and Figure 1 report only the pre-existing false warnings (J_H4 subscripts read as panel letters, and line-spacing on Figure 1C3), the same as their unmodified builds.
- `print_check` passes for S1 page 2.
- Word count: 6,990 of 7,000.

## 7. Flow of the main text (22 Oct)

The Codex compression had left many Results and Discussion sentences reading as bare statements. Step 7f of `manuscript/revise_operations.py` (texts in `manuscript/flow_edits.py`) restores transitions and the order of the argument, mostly from the authors' v1.3 wording. It keeps every later correction (FDP terms, the review fixes and the ablation sentence). A check fails the build if any number appears that is not already in v1.3 or the current text.

To stay under the cap, three further changes were made:
- **Step 7g** moves methods detail from the Results into STAR Methods verbatim. This covers read-preparation and paralog details, the simulation design and adjusted Rand index, the separate single-read test, the SW_HEL pooled host, the PacBio junction caveats, tree fit and the four-gamete test, the UMI rRNA benchmark, and the off-target p-values on frequency.
- **Step 7h** applies the cuts the author chose:
  - **B:** the Discussion's "Together, …" summary sentences, plus a shorter closing sentence.
  - **D1:** the Clair3/DeepVariant sentence in the Discussion.
  - **C:** caption phrases that the panels now state themselves.

**Word count:** 6,998 of 7,000.

## 8. V2 – for Alberto review (22 Oct)

**Output.** The current file is `SHMtrace_CellReportsMethods_V2_for_Alberto_review.docx`, with its PDF in `render/` (90 pp.). The title page and footers now read "SHMtrace V2 – for Alberto review". The earlier `…_v1.3_Codex1_ClaudeFigures2.docx`, which holds the author's edits and comment, is kept unchanged as the record.

**Author edits ported from that file.** I diffed it against a fresh build: there were no tracked changes, three edits and one comment.
- **Results, HMCES:** "After inspecting this dataset, we developed a rule that excluded it from threshold selection while retaining it for model fitting…"
- **Discussion:** "turned out to be a harder test than any internal benchmark".
- **AI declaration:** "After using these tools".
- **Comment "is this correct?"**, on the *Mkln1* paragraph. Every number matches `numbers.json` (402 traces; 8.0- vs 7.8-fold; p 0.0038/0.014; 15→29% vs 16→27%; r 0.97 and 0.62), and the study is from the same laboratory. However, the opening "gave similar agreement" (from step 7f) overstated the knockout match. It now reads "agreed on frequency but less closely on the knockout spectrum."

**Authors.** Philip Barbulescu has been added (affiliation 1), after A.A. and before A.M., for the laboratory's Sanger sequencing. He has since moved, so a present-address note was added ("Structural Genomics Consortium, Toronto, ON M5G 1L7, Canada", in the format of affiliation 1; SGC is UHN-affiliated). To pay for its words, two caption phrases the panels already show were removed: Figure 4F "with Wt and double-knockout groups shown separately" (the row labels say so) and Figure 2C "dotted line, Q20". Author contributions now read "investigation, A.A. and P.B."

**Word cap.** These edits added words, so two short items were removed:
- the clause after the Figure 2 example ("the remainder could not be distinguished confidently from background";
- "See also Figure S1." in the Figure 1 caption, since the Results already cite S1.

Word count: 6,998 of 7,000.

**Not updated:** the canonical sources in `docs/calibration/crm/manuscript/` (author list, contributions and all prose since Codex 1). See the note in section 3.

## 9. Figure consistency pass (22 Oct)

**Figure S17, page 2 (mutation trees).** The old page came from `docs/calibration/mutation_tree_figure.py`. It was drawn 13.5 in wide and placed at 162 mm, so it printed at 51%: mutation names at about 3.7 pt and the legend at about 4.3 pt, below the Cell Press 6-pt minimum. `src/figS17_trees.py` redraws it at page width (174 × 222 mm) in the Figure 6D style, using Figure 6's `draw_tree`:
- **Trees.** The same trees, rebuilt from the delivered reads with the same clone choice: Wt2 clone 0, Wt3 clone 0 and Wt4 clone 9. The Wt2 tree and the Wt3 tree (Figure 6D's) are asserted node for node against `figures/mutation_trees/trees.json`.
- **Labels.** Mutation names sit under their ticks, and read counts sit above and to the right of each node.
- **Running it.** `--compute` needs SHMtrace (`.venv`); drawing reads `scratch/s17_trees.json`.
- **Print check.** The page passes.

**Mouse and group names.** The figures now follow the text:
- the J_H4 mice are Wt2–Wt4 and KO1–KO3, and their genotype is *Aicda*⁻/⁻;
- "AID-null" stays the general term for control libraries;
- Sμ WT3 in Figure 2 keeps its deposited sample name, as in the text.

`style.mouse()` maps the data-file spellings (WT2, "Wt 2", "*Aicda*⁻/⁻ 1"). The figures changed are:
- **Figure 4A:** "3 Wt vs 3 *Aicda*⁻/⁻ mice";
- **Figure 5A:** tick labels *Aicda*⁻/⁻;
- **Figure 5C:** Wt2–Wt4;
- **Figure S1E:** J_H4 Wt2–Wt4;
- **Figure S3F:** "Wt2: 323 reads" (SciCanvas, rebuilt);
- **Figure S6:** Wt2–Wt4;
- **Figure S9:** Wt2–Wt4;
- **Figure S11C:** Wt3;
- **Figure S13:** Wt2–Wt4 and KO1–KO3, previously "Wt 2" and "*Aicda*⁻/⁻ 1".

The Figure 5 legend now reads "three Wt and three *Aicda*^{−/−} mice" (step 7j). Its word count is unchanged.

**Typesetting.** *Aicda*⁻/⁻ is italic everywhere, and its superscript no longer has mathtext's operator spacing. The scripts now write `^{{-}/{-}}`. This affects Figures 4–6, S9–S11, S13, S19 and S20.

**Label collisions:**
- **Figure 6A:** "individual Wt mice" moved clear of the lines.
- **Figure 6C:** clone counts moved off the axis line.
- **Figure S4A:** "default 5%" moved to the top of its line, away from "realised = target".

**Checks.**
- `print_check` passes on every matplotlib figure, the new S17 page included.
- Word count: 6,998 of 7,000.
- V2 was rebuilt with `revise_operations.py` and `author_docx.py` (python-docx is in `.venv`), and rendered to `render/` (90 pp.).
- The previous docx and PDF are kept in `scratch/` as `V2_before_figure_pass.*`.

**Addendum (23 Oct).**
- **Figure S9A** shows the N114 separation again: "85×, p = 0.0071". `figS_matrices_tracks.py` computes it with a Welch test on log density, entering zero-call subclones at one event, and asserts that it matches the values Figure 4 and Table S1 report. The linear axis stays. The legend gains one sentence on how the test treats zero-call subclones (step 7j). Supplementary legends do not count toward the word limit.
- **Legend typesetting.** The builder set the D of "V(D)J" in bold, as if it were a panel label. This affected the legends of Figures S3, S9 and S21. `author_docx.py` now requires a panel label to open its own parenthesis, leaving the Codex `docbuild.py` unchanged. Bold panel labels went from 110 to 107: those three cases.
- **Checks.** Word count 6,998; PDF re-rendered (90 pp.). The previous docx is in `scratch/V2_before_vdj.docx`.
- **Figure S9, page 1, restyled to match Figures 5 and 7.**
  - Panel A: points with a darker outline and a geometric-mean bar per vector, the separation in grey above, and a grey subtitle.
  - The two From → To matrices are replaced by one spectrum panel (B). Bars show the pooled shares of the pMX-*Aicda* calls by substitution type, coloured from G:C or from A:T as in Figure 5B, with each type's call count beneath it. Points show each subclone.
  - The single empty-vector call is named but not plotted.
  - The legend is updated (step 7j). The matrix version is kept in `scratch/FigureS9_n114_matrix_version.png`. Page 2, the per-library matrices, is unchanged.
- **S9 page 2, S13 (top) and S17 page 1, restyled the same way.** The From → To matrices are gone from the paper. One helper in `figS_matrices_tracks.py` (`spectrum`) draws every spectrum the same way:
  - bars show the pooled shares of the AID-active libraries, coloured from G:C or from A:T, with each type's count beneath it;
  - circles, squares and triangles show each library;
  - the AID-null libraries are named in red rather than plotted.

  The pages:
  - **S9 page 2** is now the J_H4 nanopore spectrum by mouse. Its N114 half is dropped, because page 1 now shows the same per-subclone spectra as points. Its title is now "J_H4 spectrum by mouse".
  - **S13 top** is the Wt Sanger spectrum (distinct Q≥20 changes, 471 over three mice).
  - **S17 page 1:** (A) is paired bars by type, lighter for the 529-bp window and darker for the 1,141-bp interval, with both counts beneath. (B) is the shares against each other, coloured by base pair; r = 0.989 is unchanged.
  - No text cites individual panels of these figures, so the new panel letters break nothing.

  The legends are updated (`legend()` in step 7j). The word count is unchanged at 6,998, and the PDF now runs to 89 pp. The matrix versions are kept as `scratch/*_matrix_version.png`.
- **Two redundant panels removed (author's call).**
  - **S9 page 2 (J_H4 spectrum by mouse).** Its bars were the pooled Wt shares already in Figure S17A (lighter bars) and on Figure 5B's nanopore axis. The per-mouse pattern is Figure 5C's. S9 is now one page: "(page 1 of 2)" and the pointer to the following page are gone from its title and legend.
  - **S13 top panel (Wt Sanger spectrum by mouse).** Its 12 shares equal Figure 5B's Sanger shares type for type (from `stats.json`). The per-mouse pattern is Figure 5C's, and the knockouts' 1, 1 and 2 changes are in Figure 5A and the text. S13 page 1 now holds the readable-window panel alone, and its legend fits on the same page.
  - **Build and checks.** Both removals are steps in `revise_operations.py`, and `figS_matrices_tracks.py` no longer draws either page. No text cited them. Word count 6,998; the PDF is 87 pp.
- **Figure S2 redrawn at full width.**
  - **Why.** It was drawn 85 mm wide but placed at page width, so it printed at nearly twice its size. It was also the last matplotlib figure with gridlines and no panel letters or titles.
  - **What changed.** The two panels are now side by side, at 174 × 72 mm, with bold titles and grey subtitles. In the legend, "Top" and "Bottom" became (A) and (B). The data and binning are unchanged.
  - **Checks.** The figure passes the print check. Word count 6,998; the PDF is 86 pp. The old version is in `scratch/FigureS2_half_width_version.png`.
- **Figure S12 legend:** one sentence explaining that SHMtrace's 188-fold there is a ratio of mutation frequencies over the V segment, while the 85-fold of Figures 4F and S9A is a ratio of unique-variant densities in the same libraries. Supplementary legend, so the word count is unchanged at 6,998.
- **Figure S13 title:** "(page 1 of 2)" removed. It counted the two images across two pages, and only the readable-window image is left.

## 10. Read accounting as its own figure; supplementary figures renumbered (23 Oct)

**New Figure S2, "Read accounting, related to Figure 1".** It holds what was S1E–F: the fate of every delivered read (A) and the competitive-alignment scatter (B).
- **S1** keeps the workflow schematic (A–D) and its legend.
- **The image.** Its panel letters are redrawn as A and B (`fig4_controls.py`); the file keeps its name.
- **The legend.** The E–F sentences moved to S2, along with the experimental details.

**Renumbering (step 7k, the last step of the build).** The old S2–S24 are now S3–S25, and S1E → S2A, S1F → S2B.
- **What it touches.** Only supplementary figure numbers in a list after "Figure" or "Figures", over the whole manuscript: text, captions, legends, tables and the supplemental index. Tables, the S5F model, the FACSymphony S6 sorter and other papers' "Supplementary Figure S2A" are left alone.
- **Review.** Every change is printed when the build runs, and all were checked.
- **Two drawn references.** S23 (formerly S22) drew "every library (Figure S1)", which now reads "(Figure S2)". It was rebuilt; a pixel diff shows only that character changed. The S24 legend's pointer for read preparation and competitive alignment now reads "Figures S1 and S2".

**Citation order.** As part of S1, the read accounting was cited early. On its own, its first citation fell in STAR Methods, out of order. The Results sentence that opens read preparation now cites "(Figures 1, S1 and S2)", which costs one word: the count is 6,999 of 7,000. First citations now run S1–S22 in order.
- **Still out of order, as before.** S25 (stress tests, formerly S24) is cited in the Figure 3 caption before S23 and S24 (overview and engine, formerly S22 and S23), which only STAR Methods cites. The old numbering had the same inversion.

**Other notes.** Image file names keep their old numbers; only the manuscript numbering changed. The PDF is 86 pp. The previous docx is in `scratch/V2_before_S2_split.docx`.
- **First citations now in order to the end (step 7l).** The stress tests are now S23 (were S25), the workflow overview S24 (was S23) and the calling rules S25 (was S24). The stress-test block moves ahead of the other two in Document S1.
  - **What changed.** Every citation was swapped, and each change is printed by the build. The supplemental index range "Figures S1–S25" counts figures, so it is left alone.
  - **Checks.** First citations run S1–S25 in order. Word count 6,999; PDF 86 pp. The previous docx is in `scratch/V2_before_S23_reorder.docx`.
- **Figure S25 (calling rules) shows all six stages again (author's request: the engine figure is the place for the full detail).**
  - **The figure.** Read preparation (A, read fates per library) and competitive alignment (B, identity to the target against the paralog) are back as stages 1–2. They were switched off in the first revision (`SKIP = 2`) after the critique asked to move routine preprocessing out. The statistical stages are now C–F.
  - **Fixes in the same rebuild.** The page height now follows `SKIP`, and rows are 12 px shorter so the figure (174 × 175 mm) and its legend share a page. Panel A labels the mice Wt2–Wt4, and panel D's legend sets *Aicda*⁻/⁻ in italic. Cell Press check: text 6.0–7.9 pt.
  - **Legend.** It covers A–F, with pointers to Figures S2A and S2B for the same panels, tightened so it fits on the figure's page. The title is now "Rules underlying read preparation and variant calling" (step 7m).
  - **Checks.** Word count 6,999; PDF 86 pp. The four-stage builder and image are kept in `scratch/`.
- **Figure S25 legend:** panels A and B now state that they are "the same data as Figure S2A" and "the same data as Figure S2B". Image-integrity practice is to declare data shown in more than one figure. It is a supplementary legend, so the word count is unchanged at 6,999. Figure and legend still share one page.

## 11. Graphical abstract and highlights for review (23 Oct)

**Graphical abstract.** `scicanvas/crm_ga.js` is a copy of the canonical builder, with PacBio added as a fourth platform. The paper validates four platforms; the 17 October version showed three. PacBio uses Figure 1's drawn glyph (a ring for the HiFi template, a bar for the read), because the icon set has none, and the four platforms are respaced to fit.
- **Build.** It is built with the graphical-abstract check (`GA=1`): 174 × 174 mm, text 10.9–16.8 pt, strokes 0.99–1.97 pt.
- **Numbers.** Its figures ("1–4 calls per control set", r = 0.92) still match V2.
- **Files.** The output is in `figures/GraphicalAbstract*`, with copies at the folder top as `SHMtrace_V2_GraphicalAbstract.pdf` and `.png`. The canonical abstract is unchanged.

**Highlights, eTOC and Motivation.** `manuscript/highlights_v2.py` writes `SHMtrace_V2_Highlights_eTOC_Motivation.docx`.
- **Source.** It starts from the canonical build of 21 October, whose numbers `build_crm.py` filled from the analysis.
- **The one change.** The eTOC now opens "Aderinto et al. present…": with Philip Barbulescu added there are three authors, which is the Cell Press convention.
- **Limits re-checked.** Highlights are 84, 74, 85 and 78 characters (≤ 85); the eTOC is 43 words (≤ 50).
- **Canonical sources.** `crm_text.py` still says "Aderinto and Martin". That fix belongs to porting the V2 edits back before submission.

**Read-me.** `READ_ME_FIRST.md` lists both new files and the current word count.
- **Highlight 2 reworded as a proportion (author: "1–4 calls" means nothing without the background it is out of).**
  - **Wording.** It now reads "AID-null control calls fall 64–90% vs an uncorrected rule; true mutations are kept" (82 characters). This is the reduction the Summary already reports, per assay: 1–4 control calls against 5–20 for the uncorrected Geneious rule on the same libraries.
  - **Source of the numbers.** `highlights_v2.py` fills them from `crm_numbers.compute()`, not by typing them in.
  - **Not yet changed.** The graphical abstract still says "1–4 calls per control set", and the canonical `crm_text.py` still has the old highlight.
- **Graphical abstract, "Low background" panel.** The caption now reads "64–90% fewer control calls / than an uncorrected rule" in place of "1–4 calls per control set", matching the Summary and highlight 2. The panel's plot moves up to make room for the second line, as in the other two panels.
  - **Number check.** `highlights_v2.py` now checks the abstract's typed numbers (the 64–90% and r = 0.92) against `crm_numbers.compute()`, as the canonical build does for its own abstract.
  - **Checks.** Cell Press check passes. The copies at the folder top are refreshed.
- **Highlights 3 and 4 (author, 23 Oct): SHMtrace is SHM-specific as well as error-aware.**
  - **Highlight 3** is now "Calls feed SHM readouts: AID hotspot targeting, spectra, clones and mutation trees".
  - **Highlight 4** is now "AID signatures are recovered from Sanger, nanopore and Illumina with low background".
  - **Why not PacBio.** It is left out because its dataset is class-switch junctions, not AID signatures. "Low background" rather than "suppressing error", because the Sanger route reports background rather than suppressing it.
  - **What moved out.** The re-analysis of three published studies, and the old "whole-amplicon maps", remain in the Summary, eTOC, Results and Figure 6A.
  - **Checks.** Lengths 84, 82, 82, 83.
- **Motivation rewritten (author, 23 Oct): the gap is the shape of the software, not only the cloning.** It now says that the somatic hypermutation software built since has followed B-cell receptor repertoire sequencing, starting from assembled repertoires and summarising mutations already identified, so a defined reporter locus read in bulk on any platform and judged against its own background has had no equivalent modern route; nanopore reads plus control-calibrated calling give it one. 170 words, from 157.
  - **Checked before writing.** Sible et al. 2021 is the current published J_H4 protocol (cited as reference 10); Immcantation, IgBLAST, SHazaM and MiXCR are repertoire tools, which the Discussion and Figure S13 already say. It reads "has followed", not "is only for", because MiXCR documents a nanopore preset.
  - **Deliberately not said.** Nothing about any tool's availability or uptime. SHMTool's original URL (`scb.aecom.yu.edu/shmtool`) no longer resolves and the tool has moved to `shmtool.einsteinmed.edu`; the outage could not be dated, since the Internet Archive was down. Worth confirming the live address before submission if the paper or STAR Methods prints a URL.
  - **Not counted.** The Motivation is in the separate submission file, so the 6,999-word main text is unaffected.
- **The Motivation's framing carried into the main text (step 7n).** The gap is not only that cloning is slow: the software built for SHM since has followed B-cell receptor repertoire sequencing, so a reporter locus read in bulk has had no equivalent route.
  - **Summary.** The opening sentence is now "Somatic hypermutation (SHM) is measured accurately but at low throughput by cloned-amplicon Sanger sequencing, while existing software targets B-cell receptor repertoires."
  - **Introduction.** The last sentence of the "two common analysis practices" paragraph now reads "SHM software has instead followed B-cell receptor repertoire sequencing, summarising mutations already identified… rather than using matched AID-null controls to separate candidate variants from technical background." The same three citations stay. It reads "has instead followed", not "is only for", because MiXCR documents a nanopore preset, and says nothing about any tool's availability.
  - **Discussion unchanged.** Its gap paragraph ("Repertoire toolkits… begin, however, from mismatches already taken to be real") now follows from the Introduction rather than arriving cold.
  - **Paid for by eight trims** in the same two sections, none of which changes a number or a claim: spectra/linkage phrasing, "four sequencing platforms" → "four platforms", "showed … separation" → "separated", "extended analysis to" → "covered", "proportions of" → "proportions", "mismatches that recur" → "mismatches recurring", "What all platforms, however, share is the need" → "All platforms, however, share the need" (the transition is kept), and "replicates, starting from delivered reads" → "replicates, from delivered reads".
  - **Checks.** 6,999 of 7,000 words; Summary 149 of 150; Introduction 488. PDF 86 pp. The previous docx is in `scratch/V2_before_maintext_framing.docx`.
- **Summary describes the comparator instead of naming it (author, 23 Oct).** "relative to the uncorrected Geneious rule" is now "relative to an uncorrected rule from a general-purpose sequence-analysis platform not tuned for SHM". "GUI" was dropped on reflection: SHMtrace ships a desktop app of its own, so calling a comparator GUI-based would read inconsistently; the contrast is general-purpose versus SHM-specific. One redundant word in the Introduction ("In the configurations evaluated here") pays for it. The Results, Discussion and Table S1 still name Geneious (18 mentions in the rendered PDF), so nothing is hidden; the abstract makes the general point and the body gives the specifics. The description matches that tool's own cited title, "an integrated and extendable desktop software platform for the organization and analysis of sequence data"[@kearse2012].
  - **Paid for by four Summary trims:** "mutations, recalled at observed false-discovery proportions 0.050 and 0.045," → "mutations (false-discovery proportions 0.050 and 0.045),"; "the stated zero-call convention" → "the zero-call convention"; "compare mutation profiles" → "compare profiles"; "preserved mutation linkage" → "preserved linkage"; "Further datasets covered" → "Datasets covered".
  - **Checks.** 6,999 of 7,000 words; Summary 149 of 150; PDF 86 pp. The previous docx is in `scratch/V2_before_abstract_comparator.docx`.

## 12. Corrections after the Alberto review (25 October 2026)

Four claims in V2 were wrong. The measurements behind each are in
`docs/calibration/CALLER_FIXES.md`; the Ramos numbers come from
`docs/calibration/crm/out/public/hmces` and `.../hmces_het`, the UMI numbers from
`docs/calibration/umi_karst/results.json` after the UMI grouping fix. Only the first costs
counted words.

| Where | Was | Now |
|---|---|---|
| Results, Ramos (Figure 7D) | "Control assignment thus materially affects the inferred substitution burden here", and the library called "noisy" | The excess is allelic: every AID-null clone carries a second IGH allele near half its reads, which the 50% germline filter retains and phasing stripped from two clones but not the third. Masking those positions brings all twelve clones within 0.21–1.66 per 10³ bases, from a bimodal 0.4–0.8 and 3.8–8.2 |
| STAR Methods, margin rule | the origin clone's excess left unexplained | traced to that second allele rather than to sequencing noise |
| STAR Methods, public Illumina | "flat Q30 ... so for them the error model has a single quality level" | flat Q30 as distributed, but merging writes Q60 where both mates agree, so the model sees two levels; measured on the merged reads PRJNA973112 is 93% Q60 and PRJNA918596 88% Q30, and the single-read estimate flags both |
| STAR Methods and Figure S9H | UMI consensus "57-fold (95% CI 35–99) to 0.17%", 103 molecules | 46-fold (95% CI 30–83) to 0.21%, 105 molecules |

**`mask_control_het` added to the STAR Methods** high-allele-fraction paragraph: its thresholds, that it is
off by default so no result reported here uses it, and what it removes set by set (29 sites in the Ramos VDJ,
2 in the public nanopore Sμ CH12 libraries, none in the other eight, this study's own two among them).

**Figure S9H was redrawn** (`src/figS8_downstream.py` against the current `umi_karst/results.json`), so the
panel and its legend agree. Page 1 is unchanged: its numeric labels are identical.

**The 277–413 sensitivity analysis was removed** from the Results. It was reported to support the sentence
above it, and with the cause known it is a dangling number. Its removal pays for most of the explanation.

**Word count: 7,022, which is 23 over.** Results went 3,712 → 3,735. The correction replaces a false
causal claim and cannot be made much shorter without dropping either the mechanism or the 52–99 default
result; 23 words still need to come from somewhere in the counted text.

**Preserved:** `SHMtrace_CellReportsMethods_V2_as-sent_2026-10-23.docx` is the document as Alberto received
it, byte-for-byte.

**Back under the limit (step 12 in `revise_operations.py`).** Three cuts in the Results remove exactly the 23
words: each repeats something the paper says elsewhere, so no number, claim or caveat is lost.
- S5F targeting: "ranged from 1.22 to 1.36, indicating 22–36% enrichment over the window average" → "ranged
  from 1.22 to 1.36" (an index of 1.22–1.36 is itself 22–36% enrichment).
- The ≥50% filter: "Excluded sites are therefore reported for inspection." removed; the filter's introduction
  says "every excluded site is reported", and the Limitations say "Excluded sites are reported for review".
- "The desktop and command-line interfaces share one analysis engine." removed; the Discussion says "one engine
  behind the desktop and command-line interfaces".

Word count 6,999 of 7,000; Summary 150 of 150; Results 3,712. The PDF is re-rendered (87 pp.). The docx before
these cuts is in `scratch/V2_before_limit_cuts_2026-10-25.docx`.

