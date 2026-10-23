# SHMtrace v1.3 — Claude figure revision 1 (22 Oct 2026)

A separate revision built on top of the Codex manuscript (`codex_version/…_Codex_1.docx`). Nothing in
`docs/calibration/crm/`, `codex_version/` or `crm/out/` was modified (372 files hash-checked against
`original_hashes.txt`). `codex_figures_v2/` belongs to a concurrent Codex session and was never touched.

## Deliverables

| What | Where |
|---|---|
| Revised manuscript | `SHMtrace_CellReportsMethods_complete_v1.3_Codex1_ClaudeFigures1.docx` |
| Rendered PDF (89 pp.) | `render/SHMtrace_CellReportsMethods_complete_v1.3_Codex1_ClaudeFigures1.pdf` |
| Figure exports (PNG 600 dpi, PDF, SVG, TIFF) | `figures/` (unchanged figures are copies of the Codex assets) |
| Editable figure sources | `src/` (copies of `crm/figs/*.py`; output and scratch writes redirected here) |
| Manuscript edits, reproducible | `manuscript/revise_operations.py` → `manuscript/document_operations.json` → `manuscript/author_docx.py` |
| Before images | `before/` |

Rebuild: run the changed scripts in `src/` with the project venv, then
`python manuscript/revise_operations.py && python manuscript/author_docx.py`.

## Do the critiques make sense?

Mostly yes, and the main-figure critiques were implemented. Checked against the data, three of them also
turned up real defects (see "Defects found"). Where I departed from a recommendation:

- **Figure 3 split (agree; implemented differently).** The handoff suggests a new supplementary figure. A new
  number would renumber S9–S24 (the 21 Oct renumbering touched 86 files). Instead, old 3E–3H join the
  single-read and power panels as **Figure S8, page 1** ("Downstream estimates and group tests"). The UMI panels
  become **page 2**. The paper already uses two-page figures (S9, S13, S17), and supplementary citation order is
  preserved (S7 → S8 → S9). **Author decision:** whether UMI should get its own number.
- **Figure 4I → supplement (agree).** It is now **Figure S1E** (page 2 of S1, beside read preparation). The
  Results sentence "Read preparation removed 67%…" stays.
- **Figure 1 and the SciCanvas workflow figures (S1 workflow page, S3, S21, S22, S23): not done.** The critique
  is sound, especially folding S22 into Figure 1/S1 and replacing the S21 prose cards with a table. But these
  are design decisions the authors should make first, not a reformatting job. Next step.
- **Supplements not revised in this pass:** S2 (its legend already defines shading and the count panel), S6 and
  S10 (both already use a separate zero strip), S9, S11, S13–S15, S17–S20. Their critiques are presentation-only
  and still open.

## Main figures

**Figure 2.** A and B are now one matched pair on one axis, with the verdict (called / not called, p) in
the panel. C was a separate four-group plot; it is now two quality strips aligned under A and B, titled with
the message (WT3 reads have similar quality at both sites; the controls differ). D has direct labels instead
of a legend, a marked threshold point, and the title "Empirical threshold". The legend says it is not a
guaranteed FDR. E is a clean histogram with its denominator (1,197 calls, 12 libraries). The comparator
shading and text were removed and moved to Figure 3C.

**Figure 3.** Four panels. A–B are unchanged spike-in data, with SHMtrace heavier and comparators lighter.
C now overlays the cumulative allele-fraction distribution of real SHMtrace calls, which carries the
comparator message that was in 2E. D is enlarged to full height, with the linear-0–0.1/log-above axis labelled
on the panel and a direction cue.

**Figure 4.** A–E are small multiples with one shared method-label column. Libraries with no call sit in a
shaded "0" strip left of an axis break, not at a fake small value. Group sizes are in the titles, and each
panel states that its axis range differs. F and G have colour bars (G now has one row per dataset, no gaps),
and the legend says colour is descriptive. H is larger, shows the observed held-out rate against the
"1% expected if calibrated" reference, rings the selected margin and gives the number of tests. I has moved
to S1E.

**Figure 5.** A names each platform's counting unit under its group, with the y label "per 10^5 bases
analysed". C puts the mice in the foreground: paired Sanger→nanopore points per mouse, with mean±SD drawn
lightly behind. D has thinner lines and smaller markers; the terminal-bin note moved to the legend. E labels
its reference lines as the share of eligible positions in each motif. The legend frames the claim as
"similarity of pattern, not interchangeable values".

**Figure 6.** A uses light mouse lines, a dominant mean and direct labels, plus an arrow marking the interval
read only by nanopore. B defines the targeting index on the panel and gives positions per quintile; the
legend says the five points are bins, not replicates. C separates clone composition from a small
read- vs clone-weighted plot beside it. D has a reading guide (ancestor, ticks, node numbers, branches).

**Figure 7.** Question headers group the panels: burden (A–B), signature (C), indels (D), repair (E).
Dataset accessions are in the panel titles. A and B use the shared zero strip; C prints call counts inside
the bars; D has larger points on lighter bars and the y label "indels called per clone library", with every
library at 20,000 reads. E has direct line labels and states what each line joins.

## Supplementary figures revised

- **S4:** emphasised equality line, rings at the 5% operating point, a lighter over-target area, and an axis
  break for the 0 row.
- **S5:** panel letter A restored (it was missing). The one-control BH fallback is shaded and labelled apart
  from the decoy-calibrated conditions. "clean/off/on" are spelled out.
- **S7:** the below-floor region covers both 1 and 2 reads (previously only 1) and is labelled "outside the
  calling rule". No zero-recall points are drawn there. The bar numbers are named as planted counts.
- **S8:** see Figure 3. The group-test panel separates false-positive rate (no true difference) from power.
- **S12:** absolute values come first. A row under the axis gives each route's reference, the bases it keeps
  and what it counts. Zero calls go in a strip. The ratio panel is labelled secondary and greyed except for
  SHMtrace.
- **S16** (SciCanvas, rebuilt from a copied builder): "r = 1.0000" → "r > 0.9999", and "Claims" →
  "Published vs reanalysed". The Results text was changed to match.
- **S24:** defines both outcomes (any discovery; per-library FDP = false ÷ all calls). The case-only column is
  shaded as the boundary of control-based calibration, with n per setting and the CI definition on the panel.

## Relocation list

| Was | Now | Citations rewritten |
|---|---|---|
| Fig 3E (lineages) | Fig S8A | Results ¶ lineages; STAR Methods (high-allele-fraction rule); S8 legend |
| Fig 3F (clone-weighted) | Fig S8B | Results |
| Fig 3G (quantification) | Fig S8C | Results |
| Fig S8A (single-read) | Fig S8D | Results |
| Fig 3H (group tests) | Fig S8E | Results; STAR Methods ×3 |
| Fig S8B (power) | Fig S8F | Discussion |
| Fig S8C (UMI simulated) | Fig S8G | — |
| Fig S8D (UMI real) | Fig S8H | Results; Table S6 |
| Fig 4I (read fates) | Fig S1E | Results; S1 legend |
| Fig 2E comparator text | Fig 3C overlay | Results now cites "Figures 2E and 3C" |

All moves were applied in one pass through placeholders, so no citation was remapped twice
(`manuscript/citation_remap_log.json`). Supplementary first-citation order is unchanged.

## Defects found and fixed (display only; no analysis result changed)

1. **Fig 7A/B mean bars** used fixed display floors (0.01, 0.005) for zero-call libraries, while the legend
   and the paper's test rule enter zeros at one event on each library's own exposure. The bars now follow the
   stated rule; for groups containing zeros they move down (e.g. AID-null 0.017 → 0.007 per 10^3 bp). No
   reported number changes.
2. **Fig 6D subtitle** could be read as "21 reads not drawn". Node counts include those reads
   (`mutation_tree.build` counts every carrier). The subtitle and legend now say they carry mutations of two
   branches and are counted at each.
3. **S7** shading covered only the 1-read bin; the 2-read bin is also below the 3-read rule.
4. **S5** was missing its panel letter A.
5. **S12 legend:** the asterisk in `IGHV1-47-36*01` opened an italic run over the rest of the legend in the
   Codex docx. Escaped.

## Checks

- `print_check` (Cell Press: text ≥ 6 pt, strokes 0.5–1.5 pt): every revised matplotlib figure **ok**
  (Figure 2 initially failed at 0.4-pt bar edges, which were fixed). The SciCanvas S16 build passed its own
  Cell Press check.
- Every revised figure was inspected as a PNG, and the main and changed supplementary pages were inspected in
  the rendered PDF.
- **Word count (estimate):** about 6,960 / 7,000. Legends went from 879 to 1,009 words, and the text gained 7.
  This uses the same tokenisation that reproduces Codex's 879; it was not run through
  `build_crm.journal_word_count`. There is only about 40 words of headroom.
- Not run: `build_crm.py` legend check and `renumber_supp.py check`. This revision, like Codex's, is replayed
  from recorded operations and does not go through the canonical build. Citation order was checked by hand
  (see above).

## For the authors to decide or confirm

1. Whether the UMI benchmark (S8 page 2) should get its own number (that would renumber S9–S24).
2. **Fig 5A wording:** I describe the Sanger value as "Q≥20 substitutions counted across clones per 10^5 bases".
   The script that writes `stats.json["density_sanger"]` was not located. The values (Wt3 ≈ 1,121 per 10^5)
   look like distinct variants rather than every substitution, so please confirm the unit.
3. **Fig 7A text vs display:** the text's 56-fold is a ratio of arithmetic means (`build_crm.ldha_fold`), while
   the bars are geometric means. This is not an error, but the reader sees two summaries.
4. Figure 1 and the workflow figures (S1 page 1, S3, S21–S23): the redesign in the critique is not done here.
5. Figures 4 and 6 grew in height (150 and 134 mm; Cell Press maximum is 200 mm).
