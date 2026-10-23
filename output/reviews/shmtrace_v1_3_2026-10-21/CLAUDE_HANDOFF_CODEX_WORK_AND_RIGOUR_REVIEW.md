# Claude handoff: independently inspect Codex's SHMtrace work

Date: 22 October 2026

## Request

Inspect what Codex changed in the SHMtrace manuscript and figures. Assess the changes independently, including whether the scientific claims remain supported and whether the level of rigour is appropriate for a small methods paper forming one contribution among several PhD papers.

Do not treat Codex's favourable assessment of doctoral-level rigour as an established finding. It was a provisional judgment based on the available project context and review, not a complete independent audit of the analyses. Do not assume that numerical-preservation checks prove scientific correctness or that redesigned figures are necessarily better.

Start with a written inspection report. Preserve the canonical manuscript, original figures and the Codex versions. If making subsequent corrections, create another clearly named version.

## Repository and entry points

Repository:
`/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2`

Review directory (all paths below are relative to this directory unless stated otherwise):
`/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21`

Read these first:

1. `CLAUDE_INSPECT_CODEX_FIGURE_REVISION.md` — detailed figure inspection handoff, including changes, code paths, checks, limitations and proposed inspection sequence.
2. `codex_figures_v2/README.md` — active build instructions and scope.
3. `codex_figures_v2/FIGURE_INTEGRATION.md` — panel relocation map and caption guidance.
4. Repository-relative `docs/calibration/crm/SESSION_HANDOVER_2026-10-21.md` and `docs/calibration/crm/HANDOVER.md` — canonical project context. Historical numbering and word-count statements may disagree; establish the actual current state from files.

## Deliverables to inspect

| Path | Contents |
|---|---|
| `codex_version/SHMtrace_CellReportsMethods_complete_v1.3_Codex_1.docx` | Separate Codex prose revision of the complete manuscript |
| `codex_version/document_operations.json` | Construction operations containing the applied manuscript text |
| `prepare_codex_version.py` | Manuscript preparation, prose edits, detail relocation and initial figure corrections |
| `author_codex_docx.py` | Word authoring and layout adjustments |
| `PROPOSED_REVISION.md`, `REVIEW.md`, `DETAIL_RELOCATION_LEDGER.md`, `DISCUSSION_ENDINGS.md` | Earlier editorial review and proposed changes; not substitutes for inspecting the final DOCX |
| `codex_figures_v2/SHMtrace_Codex_Figures_v2.pdf` | Separate 44-page figure review pack |
| `codex_figures_v2/index.html` | Before/after figure gallery |
| `codex_figures_v2/figure_manifest.json` | Ordered figure/companion-page manifest |
| `codex_figures_v2/figures/` | Revised SVG/PDF/PNG/TIFF assets |
| `codex_figures_v2/sources/revised/` | Active figure builders and layout orchestration |
| `codex_figures_v2/sources/originals/` | Copied pre-revision figure sources |
| `codex_figures_v2/inputs/` | Copied saved analysis resources used for the figure revision |
| `codex_figures_v2/revised_legends.json` | Revised caption drafts |
| `codex_figures_v2/manuscript/original_legends.md` | Full original legends retained for detail and comparison |
| `codex_figures_v2/validation.json`, `CHECKS.md`, `print_check.json` | Recorded checks, with their limited scope |
| `SHMtrace_Codex_Figures_v2_sources_and_exports.zip` | Sources, copied inputs and exports; TIFFs remain outside the ZIP in the local figure folder |

Important: the Word manuscript predates the larger figure revision. It has **not** been reassembled with the v2 figures. The manuscript prose revision and figure revision are separate artifacts. Do not assume their panel references or captions are already synchronized.

The copied `codex_version/sources/` manuscript files do not necessarily contain all applied prose edits. Inspect `document_operations.json`, the preparation script and the final DOCX.

## Work performed by Codex

### Manuscript revision

Codex created a separate complete manuscript with more concise Results and legends, retained methodological detail in STAR Methods, revised Discussion endings, reference/citation corrections and initial figure/layout fixes. The reported counted length was 6,823 words under the existing project counter. Recalculate this and verify what the counter includes rather than treating the number as proof of journal compliance.

Inspect particularly:

- whether compression removed necessary qualifications or changed scientific meaning;
- whether relocated material remains accessible and is actually present in the final document;
- whether positive Discussion endings remain proportionate to the evidence;
- reference metadata, citation matching and figure/table links;
- numerical placeholder bindings and consistency with source outputs;
- the distinction between mismatches, candidate variants and called mutations;
- denominators and the distinction between reads, molecules, variants, clones, libraries and animals.

### Figure revision

Codex simplified the workflow illustrations, restructured the main figures, increased emphasis on biological observations, simplified count matrices and separated some dense comparison tables and tree examples into companion pages.

Major panel moves:

- Figure 3E–G → provisional S25A–C.
- Figure 3H → S25D (false positives) and S25E (power).
- Figure 4I → provisional S26A.
- S8 estimation/power and UMI panels occupy separate pages under the same figure number.
- S14–S16 numerical comparison tables occupy companion pages and have editable CSVs.
- Figure 7E is split into E1 (calibrated active-AID calls) and E2 (original Q > 20 mismatch summaries).
- S17 selected trees are WT2 clone 0, WT3 clone 0 and WT3 clone 14, not one example per mouse.

See the detailed figure handoff for the full list. S25/S26 are review identifiers; final supplementary ordering and any table numbers remain to be settled during manuscript integration.

The active figure build runs copied builders through `sources/revised/build_data_figures.py`, then applies custom replacements through `custom_figures.py`. Inspect both stages. Some first-stage outputs are intentionally replaced by the second stage. Follow the README; do not rerun initialization scripts over finished copies.

## What verification established

The recorded final figure checks report:

- 196 protected original files unchanged;
- 70 recorded original inputs unchanged;
- all 72 copied resources matching their recorded hashes;
- 43 figure/companion pages in a 44-page review pack;
- no detected text beyond figure boundaries;
- no detected physical SVG print-minimum violations;
- 11 derived numerical summaries matching their original cached counterparts.

The PDF pack was rendered and visually reviewed. The checks support preservation and selected numerical parity. They do **not** establish correct biological grouping, all pairing identities, full inferential validity, adequate independent validation, novelty, fair caller comparisons or transfer to unseen preparations. Inspect those separately.

No new method-validation benchmark was run for the design revision. Calling defaults and statistical test definitions were not intentionally changed.

## Scientific-rigour questions

Assess the work as a bounded methods contribution, not as an entire doctorate in one paper. Separate essential validity requirements from optional expansion.

1. **Central contribution:** What specifically does SHMtrace add beyond existing approaches? Is that contribution supported without overgeneralizing about competitors?
2. **Calibration:** Do the sample caller, calibration procedure, decoys, saved panels and comparator harnesses use compatible eligibility and error rules? Identify differences that are intentional versus problematic.
3. **Scope of false-discovery claims:** Keep observed FDP, empirical decoy calibration and formal FDR guarantees distinct. Verify that figures and prose use the appropriate terms.
4. **Independent validation:** Distinguish planted events, resampling existing controls, simulations, matched Sanger comparisons and external biological datasets. Establish whether independent preparations/runs are available and what transfer claims they support.
5. **Biological replication:** Verify experimental units, group sizes, pairing, clone dependence, pooled summaries and treatment of repeated reads. Check that precision is not exaggerated by treating dependent observations as independent replicates.
6. **Benchmark fairness:** Check reference choice, depth, quality filtering, caller settings, calibration route, controls and mutation fractions. Comparisons should support the stated question rather than a broader untested conclusion.
7. **Measurement definitions:** Verify unique-variant density, read-weighted frequency and clone-weighted frequency; check depth dependence and platform-specific counting units.
8. **Agreement:** Correlation of compositional categories or nested windows is descriptive. Check for claims of absolute agreement or independent validation that these plots cannot establish.
9. **Limitations:** Assess case-only artifacts, control mismatch, limited control number, small sample sizes, quality simplification and ambiguity handling where relevant to the current code/results. Verify current status rather than assuming older notes remain accurate.
10. **Claim–evidence alignment:** For each central Results/Discussion claim, identify the supporting output and whether the evidence is direct, simulation-based, descriptive or inferred.

## Figure concerns requiring independent judgment

- Some new diagrams may be too sparse or text-heavy; absence of rounded cards is not itself better scientific communication.
- Several supplementary layouts retain much of their original structure; some recommendations were only partly implemented.
- Removing percentages from matrices aids readability but may obscure a useful comparison unless denominators and corresponding summaries remain available.
- Reconstructed tree relationships and fit denominators need checking against the saved data and construction logic.
- S14–S16 use numerical data already frozen in existing SciCanvas builders, including their existing precision; distinguish this from a new full-precision reanalysis.
- S13 coverage curves were retained from saved vector artwork, with denominator clarification and stroke normalization; raw traces were not reprocessed.
- A text-boundary check cannot detect every overlap, poor hierarchy or misleading visual encoding. Compare exports at the intended 174-mm width.
- Caption drafts may retain stale references or descriptions. Inspect every changed panel against both the revised and original full legends.

## Requested output

Create an independent Markdown report containing:

1. High-priority scientific errors or unsupported interpretations, with exact evidence and file/panel references.
2. Manuscript editing regressions, missing qualifications and reference inconsistencies.
3. Figure/design regressions, export inconsistencies and reproducibility problems.
4. A claim–evidence table for the central methodological claims.
5. A decision for each main figure: adopt, revise further or revert, with concrete fixes.
6. Essential remaining validation versus optional additions, proportionate to a small PhD methods paper.
7. A manuscript integration checklist covering panel references, legends, supplementary ordering, word count and final rendering.
8. A calibrated overall assessment: what is established, what remains uncertain and whether any unresolved point materially weakens the paper.

Do not silently change reported numbers, analysis defaults or statistical tests. Identify analytical changes separately from presentation fixes. Preserve all originals and clearly label any new version.
