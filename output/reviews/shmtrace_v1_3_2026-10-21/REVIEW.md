# SHMtrace v1.3 manuscript review

21 October 2026. Reviewed against the current working tree above `ef9240d` and the accepted text of the current v1.3 Word manuscript.

The paper has a coherent contribution: SHMtrace connects control-informed mutation detection, linked long-read profiles, and comparisons across biological replicates in a reproducible workflow. The matched-mouse comparison and independent sequence reanalyses give that contribution practical substance. The most useful next revision is to make that argument easier to follow, correct remaining source and figure inconsistencies, and move supporting detail out of the counted narrative.

The current project count is **9,139 words**. The separate editorial proposal is **6,803 words**, with a Discussion of **1,148 words versus 1,145 currently**. The proposed reduction comes from Results and legends, without shortening the Introduction, Limitations or Resource availability. No manuscript source, Word document, analysis result or figure was changed in this review.

## Review files

- [Proposed revision and STAR Methods additions](/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21/PROPOSED_REVISION.md)
- [Six proposed Discussion endings](/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21/DISCUSSION_ENDINGS.md)
- [Detail relocation ledger](/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/output/reviews/shmtrace_v1_3_2026-10-21/DETAIL_RELOCATION_LEDGER.md)

The proposal is reviewable text, not a rebuilt submission manuscript. Its source-preserving JSON and Python files are in this folder. All removed Results placeholders have a proposed Methods destination or an existing occurrence elsewhere; their original paragraphs are retained in the ledger. That is an accounting check, not a substitute for reading the destination in context before implementation. Hard-coded supporting details are also recoverable from those original paragraphs.

## What was verified

The current article DOCX matches the generated source text exactly after accepting changes and normalizing whitespace; the corresponding article portion of the combined v1.3 DOCX also matches. Neither file contains tracked insertions/deletions or comments. There are no newly discovered Word-only edits to port before a later rebuild.

The current DOCX was rendered afresh into an 85-page review PDF, because the PDF already in the build folder predates the latest DOCX. Visual review covered all seven standalone main figures, overview sheets covering all 24 supplementary figure numbers including their multi-page entries, selected supplementary figures at larger size, and selected rendered manuscript pages. This is not a complete page-by-page typography certification of all 85 pages.

Checks run on the current sources:

| Check | Finding |
|---|---|
| Journal word counter | 9,139; Summary exactly 150 |
| Legend checker | 33 legend entries pass; this includes multi-page supplementary entries |
| Supplement numbering | S1–S24 and Tables S1–S7 in citation order; no renumbering needed |
| Existing SVG print checker | All 21 listed main/supplement/protocol assets pass |
| Reference keys | 51 cited keys; none unresolved |
| Article bibliography | 51 entries; no duplicate DOIs detected |
| Combined DOCX navigation | 389 internal links; none dangling; no duplicate bookmark names |
| Author-supplied details | 14 known laboratory/deposition fields remain missing |

The reference check establishes structural integrity, not that every source supports every associated claim. Selected comparator and biological citations were checked against primary publisher/PubMed records below. No benchmarks were rerun; the earlier implementation tests are not being treated as new manuscript validation.

## Corrections to make before circulation

### 1. Two citation errors are visible to readers

In [crm_text.py](/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/docs/calibration/crm/manuscript/crm_text.py:98), `minimap2^{29}` points to Elias and Gygi's target–decoy paper in the current bibliography. The actual minimap2 entry is currently 48. Replace the number with `minimap2[@li2018]`; changing it to another hard-coded number would preserve the underlying problem.

The *Mkln1* Results paragraph at line 312 uses `@barbulescu2024`, the MAEA/FAM72A paper. Its Methods and Table S2 correctly use `@barbulescu2025`. The relevant publication is the [2025 MKLN1 study](https://academic.oup.com/jimmunol/article/214/10/2564/8239066). The proposal corrects both citations through keys, allowing numbering to regenerate.

### 2. The 95% recall sentence currently describes the wrong denominator

The planted-JH4 paragraph says “81% overall, and 95% of mutations were supported by at least three reads.” The analysis quantity is **95% recall among planted mutations with at least three supporting reads**. It is not the fraction of planted events assigned that support. The proposal makes both recall denominators explicit. This is a change of interpretation, without changing a numerical result.

### 3. Several correctly displayed numbers have incorrect source bindings

These are fragile source defects rather than currently incorrect rendered values:

| Source expression | Current display | Required correction |
|---|---|---|
| `Q>{mo_control_fold}` | Q>20, because the fold happens to be 20 | Use the quality criterion from the comparison, independently of the measured fold |
| `p = {mo_ugi_dead}.0021` | p=0.0021, because the inactive-Ugi call count is zero | Use the existing `{mo_arm_p}`, which computes and labels the paired test |
| Pim1 Geneious FDR uses `{sp_n114_gen_fdr}` | Both happen to round to 0.050 | Use `{sp_pim1_gen_fdr}` |
| “Four traces assigned to `{c2_ko_calls}` mice” | “two,” because the control mutation count is two | Bind to the actual assignment count, or describe shared assignments without borrowing a mutation-count key |

The mOrange p-value is **not presently numerically wrong**. It will become wrong if an unrelated count changes. The existing `{mo_arm_p}` currently yields “paired t, p = 0.0021.” The corresponding paired effect is 92 percentage points, whereas the pooled shares are 95% and 4%; these differ because they use different weighting, not because subtraction was done incorrectly. The unbounded t interval, 74–110 percentage points, should not be silently clipped; if a bounded interval is desired, change and report the statistical procedure.

### 4. Figure 2 points to the wrong panel twice

Both the annotation in [fig2_call.py](/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/docs/calibration/crm/figs/fig2_call.py:125) and the Figure 2 legend refer to **3H** for higher-fraction somatic-caller recall. That experiment is now **3C**. The existing legend checker passes because 3H exists; it cannot detect that the referenced panel addresses the wrong topic.

### 5. Figure S21 still contains obsolete counts and destinations

The rendered figure says “447 chromatograms from the six mice of Figure 5.” The source accounting is:

**447 archived traces → 184 assigned to the six mice → 2 excluded → 182 processed → 48 failed analysis → 134 retained.**

Figure S1's 184 input traces and the Results' 182 processed traces describe different stages and are compatible. S21's attribution of the entire archive to the six mice is not. Use the existing accounting keys and label the stage.

S21 also directs bookkeeping readers to S12, now the Immcantation comparison, and mutation-tree readers to S16, now SWHEL; mutation trees are S17. The bookkeeping validation is now in STAR Methods. This figure therefore needs a semantic cross-reference pass, not just file renaming. It also retains the pooled-category Fisher interval; either label its descriptive scope or remove it consistently with Figure 5B.

### 6. The Figure 5 legend describes an interval no longer drawn

Figure 5B displays r and slope, without a Fisher interval. Its legend still says “The displayed Fisher interval…”. Remove that sentence. The current figure builder explicitly explains why that interval was omitted. Retain the descriptive interpretation of pooled correlation, and keep the mouse-resolved panel C beside it.

### 7. Dataset counts and destinations need small corrections

- The initial real-library paragraph describes three external “loci”: Sμ in two cell contexts plus *Pim1*. These are **two loci in three assay contexts**, consistent with the twelve-locus counting convention elsewhere.
- STAR Methods still opens the engineered-diversification section with **“Four datasets”**, but now describes two retained datasets. Change this to two after the TAM/TCR removals.
- Table S2 sends the Karst UMI benchmark to Figure 7; its current destination is **S8D**. The simulation row's Figure 5 destination should also be reconsidered after removal of the bookkeeping panel; use Figure 3 and the relevant Methods/supplementary destination.
- In the control-number paragraph, make clear that the 2.8–3.1% false-discovery range belongs to **two to four controls**. Figure S5B gives approximately 0.1% for the one-control BH fallback. The current sentence can read as applying the range to all four designs.

### 8. The reviewer-package manifest is not current

[make_review_package.py](/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/docs/calibration/crm/make_review_package.py:22) still requires `FigureS7.png` and `FigureS12_bookkeeping.png`, neither present. It omits the current `FigureS7_indels.png` and `FigureS8_estimates.png`. The package would fail during copying, and its asset list does not yet represent the supplement. Deriving the list from the same figure manifest as the manuscript would reduce future drift.

The older HANDOVER also retains obsolete instructions about an S6 gap and contradictory word-count exclusions. The session handover and current code supersede those statements. Do not renumber the present supplement merely because the older note says to do so.

## Figure construction and visual priorities

The figures largely share a coherent palette and typography. Their main weakness is information density and annotation maintenance, with a small number of real construction defects. Passing the SVG checker does not establish that labels fit inside the canvas or remain large enough after embedding.

| Figure | Assessment | Recommended action |
|---|---|---|
| **1** | The read/Sanger fork communicates an essential distinction well. The page is crowded with complete explanatory sentences, and PacBio enters the input row without an equally explicit junction-analysis branch. | Keep the fork. Shorten box text to labels, retain explanation in the legend/S1, and state where PacBio junction characterization differs from calibrated substitution calling. Avoid presenting UMI consensus as an additional sequencing platform. |
| **2** | The two-candidate example is an effective explanation of why quality alone is insufficient. Panel E's long annotation distracts from its allele-fraction distribution. | Correct 3H→3C in both places. Shorten the gold note to the tested range and move detailed comparator recall values to Figure 3/legend. “Retained candidate” is safer than implying independent proof of a native biological mutation. |
| **3** | The two-row progression from calling to quantification is sensible. Panel D overlays dataset-shape and caller-color encodings, its legend occupies data space, and zero false calls are assigned an artificial x=0.05. | Use a separate no-false-call strip, a suitable linear inset or another explicitly separated representation. Reuse S6's treatment where appropriate. Move the shape legend outside the data and retain the n=3 permutation-resolution explanation. Update the legend in the same edit batch. |
| **4** | Highest visual priority. Off-target row labels in panel F are clipped at the left edge in the native PNG and in the newly rendered manuscript. The bottom row of library labels is also crowded against the canvas. Nine panels mix biological comparison, error calibration and read accounting. | Reserve a wider left margin for F, wrap or shorten its row names, and add bottom clearance. Keep the native physical size fixed rather than fixing clipping by arbitrary tight cropping. Consider moving H/I to existing Methods/supplementary workflow material to enlarge A–G; this would not require an eighth main display item. |
| **5** | The matched-mouse panel C is a material improvement over pooled correlation alone. The positional panel has a technical terminal-bin issue. | `np.arange(101, 655, 25)` produces a last midpoint at **638.5**, outside xlim 100–630; only positions 626–629 of that nominal bin fall within the 529-bp window. Explicitly show and label the short terminal bin, or use exposure-normalized bins. Recompute the plotted positional comparison if binning changes. Preserve mouse identities; subtle within-mouse connectors could improve C without changing its requested bars. |
| **6** | Whole-amplicon decline and clone weighting are clear strengths. The tree requires careful reading because node counts are nested and 21 reads are incompatible. | Keep “read-level reconstruction” in the title and legend; retain explicit nested-count semantics. Consider placing the mutation-class legend above the tree and reducing the long annotation over it. Show per-mouse uncertainty in the targeting-quintile panel if it can be added without clutter. |
| **7** | The examples have clear biological purposes, but all five panels shown are Illumina; UMI and PacBio evidence are supplementary. Panel E's few retained calls make proportions coarse. | Use an Illumina-specific main-figure title or explicitly point to the other-platform supplement. Keep paired lines in E and consider annotating each library's denominator. Report the paired experiment effect alongside the plot; do not treat “ns” as equivalence for *Ldha*. |

A further size issue affects the review document: matplotlib figures are designed at 174 mm but embedded at 165 mm, scaling a nominal 6-point label to approximately **5.69 points**. This does not mean the native submission asset fails its own checker. It means the embedded review copy is smaller, and any final-size check must use the actual intended publication width. Shrinking an already dense figure to fit a long legend is not an adequate readability fix.

### Supplementary material

- **S6:** retain the restored dilution result. Separate no-call strips and individual draws convey its detection cliff more faithfully than a fabricated logarithmic zero. The author's paragraph is preserved in the proposal.
- **S7:** the shaded one/two-read region correctly signals an algorithmic floor, distinct from measured recall. Class-specific denominators are useful; maintain the separation between insertion/deletion classes and supporting-read effects.
- **S8:** the combined power/consensus page is useful but covers distinct validation questions. Keep the real-rRNA benchmark explicitly technical and the SHM UMI simulation explicitly simulated. The implemented consensus restriction on indels remains relevant; do not generalize this to validated UMI indel calling.
- **S9–S13:** standardize matrix denominator/color-scale explanations and clarify zero-count matrices. Different scales for sparse controls can be valid, but must be explicit. S10's symlog construction actually accommodates zero; it should not be criticized as if it used Figure 3D's artificial zero.
- **S14–S16:** these provide important external comparison evidence. The compact published/reanalysis tables are valuable, but tiny footnotes carry substantial meaning. Preserve discrepancies and distinctions between tests when condensing their legends.
- **S18–S20:** keep selected-cluster identity validation separate from population-frequency/FDR validation; keep the PI3K treatment comparison descriptive; retain PacBio breakpoint ambiguity and PCR-recovery qualifications.
- **S21–S23:** three supplementary explanations now overlap with Figure 1. Decide which is the platform comparison, which the procedural detail and which the mathematical explanation. Resolve S21's stale content before considering consolidation or renumbering.
- **S24:** important evidence about conditional robustness. Case-only artifacts produced calls in every evaluated null library; control calibration cannot recognize an artifact absent from controls. This result should support a precise scope statement, not be hidden or converted into a general FDR guarantee.

## Language and scientific positioning

The expanded Discussion has a useful structure and stronger sense of purpose than a sequence of qualifications would provide. Preserve its six subjects: what the study establishes, its relationship to other tools, operating ranges, platform-specific value, reuse of published data, and practical contribution. Let the final sentence of each paragraph explain what the evidence enables.

Several expressions need narrowing within that positive structure:

| Current wording | Issue | Suggested direction |
|---|---|---|
| “the choice is no longer necessary” | Suggests the accuracy/throughput tradeoff is resolved generally | The evaluated workflow supports cloning-free SHM measurement with empirically assessed background |
| “whatever a platform delivers” | Stronger than supported input formats and preparation routes | Delivered reads or chromatograms in the supported workflows |
| “instead of trusting quality scores” | The method still uses quality information | Rather than relying on quality scores alone |
| “background is measured rather than modelled” | Background measurements fit a statistical model | A control-informed model, with distinct analysis routes |
| “every well-designed SHM experiment” includes AID-null controls | Overgeneralizes, including beyond this study's own public datasets | Where appropriate controls are available |
| ClairS/DeepSomatic variants occupy “a large share of reads” | Their intended scope includes low-frequency somatic variation | Explain the unusually low support in these assays and report the evaluated configurations |
| “unchanged” or “statistically indistinguishable” | A non-significant test does not establish equality | Similar observed means / no detected difference, with the relevant uncertainty |
| “without retuning” | Needs qualification after the Ramos-derived exclusion rule | Presets transferred, while one documented exclusion rule was developed during evaluation |

The [ClairS primary paper](https://www.nature.com/articles/s41592-026-03152-4) explicitly examines low allele fractions, including benchmarking at VAF≥0.05. The paper's strongest comparison is therefore the measured operating range on these SHM amplicons, rather than a universal characterization of what somatic callers can detect.

Use “biological replicate” for the general framework, reserving “animal” for animal experiments. N114's three subclones came from one transduction, as the table already states; inference across those subclones is not replication across independently performed transductions. Similarly, the Illumina deposit replicates with unspecified origin should not acquire stronger independence claims through generalized prose.

Keep the two error-handling operations distinct: the margin calibrates the per-candidate test; the decoy procedure selects a threshold. Keep the separate Sanger counting route and the PacBio junction module equally explicit. State the empirical FDR limitation once clearly in the main description and keep the detailed consequences in Limitations/Methods, leaving space for affirmative interpretation in the Discussion.

Minor prose corrections include lowercasing “Activation-induced” in the Introduction, adding a comma before its relative clause, removing the redundant “Because…therefore,” replacing the semicolon before “while curation” and choosing a consistent spelling convention. The proposal preserves the author's dilution paragraph and most original phrasing; it does not impose a wholesale British/US conversion.

The title's “measured false-discovery rate” is strongest for truth-known benchmarks and potentially broader than justified for native samples. An optional alternative is **“SHMtrace integrates control-calibrated mutation calling with reproducible profiling of AID-associated diversification.”** This is an editorial option, not a required title change.

## Word reduction with supporting detail preserved

| Counted component | Current | Proposal | Change |
|---|---:|---:|---:|
| Results, including subsection titles | 4,721 | 3,298 | −1,423 |
| Main figure titles and legends | 1,780 | 864 | −916 |
| Discussion | 1,145 | 1,148 | +3 |
| Introduction | 488 | 488 | 0 |
| Limitations | 355 | 355 | 0 |
| Resource availability | 259 | 259 | 0 |
| Other counted content | 391 | 391 | 0 |
| **Total** | **9,139** | **6,803** | **−2,336** |

Results and legend changes alone produce **6,800 words**. The remaining three-word increase gives the revised Discussion affirmative endings at essentially its original length. Most of the numerical compression is in duplicate caller totals, the detailed two-candidate example, calibration mechanics, the long SWHEL accounting narrative and individual PacBio comparison statistics. Their Methods destinations retain the evidence while leaving the Results to report the biological or analytical conclusion.

The proposal retains all seven Results subsection headings, the author's dilution paragraph, the full matched-mouse comparison, the distinct studies, null/weak-signal findings, repair-mutant results, and limitations on native lineage and junction inference. It does not alter analysis parameters, numerical outputs or default behavior.

The count uses the existing project function with placeholders filled from current outputs. It excludes references and STAR Methods under the project's working interpretation of the journal guidance. Publisher author-guidance pages were inaccessible during this review, so this is not a newly verified statement of the journal's current counting policy. The project counter also omits some fixed section-heading text. Keep the approximately 200-word margin and check the final submission file against the journal's requested scope after implementation and author details are completed.

## Positive Discussion closing

The manuscript should close on the demonstrated contribution, with future work immediately before it. The proposed closing is:

> SHMtrace brings control-calibrated mutation detection, whole-amplicon profiles and comparisons across biological replicates into one reproducible workflow. By preserving mutation linkage and documenting analytical choices, it makes long-read SHM data more informative and supports transparent reanalysis of existing sequence collections. This provides a practical foundation for studying how AID targeting and DNA repair shape antibody diversification.

This emphasizes a distinctive combination without asserting priority over all existing software, universal accuracy, or a formally guaranteed FDR.

## Reference metadata to finish

The major remaining problem is citation placement, not broken keys. Targeted primary-source checks support the following:

- The [MKLN1 article](https://academic.oup.com/jimmunol/article/214/10/2564/8239066) is the 2025 paper, volume 214, pages 2564–2574.
- The [ClairS article](https://www.nature.com/articles/s41592-026-03152-4) is correctly recorded as Nature Methods 23, 1507–1520 (2026).
- The [DeepSomatic publisher record](https://www.nature.com/articles/s41587-025-02839-x) gives volume 44, 1569–1578 (2026), despite online publication in October 2025. Its 2026 issue citation need not be “corrected” to 2025 just because the DOI contains 025.
- The [JoVE record](https://www.jove.com/t/61551/analysis-somatic-hypermutation-jh4-intron-germinal-center-b-cells) supplies issue 170 and article e61551, absent from the current formatted entry, which has an empty page/article slot.
- The [GigaScience record](https://academic.oup.com/gigascience/issue/10/2) supplies article **giab008** for SAMtools/BCFtools; its current entry has an empty article slot.
- The [Seija PubMed record](https://pubmed.ncbi.nlm.nih.gov/42816610/) confirms the 2026 article and DOI, but also records middle initials missing from several manually entered authors and links an October 2026 author correction. Check that correction before freezing the citation; its contents were not retrievable here, so no implication for the underlying biological results is asserted.

The mOrange and Clair3 citation-key repairs documented in the handover are retained in current sources. A final reference pass should normalize missing article identifiers and punctuation while preserving verified titles/DOIs, then regenerate numbering after the minimap2 and MKLN1 fixes.

## Recommended implementation order

1. Correct citation keys, accidental placeholder bindings, recall-denominator wording and stale dataset/figure destinations.
2. Apply the counted Results/legend proposal in Python sources, integrating retained detail into existing STAR Methods subsections. Preserve the Summary's 150-word ceiling and all numeric placeholders.
3. Adopt the Discussion revisions selectively, retaining six paragraphs and the positive final synthesis.
4. Repair Figure 4 clipping, Figure 2's cross-reference, Figure 3D's zero display, Figure 5D's terminal bin and S21's obsolete content. Regenerate affected assets and verify the filenames actually embedded in Word.
5. Rebuild into a new review version, run structural/print checks and inspect actual-size rendered figure/legend pages. Refresh the package manifest before generating a reviewer ZIP.

The 14 missing laboratory/deposition details and independent-preparation validation remain distinct submission/scientific work. They were not filled by inference, and this editorial review does not establish cross-run generalization.
