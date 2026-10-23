# SHMtrace v1.3: figure critique and revision handoff for Claude

## Task

Revise the manuscript figures using the critique below. Improve scientific clarity, visual hierarchy, consistency, and readability at publication size. These are proposed revisions, not changes already applied. Use judgment: explain any recommendation that you do not implement rather than mechanically following every suggestion.

Preserve the original manuscript and figure assets. Make the revised figures and manuscript a separate, clearly named version. Do not change underlying analysis results, counting rules, filtering defaults, or statistical tests merely to improve appearance. Any proposed analytical change must be identified separately and validated before replacing a reported result.

## Files and working context

- Repository: `/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2`
- Original session handover: `docs/calibration/crm/SESSION_HANDOVER_2026-10-21.md`
- Build instructions: `docs/calibration/crm/HANDOVER.md`
- Codex manuscript: `output/reviews/shmtrace_v1_3_2026-10-21/codex_version/SHMtrace_CellReportsMethods_complete_v1.3_Codex_1.docx`
- Figure assets inspected for this critique: `output/reviews/shmtrace_v1_3_2026-10-21/codex_version/figures/`
- Codex preparation and authoring scripts: `output/reviews/shmtrace_v1_3_2026-10-21/prepare_codex_version.py` and `author_codex_docx.py`
- Original manuscript sources and build: `docs/calibration/crm/manuscript/` and `docs/calibration/crm/build/`

Paths above are repository-relative unless absolute. Inspect the current source and build route before editing. The manuscript is generated from Python; edits made only in Word can disappear on rebuild. The Codex copy is a separate revision: do not assume rebuilding the canonical manuscript reproduces its text. Likewise, do not assume the copied files in `codex_version/sources/` contain every applied prose revision; inspect the preparation script and `document_operations.json`.

## Overall assessment

Figures 3 and 4 need the most restructuring; Figures 5 and 6 have stronger foundations. Most data figures look like conventional scientific plots. Their weaknesses concern hierarchy, crowded annotations, and how clearly the plots support the scientific claim, rather than evidence of AI authorship.

The central design principle is **one primary scientific question per figure, with the strongest evidence receiving the most space**. Aesthetic traits cannot reliably identify whether AI created a figure.

## Figure 1 and workflow illustrations: critique carried forward

The workflow design has overly uniform containers, excessive information density, mixed illustration styles, small labels, and a generic dashboard appearance. The deeper issue is that biological motivation, sequencing errors, calibration, software processing, and interpretation receive similar visual weight.

The figure should make one idea immediately clear: **matched controls inform whether an observed mutation exceeds the technical background**.

Recommended changes:

1. Give control calibration and mutation calling more space than routine file inputs and preprocessing.
2. Specify the comparator when contrasting with quality-based approaches. Do not imply that all existing methods ignore controls or error models.
3. Replace decorative miniature charts with one readable worked example of sample evidence versus background. Mark synthetic examples as schematic.
4. Make arrows distinguish data flow from statistical comparison and interpretation.
5. Distinguish observed mismatches, candidate variants, and called mutations.
6. Describe decoy calibration as empirical; do not imply a formally established general FDR guarantee.
7. Use a small number of stable colour meanings; remove redundant numbered badges and most enclosing cards.
8. Check labels at final publication width.

Suggested structure: **sample and control inputs → control-informed calibration → mutation calls and downstream summaries**.

## Main figures

### Figure 2 — Why controls distinguish mutations from background

The underlying argument is good: two observed changes can have similar quality scores but different biological credibility because their control backgrounds differ.

- **The strongest comparison is visually fragmented.** Panels A–C require readers to move between allele fractions, counts, quality scores, and significance annotations. Keep A and B adjacent with identical formatting, and place a simplified quality comparison directly below them. Explicitly connect the message: similar quality, different control background, different calling outcome.
- **Panel D is difficult to decode quickly.** Several curves, logarithmic axes, and threshold annotations compete. Directly label candidate and decoy curves, emphasize the selected threshold, and move the mathematical explanation into the legend. Describe this as empirical calibration rather than implying a general FDR guarantee.
- **Panel E mixes distributions with benchmarking commentary.** The yellow shaded region and comparator-recall text interrupt the histogram. Move comparator information to Figure 3. Retain a clean distribution with a clear denominator and only necessary annotations.

Suggested structure: a large worked comparison, followed by two smaller panels showing threshold selection and allele fractions of retained calls.

### Figure 3 — Performance and simulations

This is the most overloaded main figure. Eight panels cover caller comparisons, sensitivity, false calls, lineage reconstruction, quantification, and group testing.

- **Separate calling performance from downstream inference.** Keep spike-in benchmarks and the recall–false-call trade-off together. Move clone estimation, quantification, and group-testing simulations into a dedicated supplementary figure, retaining a main-text summary only if essential.
- **Simplify method encoding.** Multiple colours, markers, and line styles must be remembered across panels, while the somatic-caller panel introduces another set. Use one persistent SHMtrace colour, subdued comparator colours, and direct labels where possible.
- **Give Panel D more prominence.** It combines sensitivity and false calls and communicates the central trade-off better than several separate recall curves. Enlarge it. Explain that the false-call axis is linear near zero and logarithmic beyond the transition.
- **Separate calibration and power in Panel H.** The “none” conditions assess false positives; the 2× and 4× conditions assess power. Visually distinguish and label their different purposes.
- **Present uncertainty consistently.** Where repeated simulations were performed, identify what error bars represent and how many runs contributed. Do not invent uncertainty bars where there were no independent repetitions.

Suggested structure: four panels centred on sensitivity, false calls, and their trade-off; detailed downstream simulations in the supplement.

### Figure 4 — Biological separation and specificity

The breadth of evidence is valuable, but nine panels make this feel like a collection of validation results rather than one coherent figure.

- **Reduce repeated diagonal labels in A–E.** Use aligned small multiples with one shared method-label row and explicit sample sizes. If y-axis scales differ, make that conspicuous.
- **Clarify zero calls.** Triangles below logarithmic axes can look like tiny positive values. Use a clearly separated strip labelled “0 calls,” consistently throughout the paper.
- **Explain the heatmaps.** State exactly which summaries form the case/control ratio and how zeros are handled. Give both heatmaps explicit colour scales; colour alone should not imply a statistical judgment.
- **Give Panel H enough space.** The margin-versus-held-out-control result communicates an important calibration limitation. Distinguish the observed held-out rate from the reference target.
- **Move Panel I to the supplement.** Read attrition is important QC, but interrupts a figure about biological separation. Retain a concise usable-read summary in the main text.

Suggested structure: biological contrasts above; separation, control calls, and held-out calibration below. Move detailed read accounting out.

### Figure 5 — Sanger–Nanopore comparison

This is one of the strongest figures because it addresses a concrete validation question. It needs clarification more than wholesale redesign.

- **Distinguish correlation from agreement.** Spectrum and positional correlations show similar patterns; they do not establish interchangeable absolute measurements. Name the different counting units and denominators explicitly in Panel A.
- **Foreground the three mice in Panel C.** Enlarge mouse markers or use paired points within each mutation class, with mean and variability shown more lightly.
- **Reduce crossings in Panel D.** Use thinner lines and smaller markers, or aligned profiles with identical axes. Retain the terminal-bin annotation, but explain bin normalization in the legend rather than crowding the title.
- **Label Panel E reference lines.** Readers should immediately understand whether they represent motif availability or another baseline.

Suggested structure: retain the overall organization, strengthen replicate visibility, and make agreement in pattern the explicit claim.

### Figure 6 — Full-length information and clone structure

This contains distinctive material, but the connection between its parts is not immediately obvious.

- **Reduce competing profiles in Panel A.** Use lighter individual mouse lines, a clearly dominant summary, and direct labels for the shared Sanger window and additional Nanopore interval.
- **Explain Panel B normalization.** Define the targeting index plainly, state how quintiles were constructed, and identify the underlying observations. A five-point trend must not imply more independent replication than exists.
- **Simplify Panel C.** Clone proportions, junction recovery, unassigned reads, and read- versus clone-weighted frequencies are packed together. Put the weighting comparison in a small adjacent plot and leave stacked bars to communicate clone composition.
- **Add a reading guide to Panel D.** Explain mutation ticks, node sizes, numbers, and branching. Distinguish reconstructed relationships from directly observed sequences, and explain reads outside the displayed reconstruction.

Suggested structure: two explicit themes—information beyond the Sanger window, then clone-aware interpretation—with a clear visual transition.

### Figure 7 — External datasets and generality

This demonstrates breadth, but currently resembles five independent mini-figures.

- **Organize different endpoints.** A–B address mutation burden, C spectrum, D indels, and E repair-dependent changes. Group panels by the question being tested, with short descriptive titles.
- **Stabilize colour meanings.** Preserve biological meanings where possible; otherwise provide local legends rather than relying on familiarity.
- **Keep small counts prominent in Panel C.** Present spectra descriptively and avoid allowing percentages to suggest strong precision.
- **Emphasize biological observations in Panel D.** Enlarge replicate points and reduce bar weight. Clarify whether “events per clone” is comparable across conditions.
- **Clarify trajectories in Panel E.** Use small multiples or directly labelled trajectories, and explain precisely what each line connects.

Suggested structure: external validation organized around burden, mutational signature, and indels/repair, with dataset provenance visible.

## Supplementary figures

These are primarily presentation critiques. Implementing them need not change underlying results. Where multiple image assets share one supplementary number, review them together.

| Figure | Main concern | Suggested fix |
|---|---|---|
| **S1: experimental workflow** | Experimental design and eight processing stages compete for attention. | Separate experimental design from computational processing. Enlarge amplicon coordinates and move detailed thresholds into the legend. |
| **S2: quality versus observed error** | The high-quality tail can appear as reliable as well-supported bins. | Keep the control-base-count panel aligned; identify low-support bins and define uncertainty shading. |
| **S3: software outputs** | Six cards contain excessive instructional text and miniature plotting. | Show three representative outputs at readable size. Move filenames and “watch for” guidance to documentation or a table. |
| **S4: calibration sweep** | The transformed percentage axis and large shaded area require effort to interpret. | Explain the scale; emphasize the equality line and chosen operating point. Keep false calls and recall as aligned companion panels. |
| **S5: number and quality of controls** | “Clean/off/on” labels are compressed; the one-control condition changes procedure. | Spell out conditions. Visually separate the one-control fallback from decoy-calibrated conditions. |
| **S6: dilution/titration** | Four plots mix recovery, detection, and mutation counts. | Group expected-versus-observed recovery above detection/count behaviour. Explain the zero-call strip and distinguish reference lines from fitted relationships. |
| **S7: indels** | The excluded low-support region can be confused with measured poor performance. | Mark excluded conditions as outside the calling rule. Label numbers beside bars as planted-event counts. |
| **S8: estimation and UMI analysis** | Single-read estimation, power calculations, and UMI performance form three stories. | Split into estimation/power and UMI figures. Make molecule counts and uncertainty definitions prominent. |
| **S9: N114 and substitution matrices** | Dense matrices are misleading when a control has only one or two calls. | Keep absolute counts primary. Mark zero-event spectra as undefined; warn against interpreting tiny-total percentages. |
| **S10: depth sensitivity** | Numerous trajectories obscure the practical depth conclusion. | Keep sample identities consistent and emphasize representative depths. Define recovery relative to each library’s reference depth. |
| **S11: positional background and read maps** | Overlaid traces and small insets obscure positional differences. | Separate control profiles into aligned tracks; share coordinates with read maps and explain read ordering. |
| **S12: Immcantation comparison** | Ratio bars can amplify differences caused by near-zero controls. | Lead with absolute values. Treat ratios as secondary and state differences in reference, filtering, and counting rules. |
| **S13: Sanger matrices and readable window** | Readability and spectrum validation are different questions; matrix labels consume space. | Separate the window assessment from matrices. Align coverage and mutation profiles, with an explicit denominator for edge bins. |
| **S14: second Sanger cohort** | Agreement plots, distributions, and a detailed table compete. | Lead with per-mouse agreement. Move the table out of the figure and distinguish published from reanalysed sample counts. |
| **S15: Mkln1 reanalysis** | Spectrum agreement and test reproduction are crowded; low-count controls are awkward on the burden axis. | Separate biological results from the numerical table. Give controls readable space and explain differences in included clones. |
| **S16: SWHEL reanalysis** | Rounded correlations of 1.0000 and a “Claims” table look overly definitive. | Use less precision; show absolute differences or recovery counts. Rename the table “Published and reanalysed effects” and state test differences. |
| **S17: trees and full-length spectra** | Trees have tiny labels; overlapping-window correlations can be overinterpreted. | Give trees dedicated pages with a readable key. State that windows overlap; consider showing the additional interval separately. |
| **S18: engineered systems** | Different validation questions are squeezed into two panels. | Give each dataset a short local explanation. Show recovered, missed, and additional calls alongside external-list agreement. |
| **S19: public Sμ data** | Essential coverage and read-layout differences are buried in annotations. | Enlarge the coverage schematic. Keep shared-window and repeat-window results separate, with usable sample counts beside each contrast. |
| **S20: PacBio junctions** | Yield, pooled composition, per-sample phenotype, and ambiguity use different units. | Label the observational unit in every panel. Show sample-level composition; separate ambiguity from biological microhomology. |
| **S21: platform comparison** | Dense prose makes the figure resemble a miniature review article. | Replace lower cards with a compact table: measurement unit, span, calling rule, output, and limitation. Keep one diagram per platform. |
| **S22: overview** | Repeats the workflow with many boxes, icons, and small explanations. | Retain only if it serves a distinct purpose; otherwise consolidate its strongest elements into Figure 1 or S1. |
| **S23: statistical engine** | Six formula cards and miniature plots are difficult to read. | Focus on background estimation, candidate testing, and empirical threshold selection. Move routine preprocessing elsewhere. |
| **S24: artifact stress tests** | “Any discovery” and per-library false-discovery proportion can be conflated. | Define both outcomes. Explain intervals and foreground the case-only artifact result as a boundary of control-based calibration. |

## Implementation priorities

1. Restructure Figures 3 and 4 around fewer, clearer questions.
2. Consolidate overlapping workflow explanations in Figure 1, S1, S3, S21, S22, and S23.
3. Standardize colours, zero-call symbols, counting units, typography, and replicate displays.
4. Refine Figures 5–7, especially agreement claims and clone interpretation.
5. Check every figure at its intended printed width.

## Validation and deliverables

- Preserve original assets and manuscript; produce a separately named revision with editable figure sources and final exports.
- Update panel letters, legends, in-text references, and supplementary numbering together whenever panels move or figures split/merge. Run existing consistency checks; do not bypass assertions.
- Maintain all scientific detail in an appropriate destination when simplifying a figure. Record moved panels and explanations in a relocation list.
- Use explicit observational units: reads, molecules, unique variants, clones, libraries, animals, or simulation runs. Do not substitute one for another.
- Do not turn correlations into claims of absolute agreement, pooled events into biological replicates, or empirical decoy thresholds into a formal FDR guarantee.
- Check zeros, logarithmic axes, terminal-bin widths, denominators, legends, sample sizes, and uncertainty definitions against actual source data.
- Render and inspect the revised manuscript, not only standalone images. Check minimum type size, contrast, clipping, overlaps, and readability at final width using the existing project print checks.
- Provide the revised manuscript, editable figure sources, publication exports, a concise change log, and any unresolved analytical or author decisions.
