# Independent review of Codex's SHMtrace work (manuscript Codex 1 + figures v2)

Claude · 22 October 2026 · read-only review; no file outside `claude_figures_v1/` was changed.

**What I inspected:**
- the final Codex DOCX operations (`codex_version/document_operations.json`) against the canonical filled text (`filled_blocks.json`);
- all 43 pages of `codex_figures_v2` and both stages of its build (`build_data_figures.py`, `custom_figures.py`);
- the saved analysis outputs behind specific claims (`crm/out/`, `docs/calibration/figures/`);
- the generating scripts where a claim depended on them (`mutation_tree.py`, `mutation_tree_figure.py`, `spikein.py`).

I also ran a parallel figure revision (`claude_figures_v1/`), so where a problem also affects my version, I say so.

**Bottom line.** Codex's numerical preservation is real: the hash checks hold, and none of the reported numbers I compared has changed. Most figure redesigns are improvements. However:
- one figure change is a **data-selection regression** (S17);
- the workflow supplements **lost substantive content**;
- the prose compression **dropped three qualifications that matter**.

Separately, several overstatements already in v1.3 survived both revisions. **None of these invalidates a reported result.** The largest scientific gap is the one the Limitations already state: there is no independent preparation or run.

---

## 1. High-priority scientific errors or unsupported interpretations

| # | Issue | Evidence | Origin | Fix |
|---|---|---|---|---|
| H1 | **S17 trees: wrong selection, and an incorrect correction.** The original S17C is "Wt4, clone 9 · 16 reads · 100% of mutated reads fit". Codex replaced it with WT3 clone 14, and swapped WT3 clone 0 for the same clone already shown in Fig 6D. Its manifest says this "preserves the actual third example rather than incorrectly describing it as WT4". The README says not to "introduce a WT4 reconstruction absent from the saved selected-tree data". Both statements are false. | `crm/figures/FigureS17_mutation_trees.png` is a symlink to `figures/mutation_trees/mutation_trees_by_mouse.png`, built by `mutation_tree_figure.py:233` ("One tree per Wt mouse: its most branched clone"). Codex instead read `trees.json["trees"]` (WT2:0, WT3:3, WT3:14), which is the separate *contrast* figure set (`mutation_tree_figure.py:286`). | Codex v2 | Revert to per-mouse selection (WT2:0, WT3:0, WT4:9). The original legend ("most branched clone of one Wt mouse") was correct. |
| H2 | **The Sanger/nanopore "different units" sentence misdescribes the Sanger measure.** Results ¶46 (v1.3 and Codex) says Sanger "counts mutations in individual cloned products, whereas the nanopore analysis counts distinct variants". The Fig 5A Sanger values are distinct variants per 10^5 clone bases. | Wt2: 135 unique / 18,964 bases = 711.9 per 10^5, exactly `stats.json density_sanger.WT2`. Counting every occurrence gives 885.9. The Wt3 figures are 1,121.6 and 1,892.0 respectively. | v1.3, kept by Codex. **My v1 Fig 5A label ("Q≥20 substitutions in clones") repeats it.** | Both are unique-variant densities. They differ in the calling rule (Q≥20 counting with background reported, versus control-calibrated calls) and in sampling: 26–40 single molecules versus 500 reads with clonal redundancy, which changes how many distinct variants are discovered per base. Reword the Results, Fig 5 legend and Fig 5A labels. |
| H3 | **"The same AID-null contrast"** (Discussion ¶69). Sanger gives 38-fold separation; nanopore gives 141-fold, and that figure depends on the one-event floor for zero-call controls. | `stats.json`; Fig 5A | v1.3, kept | Say "the same direction of genotype contrast" or "both separated the genotypes". |
| H4 | **Comparator cherry-pick in the Discussion.** Discussion ¶67 says counting every mismatch "separated AID from empty-vector N114 libraries by 2-fold where SHMtrace separated them by 188-fold". The 2-fold uses the wrong (IMGT) germline. With the parental germline and Q20 masking the same route gives 6- and 32-fold. SHMtrace's 188-fold depends on one-event floors for zero-call empty-vector libraries. | Results ¶40; S12 legend | v1.3, kept | Quote the 6–32-fold range with the parental germline, and say that 188-fold is floor-dependent. |
| H5 | **Lost qualification: the same controls choose the margin and the decoys.** The v1.3 sentence ("This is an empirical estimate targeting 5%, not a guarantee: it relies on a few controls … and the same controls also choose the margin") has no counterpart anywhere in the Codex text, including STAR Methods. This is the main reason the decoy estimate is conditional. | Phrase search over `codex_text` | Codex compression | Restore one sentence in Results ¶21 or STAR Methods. |
| H6 | **Lost qualification: "ns" with n = 3.** v1.3 said "With three libraries per group these tests are weak, so a caller marked ns has not been shown to fail". This is gone from the Codex text. Without it, Fig 4's ns marks read as evidence against comparators. | Phrase search | Codex compression | Restore it in the Fig 4 legend. It costs about 20 words. |
| H7 | **The power requirement is understated.** The Discussion says a two-fold difference "would require approximately twenty" animals. The S8 legend says "about 24". In `power_s3.json` the power is 0.73 at 20 animals and 0.91 at 30, so 80% power falls near 24. The pilot CV's 95% interval is 0.44 to an unbounded value, so the uncertainty is larger than "considerable" implies. | `out/power_s3.json` | v1.3, kept | Write "about 24 (and many more if the CV is underestimated)", consistent with the S8 legend. |
| H8 | **FDR and FDP used interchangeably.** The Summary has "measured false-discovery rates of 0.050 and 0.045"; the Results have "observed FDR 0.050" and "FDR 0.083"; the Discussion has "observed FDR 0.18". These are pooled realised proportions over nine libraries. Codex's own Fig 3 correctly relabels them "Observed FDP". | Results ¶26–28; Summary | v1.3; Codex fixed the figure only | Use "observed false-discovery proportion (FDP)" for realised values. Reserve "FDR" for the target. |
| H9 | **Fig 7A/B bars don't follow their legend.** The bars are geometric means with fixed display floors (0.01 and 0.005), while the legend says "one-event zero floors". For J_H4 AID-null this is 0.017 against 0.007 per 10^3 bp. | `codex_figures_v2/sources/revised/fig7_platforms.py:40,54,65` | v1.3, kept by Codex | Fixed in my v1. Codex v2 still has it. |

Lower-priority inconsistencies that should still be corrected:
- **"r = 1.0000"** is still in Discussion ¶70 (SW_HEL). Codex removed it from S16; my v1 removed it from the Results and S16 but missed the Discussion.
- **S10 wording.** Results ¶39 says J_H4 and N114 reached 90% and 59% "of discoveries at the largest tested depth". The S10 legend defines the denominator as the maximum over every depth and draw, so one of the two is wrong.

---

## 2. Manuscript editing: regressions, missing qualifications, references

**Checks that passed:**
- Placeholder values: `prepare_codex_version.py` asserts that the filled numbers equal the reviewed snapshot.
- References: 51 cited keys, 0 unknown, 0 broken internal links (`reference_and_package_checks.json`).
- The reference metadata fixes are correct: `sible2021`, `danecek2021`, the `seija2026` initials, the Pérez-Durán and Clair3 citations.
- Moving material into STAR Methods mostly worked. The "data-dependent candidate selection" rationale, the one/two/three-control behaviour and the benchmark-vs-`run_project` distinction are all present.

**Regressions and gaps:**
1. H5 and H6 above. A third qualification is also gone: v1.3 said the margin factors "do not imply that sequencing errors are universally five or ten times higher" and are separate from the decoy threshold. Without it, a reader may take the margin as a measured error inflation.
2. The pre-adjustment calibration evidence ("5.4% of J_H4 and 9.6% of N114 held-out tests fell below this threshold" without the margin) is gone from the text. It is still visible in Fig 4H, so this is acceptable if the legend points to it.
3. The Discussion endings are mostly proportionate. Two exceptions: "Agreement between platforms is itself evidence" (soften it, see H3), and "Both replacements held" (fine, but conditional on the one-run design).
4. **The Word manuscript was not rebuilt with the v2 figures.** The Codex 1 DOCX embeds the pre-v2 images. The v2 legends (`revised_legends.json`) contain stale or incorrect references:
   - S25: "Former Figure S25A–G" and "S25D–E" should both refer to Figure 3E–H.
   - S26: "Former Figure S26A" should be Figure 4I.
   - S4: says both "a pooled proportion of zero is drawn on the 0 row" and "zero is retained … on a symmetric-log axis".
   - Fig 4: still says "triangles, no calls", alongside the new strip.
   - S8: the full four-panel legend is duplicated on both pages.
   - S17 trees: titled "Mutation trees in each mouse" while stating they are not one per mouse.
5. The word count (6,823 by `journal_word_count`) covers the Introduction through Limitations plus front matter, the AI declaration, the SI index and **main** legends. Integrating the v2 legends changes it. My v1 legend set takes it to about 6,960.

---

## 3. Figure and design regressions, export and reproducibility problems

**Content regressions:**
- **S17**: see H1.
- **S1, S21, S22 lost explanatory content.** The v2 S1 shows only experimental design and says "Computational processing is shown in Figure 1". But the v2 Figure 1 no longer shows read preparation (length and expected-error filters, competitive alignment, the paralog-allele filter, depth matching) or the Sanger route steps. The processing pipeline now appears in no figure. The v2 S21 and S22 are mostly empty pages carrying a few lines of text. Either:
  - restore a compact processing panel (S1 B–D in the original) and drop S22, or
  - merge S21 and S22 into one table.

**Main-figure presentation problems:**
- **Figure 1 is too sparse to carry the paper's central idea.** It drops the worked contrast with a quality-only decision (original panel A). Its rule, "Candidate absent from controls → evidence can exceed background", reads as if absence from controls decides a call, which the legend then has to deny.
- **Zero strips (Fig 4A–E, Fig 7A–B)** are drawn on a continuous log axis, with real decade labels inside or next to the strip (Fig 7A's 10^-2 tick lies inside it). The "0 calls" text overprints the triangles. This is the defect the session handover warns about ("zero is not a point on a logarithmic axis"). Use a separate axis segment with a break, as in S6, S10 and my v1 Figs 4 and 7.
- **Other main-figure fixes:**
  - Fig 4: the method key takes a whole panel slot, and H's "1% target" should read "1% expected if calibrated".
  - Fig 3D: the locus legend covers the Pim1 LoFreq-permissive point.
  - Fig 5E: the legend overlaps the y-axis tick.
  - Fig 6C: the legend markers overlap the Wt3 point.
  - Fig 2E: the 50% line is still labelled "heterozygous germline". The rule is "≥50% excluded".

**Reproducibility:**
- `build_data_figures.py` lays figures out by intercepting `save` and repositioning axes by **creation index**: "Creation order in the preserved source is A,B,D,E,F,G,H,C", `keep=[0,1,7,2]`.
- `clean()` deletes any text matching `[A-I]` exactly.
- `zero_strip()` finds triangles by counting path vertices.
- `size()` raises every font below 6.5 pt to 6.5 pt after layout.

Any edit to a builder can silently scramble panels. This is acceptable for a review pack but not as the maintained source. Port adopted layouts into the builders themselves.

**Export consistency:** I found no SVG/PDF/PNG mismatches. TIFFs are outside the ZIP, as documented.

**Things Codex got right that are worth keeping:**
- Fig 7E1/E2 split: I verified its experiment pairing against `fig7_morange.json`.
- Fig 5C: paired mice.
- Fig 6C: weighting separated from composition.
- S9: the count matrices and the pMX-Empty grouping.
- S20: per-sample composition.
- S24: the definitions.
- S16: precision removed.
- S8: page split.
- S25 D/E: false-positive results split from power.

---

## 4. Claim–evidence table (central methodological claims)

| Claim | Supporting output | Kind of evidence | Supported as worded? |
|---|---|---|---|
| Control calibration separates AID mutations from recurrent background where quality cannot | Fig 2, S2, S23 | One worked example plus per-locus error curves | Yes, as illustration. The general case rests on Fig 3/4. |
| The decoy threshold delivers about 5% FDP at the default | Fig 3A–B, S4 | Planted substitutions in real control reads, 9 libraries/locus | Yes, **pooled**. Single libraries reach 10.6%, and targets below 5% are not delivered (both stated). It is conditional on representative controls (H5). |
| The calibration survives realistic nulls and noisy controls | S5, S24 | Resampling three existing null libraries from one run | Conditionally. The case-only artifact is a stated boundary. There is no independent preparation. |
| Recall at low support exceeds quality-based callers at comparable false calls | Fig 3, Table S1 | Planted (simulation-on-real-reads) | Yes for the configurations tested. Geneious recalls more single-read events at more false calls (stated). |
| Lower background and better separation than comparators on real data | Fig 4 | Direct, n = 3 per group, harness calibration | Yes descriptively. Folds are floor-dependent (stated), ns is not evidence of failure (H6), and comparators were not given a control-based filter (see §6). |
| Calls behave like real AID mutations (diluted real SHM) | S6 | Mixtures of real libraries | Yes. Detection is stepwise, and the draws are coarse (stated). |
| Nanopore reproduces Sanger patterns in the same mice | Fig 5, S13 | Descriptive correlations, 3 Wt mice | Pattern similarity, yes. Equivalence, no (stated). H2/H3 wording needs fixing. |
| Long reads add whole-amplicon profiles and read-level linkage | Fig 6, S17 | Descriptive | Yes. The trees are read reconstructions (stated). S17 must be restored (H1). |
| Animal-level testing is required; pooled reads inflate false positives | S8E (old 3H) | Simulation, 300 per cell | Yes, strongly. |
| Transfers to public datasets and platforms | Fig 7, S18–S20, S12 | External biological expectations | Largely. The R174E/R178E artifact and the Ramos noisy control are disclosed. The exclusion rule was developed post hoc (disclosed in the Limitations). |
| Sanger reanalysis reproduces published results | S14–S16, Table S4 | Counting fidelity against published tables | Yes, for counting and spectra. Two of the three studies come from the authors' laboratory (stated). |

---

## 5. Decision per main figure (Codex v2)

| Figure | Decision | Concrete fixes |
|---|---|---|
| 1 | **Revise further** | Restore a worked quality-only vs control-calibrated contrast and a compact processing strip. Keep the platform-route row. Replace "absent from controls → can exceed background" with a support-vs-background comparison. |
| 2 | **Adopt, minor fixes** | Label the 50% line "≥50% excluded (germline filter)". Optional: direct labels in D. |
| 3 | **Adopt, minor fixes** | Move D's locus legend off the data. Put the linear/log note inside the axis. Fix the SHMtrace label overlap in C. |
| 4 | **Revise** | Put zero calls on a separate axis segment. Remove the key panel (use labels, or one shared method column). Rename the H reference line. Restore the n = 3 caveat in the legend. |
| 5 | **Adopt, fixes** | Correct the unit wording in A (H2). Move E's legend. |
| 6 | **Adopt, minor fixes** | Fix the legend overlap in C. Reword "21 do not" as "21 carry mutations of more than one branch". Node counts include them (`mutation_tree.build`). |
| 7 | **Adopt E1/E2; revise A/B** | One-event geometric means (H9). Proper zero segment. |
| S17 | **Revert selection** | WT2:0, WT3:0, WT4:9 (H1). |
| S1/S21/S22 | **Revise or merge** | Restore the processing content. Consolidate the two text-only pages. |

My v1 already contains fixes for Fig 4's zero segment, Fig 7A/B, Fig 6D's wording and S7's shading. It has its own errors, H2 and the Discussion "r = 1.0000", so neither version should be adopted wholesale.

---

## 6. Remaining validation: essential versus optional

The yardstick is a bounded methods paper that is one of several PhD papers.

**Essential (no new wet lab):**
1. The corrections above (H1–H9) and the restored qualifications.
2. **Report, or explain the absence of, the ablations that already exist.** `spikein.py` ran `v0.1` (per-position null) and `pbinom` (quality model only) on J_H4. Both give 55.4% recall and 2.78 false calls per library, against 80.6% and 2.56 for SHMtrace. This directly answers "what does the control-calibrated test add?". **The two modes give identical numbers in all nine libraries.** Before quoting them, check that they are not the same code path.
3. Integrate one figure set into the DOCX, renumber the supplements by first citation, and rerun the word count.

**Strongly recommended (computational, low cost):**
- A naive control baseline: the Geneious rule minus any site seen in a control. This separates "using controls" from "SHMtrace's statistical calibration", which is the comparison a reviewer will ask for. It can probably be computed from the existing per-library call sets.

**Optional, or future work (already acknowledged as limitations):**
- An independently prepared run with new mice, the single most valuable addition.
- UMI-tagged SHM libraries.
- More animals per group.
- Biological nanopore indels.

None of these is required to support the claims **as currently qualified**.

---

## 7. Manuscript integration checklist

- [ ] Choose one figure set per figure (§5) and rebuild the DOCX from the operations route (`claude_figures_v1/manuscript/` or Codex's).
- [ ] Panel moves:
  - Fig 3E–H → S8 page 1 (my v1) or S25 (Codex); choose one.
  - Fig 4I → S1E or S26.
  - Update the Results ¶33–36 and ¶44, STAR Methods (high-allele-fraction rule, group tests ×3), Table S6 (Karst → S8D/H) and the S1 legend.
- [ ] Renumber S25/S26 by first citation. S25's content is first cited at Results ¶33, after S7, so it lands at S8 and shifts S8–S24. Check with `renumber_supp.py`, minding its two known traps.
- [ ] Replace every legend for Figs 1, S1, S3, S11, S21–S23 (new content). Fix the stale "Former …" references and the duplicated S8 legend.
- [ ] Apply H2–H8 and the S10 wording in the Results/Discussion. Restore H5 and H6.
- [ ] Use "FDP" for realised values throughout the text, legends and figure labels.
- [ ] Rerun `journal_word_count` (Codex 6,823; my estimate with v1 legends about 6,960, against a 7,000 cap).
- [ ] Render to PDF. Inspect every figure page at 174 mm and run `print_check` on the final assets.
- [ ] Hash-check the originals again after integration.

---

## 8. Calibrated overall assessment

**Established:**
- On these assays, a control-calibrated test with an empirical decoy threshold recovers planted mutations at low read support with pooled FDP near the 5% target.
- It produces fewer calls in AID-null controls than the quality-based configurations tested.
- Its calls behave like AID mutations in diluted real libraries and in public datasets.
- Animal-level testing avoids the gross false-positive inflation that pooling reads produces.
- The Sanger reanalyses show the counting is faithful.

**Uncertain:**
- Transfer to independently prepared runs: every calibration result is conditional on three controls from one run, and the same controls set both the margin and the threshold.
- How much of the gain comes from using controls at all versus from the specific statistical model (§6.2).
- Absolute comparability across platforms, which is not claimed.
- Precision of power estimates from three mice.

**Does anything unresolved materially weaken the paper?** No result I checked is wrong. H1 is a figure regression in Codex's package, not in the manuscript. The prose problems are overstatements and wording errors that a few sentences fix. They do matter, because they are the kind a referee quotes back. The independent-run gap is real and already disclosed. For a methods paper, the ablation/baseline point is the one I would close before submission.
