# SHMtrace V2 — session log and handover (21–25 October 2026)

A record of everything done in the Claude session that produced V2, what was decided and why, the state
the project is in now, and what is still open. It supersedes `SESSION_STATUS_2026-10-22.md`. The detail of
every manuscript and figure change, with reasons, is in [`CHANGELOG.md`](CHANGELOG.md); the map of the whole
review folder is in [`../README.md`](../README.md).

---

## 1. State on 25 October (checked, not remembered)

> **Resolved later on 25 October.** The 25 October corrections (section 4) had left the text 23 words over
> (7,022) and uncommitted, with a stale PDF. Three cuts brought it back to 6,999 (changelog section 12); the PDF
> was re-rendered; and the corrections, the cuts and this log were committed together (not yet pushed).

| | Now |
|---|---|
| Current manuscript | [`SHMtrace_CellReportsMethods_V2_for_Alberto_review.docx`](SHMtrace_CellReportsMethods_V2_for_Alberto_review.docx) (25 Oct) |
| As sent to Alberto | `SHMtrace_CellReportsMethods_V2_as-sent_2026-10-23.docx` (byte-for-byte) |
| PDF | `render/SHMtrace_CellReportsMethods_V2_for_Alberto_review.pdf` — 87 pp., 25 October |
| Word count | 6,999 of 7,000 (limit counts everything but references and STAR Methods); Summary 150 of 150 |
| Display items | 7 main figures; 25 supplementary figures, numbered in order of first citation; 51 references |
| Graphical abstract | [`SHMtrace_V2_GraphicalAbstract.pdf`](SHMtrace_V2_GraphicalAbstract.pdf) / `.png` |
| Highlights, eTOC, Motivation | [`SHMtrace_V2_Highlights_eTOC_Motivation.docx`](SHMtrace_V2_Highlights_eTOC_Motivation.docx) |
| Git branch | `v2-alberto-review-2026-10-23`; the latest commits are local until pushed |
| Changelog | sections 1–12 |

---

## 2. Where everything is

| What | Where | Safe off this Mac? |
|---|---|---|
| V2 manuscript, figures, scripts, build steps | `output/reviews/shmtrace_v1_3_2026-10-21/claude_figures_v2/` | Yes, including the 25 Oct edits |
| Earlier versions and review notes | `output/reviews/shmtrace_v1_3_2026-10-21/` (index: `README.md`) | Yes (`3021c51`, pushed) |
| Demo dataset | `examples_demo/SHMtrace_demo_JH4/` (+ zip) | Yes (`d65e6a3`, pushed) |
| SHMtrace code | `src/shmtrace/` | Yes, through `9b4c994` |
| Canonical manuscript and figure sources | `docs/calibration/crm/` | **No** — 91 staged and 88 modified files uncommitted |
| Analysis outputs every figure reads | `docs/calibration/crm/out/` (1.9 GB) | **No** — gitignored |
| Canonical build outputs | `docs/calibration/crm/build/` (743 MB) | No — gitignored, quick to regenerate |
| Raw data (J_H4 Sanger, J_H4 and N114 nanopore) | iCloud Drive: `~/Documents/Projects/2026/…/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/` (314 MB, 591 files, with README and integrity manifests) | Yes (iCloud only; no second copy) |
| SciCanvas | `~/SciCanvas-merge` (own GitHub repository) | Yes |

This Mac has **no Time Machine backup** and its disk is **99% full** (17 GB free of 926 GB).

---

## 3. What this session did (21–23 October)

### 3.1 Before V2

- **Figure critique** (`../CLAUDE_FIGURE_REVISION_HANDOFF.md`): implemented as separate Claude versions
  (`claude_figures_v1/`, then `claude_figures_v2/`), never editing Codex's work in place.
- **Independent review of Codex's work** (`../CLAUDE_INDEPENDENT_REVIEW_OF_CODEX_WORK.md`). The author
  preferred the Claude figures and asked to bring in the numerical corrections Codex caught.
- **Identical v0.1 and pbinom spike-in results** explained: the per-position control null takes precedence
  over pbinom. A pbinom note and a real quality-only baseline were added; the two baselines went into
  Table S1, the Results and the STAR Methods (quality-only 5.7–46.6 false calls per library; control null
  with BH 1.3–2.8, recall 27–64%; SHMtrace recall 63–90%). Five canonical files committed as **`7606a03`**
  on `audit-fixes-2026-10-19`; 384 tests pass.
- **Author's docx comments** addressed (Figure 4 orientation and clashes, OT1 as a CH12 AID off-target,
  Figure 5C bars, "AID-null" wording, S1 split, S9 legends, r² on S17C, wrapping).
- **SciCanvas figures** revised; Figure S22 (now S24) kept as it was, at the author's request.
- **Flow**: transitions restored in the main text; the author chose cuts B, C and D1 to stay under 7,000.
- **V2 – for Alberto review**: the author's own edits ported; **Philip Barbulescu added as an author** for
  the laboratory's Sanger sequencing, with a present-address note (Structural Genomics Consortium, Toronto).
- **Demo dataset** built and verified from a fresh unzip: nanopore 221/344/73 calls in Wt mice against
  0/2/0 in controls (p = 0.0027); Sanger 134 clones with 33 bookkeeping flags, left in on purpose.

### 3.2 The SSRN preprint (ssrn-7231542)

*High-Confidence Profiling of SHM … Using Oxford Nanopore R10.4 Sequencing* (Aderinto & Martin) is an
earlier version of this study, built on Geneious calling. Overlap with V2:

- **Text:** about 3% of V2's words are in shared runs, mostly references and PCR methods.
- **Data:** complete — all twelve libraries and the Sanger clones of the same six mice.
- **Results:** every preprint result has a counterpart in V2, recomputed with SHMtrace (for example N114
  ~18-fold → 85-fold). The r = 0.989 in both was recomputed from SHMtrace's own calls and matches by
  coincidence (638/817 calls against the preprint's 571/739).

V2 does not cite it. **Recommendation given:** do not have both manuscripts in review at once; withdraw the
preprint at a major-revision decision (or sooner), disclose it as an earlier version, and consider a STAR
Protocols paper for the bench workflow. **Alberto's decision.**

### 3.3 Figures

All matplotlib figures pass the Cell Press print check (text ≥ 6 pt, lines 0.5–1.5 pt). Numbers below are
the **final** numbering; old numbers in brackets.

- **S18 page 2, mutation trees** [S17]: redrawn at page width (`src/figS17_trees.py`); the old page printed at
  about 3.7 pt. The same trees, rebuilt from the delivered reads and checked node for node.
- **Naming**: J_H4 mice are Wt2–Wt4 and KO1–KO3 in every figure; *Aicda*⁻/⁻ italic with tight superscript.
- **Label collisions** fixed in Figure 6 (A, C) and S5A [S4A].
- **S10, N114** [S9]: linear axis kept, "85×, p = 0.0071" restored (computed and checked against Figure 4);
  restyled to match Figures 5 and 7; the From → To matrices replaced by a spectrum panel.
- **Matrices removed everywhere**: per-library and Sanger spectra restyled as bars; then **S10 page 2 and
  the S14 [S13] top panel removed** as exact repeats of Figure 5B and S18A. S14 kept, with its trace accounting.
- **S3, base quality** [S2]: redrawn full width; it was drawn half-width and printed at nearly twice size.
- **S13 legend** [S12]: 188-fold there is a ratio of mutation frequencies; the 85-fold of Figures 4F and S10A
  is a ratio of unique-variant densities.
- **V(D)J** no longer has a bold D in legends (the builder took "(D)" for a panel label).
- **S2, read accounting**: the former S1E–F became its own figure; S2–S24 renumbered S3–S25 everywhere, the
  Results now cite S2, and S23–S25 were reordered so first citations run S1–S25 in order.
- **S25, calling rules** [S23]: all six stages restored at the author's request (read preparation and
  competitive alignment as A–B); rows tightened so figure and legend share a page; legend declares that
  A–B show the same data as Figure S2A–B.
- **Figures are deterministic**: rebuilding Figure 4 and the outputs figure gave byte-identical images.

### 3.4 Checks

- **References**: all 51 exist; every DOI resolves on Crossref and matches author, year, title, volume and
  pages; no retractions; author corrections exist for Seija et al. 2026 (8 Oct 2026) and Virtanen et al. 2020.
- **Cell Press standards**: word limit, seven display items, figure size, font and line specifications, and
  supplementary titles and order all met as checked. The official author pages could not be fetched.
- **Ethics**: data handling sound. The AI-use declaration understates the AI's role; wording was drafted, and
  **the author chose to leave it for now**.

### 3.5 Text for submission

- **Highlights** (84, 82, 82, 83 characters):
  1. Candidate mutations in reads are tested against errors of matched AID-null libraries
  2. AID-null control calls fall 64–90% vs an uncorrected rule; true mutations are kept
  3. Calls feed SHM readouts: AID hotspot targeting, spectra, clones and mutation trees
  4. AID signatures are recovered from Sanger, nanopore and Illumina with low background
- **eTOC**: "Aderinto et al. present…" (43 words).
- **Motivation** (170 words): adds that SHM software has followed B-cell receptor repertoire sequencing, so a
  reporter locus read in bulk has had no equivalent route. Fact-checked: SHMTool's original URL no longer
  resolves (it moved to `shmtool.einsteinmed.edu`; the outage could not be dated); MiXCR documents a nanopore
  preset, hence "has followed", not "is only for". **Nothing about any tool's availability is in the paper.**
- **Graphical abstract**: PacBio added as a fourth platform; "Low background" now reads "64–90% fewer control
  calls than an uncorrected rule". Its numbers are checked against the analysis by `manuscript/highlights_v2.py`.
- **Main text** (step 7n): the Summary opens with the repertoire-software gap, and the Introduction's software
  sentence says the same; the Summary describes the comparator as "a general-purpose sequence-analysis
  platform not tuned for SHM" ("GUI" dropped: SHMtrace ships its own desktop app). Geneious is still named
  18 times in the body.

### 3.6 Saving and storage

- Commits on branch `v2-alberto-review-2026-10-23`, both pushed to the private GitHub repository:
  **`d65e6a3`** (V2 folder and demo, 494 files) and **`3021c51`** (the rest of the review folder, 1,078 files,
  with a `README.md` index). The 91 staged canonical changes were deliberately left out of both.
- Storage audit (section 2); advice given on Time Machine to a NAS, and on keeping a second copy of the
  raw data there (it is only 314 MB).

---

## 4. Work after this session (24–25 October)

Made outside this conversation; recorded here so the log is complete.

| Commit | Date | What |
|---|---|---|
| `a6e322c`, `c0e2b28`, `eb91022` | 24–25 Oct | UMI grouping fix and sweep; manuscript Karst UMI numbers and Methods |
| `f312ce0`, `0561037` | 24–25 Oct | `umi_sim` fix; merge of the UMI grouping branch |
| `d5e1bd1` | 25 Oct | Accuracy measured on every public set with three or more controls |
| `9b4c994` | 25 Oct | Four opt-in caller options (`docs/calibration/CALLER_FIXES.md`, `ACCURACY_BREADTH.md`) |
| `3c4e451` | 25 Oct | Manuscript: the Ramos result is allelic, and what merged qualities are — **not pushed** |

**Uncommitted, in this folder (25 Oct, 15:52):** "Corrections after the Alberto review" — CHANGELOG's second
section 11. In short:

- the Ramos excess is a second IGH allele, not a noisy library;
- merged Illumina reads carry two quality levels, not one;
- the UMI consensus is now 46-fold to 0.21% (105 molecules), not 57-fold to 0.17%;
- `mask_control_het` is documented;
- Figure S9H is redrawn;
- the 277–413 sensitivity number is removed.

These corrections are what put the text **23 words over**. The document as sent to Alberto on 23 October is
preserved as `…_V2_as-sent_2026-10-23.docx`.

---

## 5. Open items

### Now
1. **Close Word without saving.** Its lock file is still present.
2. The text is at 6,999 of 7,000 and the Summary at 150 of 150: any addition needs a matching cut.

### For Alberto to decide
- How to handle the **SSRN preprint** before SHMtrace is submitted (section 3.2).
- **Philip's authorship**: he agrees to it and approves the final version.
- The **scope of the AI-use declaration**.

### Before journal submission
- Fill the 16 placeholders: accessions, DOIs, catalogue numbers, RRIDs, sex of each J_H4 mouse, N114
  authentication and mycoplasma status, buffer and volume details.
- Port V2's edits into the canonical sources (`docs/calibration/crm/manuscript/`): text, figure numbering,
  highlights, eTOC and Motivation.
- Rewrite the AI-use declaration to cover the full scope, as drafted.
- Build some word margin, and check whether Cell Reports Methods caps the number of supplementary figures (25 here).
- Confirm SHMTool's current address in a browser if any URL is printed.
- Check whether the Seija et al. correction touches the claim cited.
- Decide where to host the demo dataset (a pre-publication data release).

### Data safety
- Commit the canonical changes in `docs/calibration/crm/` (91 staged, 88 modified) on their own branch.
- Copy `docs/calibration/crm/out/` (1.9 GB) somewhere off this disk.
- Set up Time Machine (to the NAS, with an SMB Time Machine share, about a 2 TB quota and encryption).
- Copy the raw data deposit to the NAS, and verify it against `04_Manifests_and_Integrity`.
- Free 50–100 GB of disk space.
- Before the GitHub repository is made public, remove this branch: it holds the unpublished manuscript,
  commented drafts and unpublished data.

---

## 6. Decisions the author made

| Decision | Date |
|---|---|
| Separate Claude version, not an edit of Codex's | 21 Oct |
| Keep Figure S22 (now S24) as it was | 22 Oct |
| Cuts B, C and D1 for the word limit | 22 Oct |
| Philip Barbulescu as author, with a present-address note | 22 Oct |
| Leave the demo's bookkeeping flags in | 22 Oct |
| Leave the AI declaration for now | 23 Oct |
| Keep Figure S14 [S13] and its readable-window panel as it is | 23 Oct |
| Remove S10 page 2 and the S14 top panel as redundant | 23 Oct |
| Read accounting as its own Figure S2; renumber; first citations in order | 23 Oct |
| Restore all six stages to the calling-rules figure | 23 Oct |
| Describe the comparator in the Summary instead of naming it; drop "GUI" | 23 Oct |
| Frame the Motivation and main text on the repertoire-software gap; leave SHMTool's outage out | 23 Oct |

---

## 7. Rules followed throughout

- Never modify `docs/calibration/manuscript_corrected/`.
- Do not redistribute the SHazaM HH_S5F model (AGPL) or the Karst reads (GPL-3).
- Never invent catalogue numbers, RRIDs, accessions or DOIs; never type a number in place of a placeholder or
  change a value. Numbers come from the analysis outputs.
- Keep the word-limit check; every edit priced with `manuscript/word_count.py`.
- Preserve originals and earlier versions; make new, clearly named ones.
- Commit only when asked.

---

## 8. How to rebuild V2

From `claude_figures_v2/`, with `V=/Users/abdulmateenaderinto/Nanoshm/nanoshm_v0.2/.venv/bin/python`:

```bash
python3 manuscript/revise_operations.py
$V manuscript/author_docx.py
$V manuscript/word_count.py manuscript/document_operations.json
/Applications/LibreOffice.app/Contents/MacOS/soffice --headless --convert-to pdf --outdir render SHMtrace_CellReportsMethods_V2_for_Alberto_review.docx
$V manuscript/highlights_v2.py
```

Figures: `cd src && python3 <script>.py`, then `python3 print_check.py <FigureName>`. SciCanvas figures and the
graphical abstract build from `docs/calibration/crm/scicanvas/` with `node build.js <builder.js> <out-stem>`
(`GA=1` for the graphical abstract); copy each `_preview.png` to the name the manuscript expects. The build
reads `../codex_version/` and `docs/calibration/crm/out/`, so both must be present.
