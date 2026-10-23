# Session status — 2026-10-22

State of the SHMtrace Cell Reports Methods revision at the end of this session.

## Manuscript

`SHMtrace_CellReportsMethods_V2_for_Alberto_review.docx` — the current version, labelled
**V2 – for Alberto review**. 6,998 of the 7,000 allowed words (the cap excludes STAR Methods,
supplemental legends and references). Rendered to `render/SHMtrace_CellReportsMethods_V2_for_Alberto_review.pdf`,
90 pages.

What V2 carries over V1.3:

- the figure revision from `CLAUDE_FIGURE_REVISION_HANDOFF.md`, plus the numerical corrections
  Codex caught that the first Claude pass had missed;
- the author's docx comments: Figure 4 orientation and label clashes, OT1 read as a CH12 AID
  off-target, Figure 5C bars and the "AID-null" wording, the Figure S1 split with S1E kept as
  page 2, S9 legends restored, r² added to S17C, and the text wrapping fixed;
- restored transitions in the main text, with cuts B, C and D1 taken to stay under the cap;
- Philip Barbulescu added as an author for the lab's Sanger sequencing, with his contribution
  noted and a present-address footnote (Structural Genomics Consortium, Toronto, ON M5G 1L7, Canada);
- the two new spike-in baselines in Table S1, the Results and the STAR Methods.

Preserved, not to be edited: the author's own copy
(`SHMtrace_CellReportsMethods_complete_v1.3_Codex1_ClaudeFigures2.docx`), its comment-bearing
twin, and `manuscript_corrected/` in the calibration tree.

`CHANGELOG.md` has the full per-section account; `READ_ME_FIRST.md` is the orientation note.

## How the manuscript is built

`document_operations.json` is replayed by `manuscript/author_docx.py` against
`codex_version/sources/docbuild`. Edits are scripted steps in `manuscript/revise_operations.py`
(steps 7 and 7b–7i) and `manuscript/flow_edits.py`; step 7f checks that every number in a flow
edit already appears in v1.3 or the current text. `manuscript/word_count.py` mirrors the journal's
counter, and the word-limit assertion stays in.

Figures: matplotlib scripts in `src/` write into `figures/`; SciCanvas builders in `scicanvas/`
are built with `docs/calibration/crm/scicanvas/build.js`, after which `_preview.png` is copied to
the name the manuscript expects. `print_check` enforces text ≥ 6 pt and strokes of 0.5–1.5 pt.
Figure S22 was left as it was, deliberately.

## Code

Committed as `7606a03` on branch `audit-fixes-2026-10-19` — five files only, the rest of the
pre-existing uncommitted work left unstaged:

- `src/shmtrace/variants.py` — note that pbinom-with-controls uses the control rate;
- `tests/test_error_models.py` — the test for it;
- `docs/calibration/crm/spikein.py` — the `ctlnull_bh` and `qualonly_pbinom` ablation profiles;
- `make_tables.py` and `figs/style.py` — corrected labels.

384 tests pass.

The v0.1 and pbinom spike-in modes were identical because the per-position control null takes
precedence over pbinom; the real baselines are now quality-only at 5.7–46.6 false calls per
library, `ctlnull_bh` at 1.3–2.8 with 27–64 % recall, against SHMtrace's 63–90 % recall.

## Demo dataset

`nanoshm_v0.2/examples_demo/SHMtrace_demo_JH4/` and its 51.8 MB zip, uncommitted. Six mice
(WT2–WT4, KO1–KO3) read two ways: nanopore J_H4 reads and Sanger clones of the same amplicon,
with the references, the sample sheets and an `expected/` tree.

Verified from a fresh unzip: outputs match `expected/` byte for byte. Nanopore 221/344/73 calls
in the wild-type mice and 0/2/0 in the controls, p = 0.0027. Sanger 134 clones with 33 QC flags
(2 duplicate_clone, 27 label_mismatch, 4 load_outlier) and 2 sheet exclusions.

The bookkeeping flags were left in place on purpose — the README describes them, and showing
rather than silently resolving them is the point.

## Open items

Offered and not yet taken up:

- point `docs/TRIAL.md` and the main README at the demo;
- port the V2 edits back into the canonical sources before submission;
- fill the remaining placeholders — accessions, DOIs, catalogue numbers, RRIDs — none of which
  should be invented;
- decide where to host the demo, given that it is a pre-publication data release. That is the
  author's and PI's call.

## Standing constraints

Never modify `manuscript_corrected/`. Do not redistribute the SHazaM HH_S5F model (AGPL) or the
Karst reads (GPL-3). Do not invent catalogue numbers, RRIDs, accessions or DOIs. Never hard-code
a number in place of a placeholder or change a value. Keep the word-limit assertion. Preserve the
originals and the Codex versions, making new, clearly named versions instead. Commit only when asked.
