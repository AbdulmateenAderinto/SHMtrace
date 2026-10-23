# SHMtrace v1.3 review round (21–23 October 2026) — index

This folder holds one review round of the SHMtrace Cell Reports Methods manuscript: a review of v1.3, a
revision by Codex, and two revisions by Claude. **The current version is V2, in `claude_figures_v2/`.**
Everything else is the record of how it was reached, kept so each step can be inspected and rebuilt.

**Start here:** [`claude_figures_v2/READ_ME_FIRST.md`](claude_figures_v2/READ_ME_FIRST.md), then
[`claude_figures_v2/CHANGELOG.md`](claude_figures_v2/CHANGELOG.md) for every change and the reason for it.

## How the versions follow from each other

```
v1.3 (canonical sources, docs/calibration/crm/, commit ef9240d plus uncommitted changes)
 └─ REVIEW.md + PROPOSED_REVISION.md ............ review of v1.3 and an editorial proposal (21 Oct)
     └─ codex_version/ .......................... Codex 1: the revised manuscript (22 Oct)
         ├─ codex_figures_v2/ ................... Codex's figure revision (22 Oct)
         └─ claude_figures_v1/ .................. Claude figure revision 1, built on Codex 1 (22 Oct)
             └─ claude_figures_v2/ .............. V2 – for Alberto review  ← CURRENT (22–23 Oct)
```

Each revision was made as a new folder; no earlier version was edited in place.

## What V2 needs to rebuild

The V2 build (`claude_figures_v2/manuscript/revise_operations.py`, then `author_docx.py`) reads from outside
its own folder, so these must stay where they are:

| Needed | Path | In git? |
|---|---|---|
| Document builder, reference list | `codex_version/sources/docbuild.py`, `codex_version/sources/refs.json` | yes |
| The Codex 1 document that V2's steps are applied to | `codex_version/document_operations.json` | yes |
| Analysis outputs the figures are drawn from | `docs/calibration/crm/out/` (1.9 GB) | **no** — gitignored; regenerate from the analysis scripts or keep a separate backup |
| Per-library calls, trees and Sanger calls | `docs/calibration/figures/`, `docs/calibration/manuscript_corrected/sanger_calls/` | yes |
| SciCanvas builder | `docs/calibration/crm/scicanvas/build.js`, plus the SciCanvas app (`~/SciCanvas-merge`, its own repository) | yes |
| Python with SHMtrace and python-docx | the repository's `.venv` | no (environment) |

Paths are absolute, under `/Users/abdulmateenaderinto/`, in two places: `claude_figures_v2/src/style.py`, and the
image paths stored in the recorded document operations (`codex_version/document_operations.json` and
`claude_figures_v2/manuscript/document_operations.json`). On another machine, edit `style.py` and rewrite those
image paths before rebuilding.

## Contents

### Version folders

| Folder | What it is |
|---|---|
| `claude_figures_v2/` | **Current.** V2 – for Alberto review: manuscript (docx, PDF in `render/`), figures and their scripts (`src/`, `scicanvas/`), the build steps (`manuscript/`), graphical abstract, highlights/eTOC/Motivation, `READ_ME_FIRST.md`, `CHANGELOG.md`. `scratch/` holds intermediate versions referred to in the changelog. |
| `claude_figures_v1/` | Claude's first figure revision, built on Codex 1. Superseded by V2; has its own `CHANGELOG.md`. |
| `codex_version/` | Codex 1: the revised manuscript (`…_v1.3_Codex_1.docx`), the recorded document operations V2 starts from, and `sources/` (document builder and references). Page contact sheets (`contact-*.jpg`) for review. |
| `codex_figures_v2/` | Codex's figure revision: figures, inputs and their hashes, `README.md`, `CHECKS.md`, `FIGURE_INTEGRATION.md`, a review PDF and `index.html`. |
| `render/` | v1.3 as reviewed: the rendered PDF and one image per page. |

### Review notes

| File | What it is |
|---|---|
| `REVIEW.md` | Review of v1.3 against the accepted Word manuscript (21 Oct). |
| `PROPOSED_REVISION.md` | An editorial proposal taking v1.3 from 9,139 to 6,803 words; not applied as written. |
| `DETAIL_RELOCATION_LEDGER.md` | Every numerical placeholder the proposal removed, and where it is still covered. |
| `DISCUSSION_ENDINGS.md` | Proposed final sentences for the Discussion paragraphs. |
| `CLAUDE_FIGURE_REVISION_HANDOFF.md` | The figure critique that Claude figure revisions 1 and 2 respond to. |
| `CLAUDE_INSPECT_CODEX_FIGURE_REVISION.md` | Handoff asking Claude to inspect Codex's figure revision. |
| `CLAUDE_HANDOFF_CODEX_WORK_AND_RIGOUR_REVIEW.md` | Handoff asking Claude to review Codex's manuscript and figure work. |
| `CLAUDE_INDEPENDENT_REVIEW_OF_CODEX_WORK.md` | Claude's read-only review of the Codex work (22 Oct). |

### Scripts behind the review (read-only; none edits the manuscript sources)

| File | What it does |
|---|---|
| `audit_snapshot.py` | Manuscript audit; writes `audit_snapshot.json`. |
| `evaluate_proposal.py`, `proposed_edits.py` | Count and export the editorial proposal against frozen v1.3 sources. |
| `prepare_codex_version.py`, `author_codex_docx.py` | Prepare and author the isolated Codex 1 manuscript. |
| `reference_and_package_checks.py` | Structural checks on references and the package; writes `reference_and_package_checks.json`. |
| `render_codex_svg.cjs` | Renders Codex's SVG figures. It loads `sharp` from a Codex runtime path on this machine, so it is not portable. |

### Data captured during the review

| File | What it is |
|---|---|
| `numbers.json` | The analysis-derived numbers v1.3's placeholders were filled from. |
| `references.json` | The v1.3 reference list, numbered. |
| `source_blocks.json`, `filled_blocks.json` | v1.3 text blocks with placeholders, and filled. |
| `proposal_counts.json`, `proposal_source_blocks.json` | Word counts and source blocks for the editorial proposal. |
| `audit_snapshot.json`, `reference_and_package_checks.json` | Outputs of the two check scripts. |
| `accepted_complete.txt`, `pdf_pages.txt` | Text of the accepted v1.3 manuscript, whole and page by page. |
| `source_vs_article.diff`, `source_vs_complete.diff` | Empty: the comparisons they record found no differences. |
| `supp_contact_1–5.jpg` | Contact sheets of the supplementary pages, used for review. |
| `SHMtrace_Codex_Figures_v2_sources_and_exports.zip` | A packaged copy of `codex_figures_v2/` (330 files). |

## Not in this folder

- **Raw data** (J_H4 Sanger traces, J_H4 and N114 nanopore reads): the data deposit
  `ONT_SHM_DataDeposit_Ver8_2026-07-30` in iCloud Drive (`~/Documents/Projects/2026/…/submission_v4/`), with its
  own README and integrity manifests. `~/Nanoshm/n114_fastq` links into it.
- **Canonical manuscript sources:** `docs/calibration/crm/manuscript/`. V2's edits have not yet been ported
  back there.
- **Demo dataset:** `examples_demo/SHMtrace_demo_JH4/`.

## Where this is saved

Committed on branch `v2-alberto-review-2026-10-23` and pushed to the private GitHub repository
(`AbdulmateenAderinto/SHMtrace`), 23 October 2026. Before that repository is made public, remove this branch
or these files from it: they include the unpublished manuscript, drafts with comments, and unpublished data.
