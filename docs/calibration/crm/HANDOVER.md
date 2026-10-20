# Handover: the SHMtrace Cell Reports Methods submission

Written for whoever picks this up next. It says where things are, how to rebuild them, what must not be
touched, and what is still open. `README.md` in this folder lists the analyses; this file covers the parts
that are easy to get wrong.

Repository root is `~/Nanoshm/nanoshm_v0.2`. Everything below is relative to `docs/calibration/` unless
said otherwise. Updated after the 17 Oct co-author review was adopted (patch `3317fd9` and follow-up commits).

---

## 1. Use the project virtualenv, not the system Python

```bash
~/Nanoshm/nanoshm_v0.2/.venv/bin/python
```

This matters. The venv has **matplotlib 3.11**; the system interpreter (miniforge, 3.13) has **3.10**,
which cannot render `\underline` and so **fails to build Figure 2**. A failed `savefig` can also leave a
truncated PDF in `crm/figures/`, which looks like a corrupt committed file. If that happens:
`git checkout -- docs/calibration/crm/figures/<name>.pdf`.

Most scripts also want `-W ignore`; scipy emits a harmless precision warning on near-identical data.

```bash
PY=~/Nanoshm/nanoshm_v0.2/.venv/bin/python
$PY -W ignore crm/manuscript/build_crm.py crm/build
```

---

## 2. Layout

| Path | What is there |
|---|---|
| `crm/manuscript/` | the article. `crm_text.py` (Results/Discussion), `crm_methods.py` (STAR Methods **and all figure legends**), `crm_numbers.py` (every `{placeholder}`), `build_crm.py` (assembles, holds the supplementary figure list and Tables S3–S5), `lab_details.py` (author-only details), `docbuild.py` (docx primitives) |
| `crm/figs/` | the matplotlib figures: `fig2_call.py` … `fig7_platforms.py`, `figS2`, `figS3`, `figS5`, `figS6`, `figS8`, `figS9`, `figS15`. `style.py` holds the shared palette, `FULL`/`HALF` widths and `save()` |
| `crm/scicanvas/` | the SciCanvas figures (section 4) |
| `crm/out/` | every analysis result as JSON. **Gitignored** — not in the repo, lives only on this machine |
| `crm/figures/` | built figures (PNG/PDF/SVG/TIF), committed |
| `crm/build/` | the built `.docx` files and the reviewer zip |
| `crm/REVIEW_RESPONSE.md` | the issue-by-issue record of the external review. Read this first |
| `crm/ARCHIVE.md`, `ARCHIVE_ENVIRONMENT.txt` | which command produced each result, and the machine/versions |
| `manuscript_corrected/` | **the earlier JIM submission. Do not modify** (section 3) |

The text and the numbers cannot drift: every result in the prose is a `{placeholder}` filled by
`crm_numbers.py`, and the build aborts on any unfilled one.

---

## 3. Hard constraints

1. **Never modify `manuscript_corrected/`.** That is a separate, earlier submission (the J Immunol paper)
   whose figures must stay as submitted. Supplementary Figures S4, S5 (components), S7 and S11 in *this*
   paper reuse images from it. When one of those needed changing, the fix was to **copy the builder into
   `crm/` and edit the copy** — see `crm/scicanvas/figS1_workflow.js`, which is this article's own copy of
   the detailed workflow figure. Do the same rather than editing the original.
2. **The main text is capped at 7,000 words.** Since the 18 October change to `journal_word_count` (legends and front matter now count) the build reports **7,736**, about 740 over (20 October: 7,538 after the TAM and TCR studies were dropped, +90 for the Figure S19 paragraph, +108 for Figure S20); it writes both documents and then fails the assertion. `build_crm.py` asserts this
   and fails the build if exceeded. The count covers Introduction → Limitations only. **Figure legends,
   STAR Methods, tables and the Motivation are not counted**, so detail belongs there.
   Summary ≤ 150 (now 139), eTOC ≤ 50 (now 43), four highlights ≤ 85 characters each. The documents are
   **US Letter**; `docbuild.figure_block` reads the page height from the document, so do not hard-code one.
3. **Do not invent catalogue numbers, RRIDs, accessions or DOIs.** Fourteen such details are unknown and
   are tracked in `lab_details.py`; the build prints them as "submission blockers" on every run. Leave them
   as placeholders.
4. **Keep statements about the original studies neutral.** Discrepancies found in Phil's / the published
   data go to `~/Documents/Genesifter JH4 data/CATALOGUE.md` and to the authors, not into the paper as
   accusations. The paper states what does and does not match, and why, without blame.
5. **Do not redistribute** the SHazaM HH_S5F model (AGPL) or the Karst et al. reads (GPL-3). They are
   fetched from their own sources; `ARCHIVE.md` says so.
6. Commit trailer in use: `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.

---

## 4. SciCanvas

SciCanvas 2.1 is a desktop figure tool at `~/SciCanvas-merge` (Electron; `npm start` opens the app).
Figures here are built **outside** the app by a Node script, so they are reproducible from source.

### Building

```bash
cd ~/Nanoshm/nanoshm_v0.2/docs/calibration
node crm/scicanvas/build.js crm/scicanvas/<builder>.js crm/figures/<out-stem>
```

Environment variables `build.js` reads:

- `SCICANVAS` — the app checkout. Defaults to `~/SciCanvas-merge`.
- `PLOT_PYTHON` — the interpreter for code plots. Defaults to the project `.venv`.
- `GA=1` — relaxes the Cell Press limits for a graphical abstract (text up to 18 pt, strokes up to 2.5 pt).
  **The graphical abstract will not build without it.**

Each build writes four files from one stem: `<stem>.scifig` (reopen in the app to edit), `<stem>.svg`,
`<stem>.pdf` (vector, fonts embedded, via headless Chrome) and `<stem>_1000dpi.tif`. It also writes
`<stem>_preview.png` — see the gotcha below.

### Builders and their outputs

| Builder | Figure | Data injected by |
|---|---|---|
| `crm_figure1.js` | `Figure1` | — |
| `crm_ga.js` | `GraphicalAbstract` (needs `GA=1`) | — |
| `figS1_workflow.js` | `FigureS1_workflow` | read accounting, by `build.js` |
| `figS8_cohort2.js` | `FigureS8_sanger_cohort2` | `cohort2_data.py` |
| `figS9_mkln.js` | `FigureS9_sanger_mkln` | `mkln_data.py` |
| `figS10_swhel.js` | `FigureS10_swhel` | `swhel_data.py` |
| `figS16_switch_pacbio.py`, `figS17_public_smu.py` | Figures S16 and S17 | matplotlib, not SciCanvas; see `PUBLIC_VALIDATION.md` |
| `figS12_platforms.js` | `FigureS12_platforms` | `s12_reads.py` (prints a READS array to paste); the Illumina band's numbers are from `s12_illumina.py` (→ `out/s12_illumina.json`), typed into the builder |
| `figS13_overview.js` | `FigureS13_overview` | — |
| `figS14_engine.js` | `FigureS14_engine` | `engine_data.py` |
| `star_figure1.js` | the STAR protocol's figure | — |


**Print minimums for the matplotlib figures.** Every figure drawn by `figs/*.py` meets Cell Press's production
rules: text at least 6 pt (super- and subscripts included: `style.setup()` keeps them at their text's size),
lines 0.5-1.5 pt. `python figs/print_check.py` reads each figure's SVG at print size and lists any breach; run
it after changing a figure (the STAR protocol's Figures 3 and 4 are checked too). Those sizes are set in the
scripts themselves, so the figures stay generated from the analysis outputs.

**A final hand-finishing pass in SciCanvas.** Once the numbers are frozen, `scicanvas/mpl_to_scicanvas.py
<script>` turns a `figs/*.py` figure into a SciCanvas builder: one code plot per panel (the numbers written into
its code), the panel letters and titles as SciCanvas text, fixed margins so the axes sit where matplotlib put
them. Build it with `node scicanvas/build.js scicanvas/<script>.js figures/<Figure>_scicanvas` and edit the
`.scifig`. Do this last: the converted copy no longer follows the data, and re-running the converter after new
analyses replaces anything edited in SciCanvas.

### How a builder works

A builder is a single `() => page` arrow function written against SciCanvas's own `Make` API — the same
form the app runs from its console, so it can be pasted straight into SciCanvas. The API used here is
`Make.text`, `Make.rect`, `Make.ellipse`, `Make.icon`, `Make.image`, `Make.connector`, `Make.badge`.
Builders are **sandboxed and read no files**, which is why numbers are injected into the source.

**Injecting numbers.** `inject.py` writes a JSON blob into the builder between `/*DATA*/ … /*END*/`
markers. The data scripts (`cohort2_data.py` etc.) compute the numbers from `crm/out/` and call it, so the
figure and the manuscript quote the same values. To refresh a figure's numbers, run its data script, then
`build.js`.

**Code plots.** `PLOT({ code, x, y, w, h, name })` embeds a matplotlib plot. The code runs through
SciCanvas's own plot worker at printed size, after the shared prelude `plots_common.py` (palette, rcParams),
which `build.js` saves as the document's plot prelude. The result is placed as an SVG asset *plus its code*,
so reopening the `.scifig` in SciCanvas keeps every plot re-runnable and every label editable. Only
`figS9_mkln.js` and `figS10_swhel.js` use this.

**Cell Press checks run on every build and a breach stops it**: 174 mm wide, ≤ 200 mm tall, Arial, text
6–8 pt, strokes 0.5–1.5 pt, grey fills 10–80%. This is a feature — if the build fails on "text 5.40 pt",
the figure would have been rejected in production. Fix the figure, do not raise the limit.

### Two gotchas

1. **`build.js` writes `<stem>_preview.png`, but the manuscript references different names for the older
   figures** — e.g. `build_crm.py` points at `FigureS10_swhel.png`, not `FigureS10_swhel_preview.png`.
   After rebuilding, make sure the file the manuscript references is the one you just made. For
   `FigureS12_platforms` and `FigureS13_overview` that means copying `<stem>_preview.png` to `<stem>.png` and
   saving `<stem>.tif` from it (1,000 dpi label, LZW); on 17 Oct both had been rebuilt but the documents still
   embedded the old images until this was done. The existence check below does not catch a stale file:

   ```bash
   $PY -W ignore -c "
   import sys; sys.path.insert(0,'crm/manuscript')
   from pathlib import Path
   import build_crm as b
   for t, img, *_ in b.SUPP_FIGS:
       for p in (img if isinstance(img, list) else [img]):
           if not Path(p).exists(): print('MISSING', t.split('.')[0], p)
   print('checked', len(b.SUPP_FIGS))"
   ```

2. `crm/README.md` still says the SciCanvas figures are built with
   `manuscript_corrected/figure1_scicanvas/build.js`. There is now a local `crm/scicanvas/build.js` (the
   same script plus code plots). **Use the local one**; `ARCHIVE.md` has the correct command.

---

## 5. Rebuilding

```bash
PY=~/Nanoshm/nanoshm_v0.2/.venv/bin/python
cd ~/Nanoshm/nanoshm_v0.2/docs/calibration

# a matplotlib figure
cd crm/figs && $PY -W ignore fig3_calibration.py       # writes crm/figures/Figure3.{png,pdf,svg,tif}
# Figure 5E reads out/fig5_per_mouse.json: $PY crm/fig5_per_mouse.py first (it refuses to write unless the
# per-mouse calls pool back to Figure 5B exactly). Figure S7 page 2 is figs/figS7_bookkeeping.py.

# the documents (always rerun after touching text, numbers or legends)
$PY -W ignore crm/manuscript/build_crm.py crm/build

# the archive record and the reviewer package
$PY -W ignore crm/make_archive_manifest.py             # ARCHIVE.md + ARCHIVE_ENVIRONMENT.txt, 21 analyses
$PY -W ignore crm/make_review_package.py               # crm/build/SHMtrace_review_package.zip

# checks
cd ~/Nanoshm/nanoshm_v0.2 && $PY -m pytest -q          # 292 tests
$PY -m ruff check src                                   # clean
```

`ruff check src` is clean and must stay that way. `ruff check .` reports ~466 errors: 320 are in
`manuscript_corrected/scripts/` (the JIM submission, which must not be modified) and the rest are
pre-existing in `docs/calibration/` — mostly unsorted imports and unused loop variables in the older figure
scripts (`fig4_controls.py`, `fig5_sanger.py`, `fig7_platforms.py`, `star_figs.py`, `resolve_refs.py`).
None are new. The working rule: lint `src` plus the files you touch, and leave the rest alone rather than
reformatting files you have no reason to change.

```bash
$PY -m ruff check src <the files you edited>
```

A figure is regenerated with its own script; the long analyses in `crm/out/` do **not** need rerunning
unless the method changed. `ARCHIVE.md` lists the command behind every one.

---

## 6. Where the review stands

`REVIEW_RESPONSE.md` is the authority. Summary after the 17 Oct co-author review:

- **All 12 major issues**: addressed. Issues 3 and 4 ended with `crm/fdr_validation.py` + `fdr_estimator.py`
  + `fdr_null_reps.py`, 220 runs of `run_project` on resampled AID-null reads, reported as **Figure S15**. The
  co-author review found four defects in that harness; they are fixed and every run was regenerated. The
  numbers in the text come from `validation_summary.py`, not from prose.
- **All 31 numbered inconsistencies (C01–C31)**: addressed, including C29 (Table S2's About sheet).
- **All 15 minors**: closed, except what only the authors can supply.
- **Figure-by-figure plan**: complete, Figures 1–7 and S1–S15.

### Genuinely still open

1. **One independently prepared validation run.** The co-author review's main remaining suggestion: the
   stress tests resample three existing AID-null libraries from one run, so they test the procedure
   conditionally, not its transfer to a new preparation.
2. **Issue 11** — biological QC flags suspicious calls but does not prove origin; a clone-aware recovery
   step was discussed but not built. (line ~153.)
3. **The 14 submission blockers** in `lab_details.py`. Printed on every build. Four are deposit
   identifiers (SRA BioProject, data DOI, Zenodo DOI, companion citation); the rest are catalogue numbers,
   RRIDs, bench details, per-mouse sex and the N114 authentication/mycoplasma status.
4. **The AI-use declaration** now names Anthropic Claude and OpenAI Codex and says the assistance covered
   "drafting the analysis and figure code, language editing and figure layout". That is the authors' statement to make — it should be confirmed, not
   inherited.
5. **`switch.py`'s ambiguity flag is inert on real data** (17 Oct, confirmed defect, deliberately not
   patched). `summarise()` reports `ambiguous_pct`, and it reads 0.0% on all twelve PacBio libraries of
   PRJNA831666 — not because breakpoints are unambiguous, but because the check needs a secondary alignment
   that minimap2 does not return for these reads. Measured directly, 19–47% of unique junctions place equally
   well at more than one acceptor position. `tests/test_switch.py` does not cover the flag. The manuscript
   never quotes the built-in figure; Figure S16D reports the direct measure instead. A fix is proposed in
   `PUBLIC_VALIDATION.md`, not applied.
6. **Why Figure 5E is not per-molecule mutation load** (17 Oct, a note rather than a problem). Both platforms
   rank the wild-type mice the same way by mean mutations per molecule (Wt3 > Wt4 > Wt2: Sanger 9.7, 5.8, 4.5
   per clone; nanopore 11.9, 2.6, 1.5 per read), but the distributions differ: the unmutated share is 51%,
   37% and 42% of Sanger clones against 75%, 8% and 63% of nanopore reads. Reads sample clones in proportion
   to their size (Wt2's largest clone is largely unmutated; Figure 6C), and 26-40 clones per mouse sample
   each mouse thinly, so a per-molecule load panel would mostly show sampling. Substitution classes are a
   share of each mouse's own calls and do compare, which is why 5E shows those.
7. **Two public Illumina datasets carry a flat Q30.** PRJNA973112 and PRJNA918596 as distributed by SRA
   have Q30 on every base, so their error model has one quality level. STAR Methods and the Figure S12 legend
   now say so; no result changed. If the original qualities can be obtained (the submitters' files rather
   than SRA's simplified ones), rerunning those two would test the quality strata on Illumina.
8. **The reviewer's proposed replacement Summary** (in the review, front-matter section) was deliberately
   not adopted; the current Summary is 149 words and already carries the two-route fix. The authors' call.

### Two things to know about the method, because the text depends on them

- **There are two calling routes.** Reads go through the control-calibrated test with a target–decoy FDR
  estimate; **Sanger clones are counted at Q ≥ 20 with the background reported beside them, never
  subtracted**. Figure 1C is a fork for this reason. Any sentence implying one procedure for all platforms
  is wrong.
- **`run_project` and `benchmark_real.py` are different procedures, and both select a control-based
  threshold.** `run_project` is the method and is authoritative for biology: one project-level margin and
  decoy collection, the tested control's own decoy removed, and the noisy-decoy exclusion rule.
  `benchmark_real.py` is a caller-comparison harness that calibrates the margin and threshold within each
  library's peer-control set and applies no noisy-decoy exclusion. Figure 4's panels come from the latter and its
  legend says so. Do not mix their numbers.

---

## 7. Conventions worth keeping

- Figure legends live in `crm_methods.py` (main figures) and `build_crm.py` (supplementary). They are not
  word-capped, so when the main text is near its 7,000-word cap, detail goes to the legend.
- Every figure change is verified by *looking at the rendered PNG*, not by assuming the code is right.
  Several real problems here were only visible that way: labels overlapping data, a clipped annotation, a
  node count drawn inside its own circle, a legend running off the page.
- When a claim is checked against the code or the data and turns out to be wrong, the fix is to correct the
  claim, not to soften it. Example: "did not reach significance" (the result was significant, p = 0.002).
- Check a negative claim with a search broad enough to find the thing under any name. "Repair-pathway
  indices exist nowhere in the code" was asserted from a search for the word "pathway" in two files;
  `src/shmtrace/repair.py` exists and calls them repair-route indices. It is not used by the project
  workflow, which is the accurate reason the paper does not report it.
- Check assumptions about the environment before building on them: pagination was sized for A4 when the
  documents are US Letter. Read such values from the source rather than typing them.
- A test harness gets the same scrutiny as the method. The FDR stress test shipped with four defects that a
  regression test per behaviour would have caught; `tests/test_crm_review.py` now holds those tests.
