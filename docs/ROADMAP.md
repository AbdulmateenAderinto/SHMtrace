# Building an open-source SHM analysis platform — development roadmap

**Goal.** A web platform that performs somatic-hypermutation (SHM) variant
analysis on Oxford Nanopore amplicon reads, using open-source tooling
(minimap2) instead of Geneious, with B-cell SHM as the default workflow (JH4
intron or immunoglobulin V region) and user-adjustable parameters. Users upload
FASTQ reads, a reference FASTA, and optional contaminant FASTA(s); the platform
aligns, filters, calls variants under the paper's rules, and reports mutation
frequency / density and substitution spectra.

**Decisions carried into this plan** (from our discussion): build **local-first
but architected to be hosted**; target the **wider B-cell/SHM community** (so:
strong defaults, guardrails, documentation); deliverable now is a **roadmap plus
a working prototype scaffold** (the prototype lives alongside this file).

---

## 1. Why this is feasible

Everything Geneious does in this workflow has a mature open-source equivalent.
The scientific content is not the alignment — it is the *filtering and counting
rules*, which are already fully specified in the manuscript and reimplemented in
the prototype `src/shmtrace/`.

| Step | Geneious (current) | Open-source replacement |
|---|---|---|
| Map reads to reference | Geneious mapper (high sensitivity) | **minimap2** `map-ont` (via `mappy`) |
| Competitive contaminant removal | manual multi-reference map | index target + contaminant together, keep best hit |
| Chimera / structural filtering | manual | drop supplementary alignments |
| Expected-error / quality filter | Geneious trim | Phred → expected-error rate ≤ 0.01/base |
| Depth normalisation | manual subsample | deterministic random downsample |
| Variant calling | Geneious "Find Variations" | pileup + per-column statistical test |
| Germline / high-freq removal | variant frequency < 50% | same rule, in `variants.py` |
| p-value filter | Geneious approximate p ≤ 0.05 | quality-aware Poisson/binomial test |
| Spectrum / density | Prism / manual | `metrics.py` |

The one genuinely hard, research-grade problem is the **variant-calling error
model** — separating true low-frequency SHM from systematic Nanopore error. The
prototype has a defensible first version; §6 and Phase 1 make calibrating it the
first serious milestone.

---

## 2. Architecture (local-first, hostable)

Three layers, deliberately decoupled so the same core serves a laptop tool and a
public website:

```
        ┌──────────────────────────────────────────────┐
        │  shmtrace/ (pure Python library, no web deps)  │  ← the science
        │  align · filter · call variants · metrics     │
        └──────────────────────────────────────────────┘
              ▲                              ▲
              │ import                       │ import
     ┌────────┴────────┐            ┌────────┴─────────────┐
     │  CLI            │            │  web layer (FastAPI) │
     │  shmtrace.pipeline│            │  upload + params     │
     └─────────────────┘            └──────────────────────┘
                                             ▲
                                    ┌────────┴─────────┐
                                    │ browser SPA      │ ← sliders, file inputs,
                                    │ (static HTML/JS) │   results tables/plots
                                    └──────────────────┘
```

The rule that makes "local now, hosted later" cheap: **the web layer never does
science.** It only marshals uploads + parameters into an `shmtrace.pipeline.run()`
call. Going from a single-user local server to a multi-user hosted service means
inserting a job queue *between* the web layer and the engine — the engine code
does not change.

- **Local mode (Phase 2):** `uvicorn webapp.main:app` on the user's machine or a
  lab server. Files never leave the network. Big FASTQs are fine — it is the
  user's own disk/CPU. This is the recommended default for a lab.
- **Hosted mode (Phase 3):** the same app behind a queue (Redis + a worker), an
  object store for uploads (S3-compatible), and per-job temp isolation. Needed
  for a public community tool where you cannot trust the client's compute.

---

## 3. Recommended tech stack

- **Engine:** Python 3.11+, `mappy` (minimap2), `numpy`; optional `pysam` for
  BAM/CRAM input. No framework lock-in — it is a plain library.
- **Web backend:** **FastAPI** + `uvicorn`. Async, typed, tiny, great file-upload
  and OpenAPI support; trivially containerised.
- **Frontend:** start with a single static HTML page (as in the prototype) —
  zero build step, easy to reason about. Graduate to **React + Vite** only when
  the UI grows (saved runs, multi-sample projects, richer plots).
- **Job queue (Phase 3):** **Redis + RQ** (simplest) or Celery; a Dramatiq
  worker is also fine. One worker process per CPU core.
- **Storage (Phase 3):** S3-compatible object store (AWS S3, or self-hosted
  MinIO) for uploads/results; Postgres for job metadata.
- **Packaging:** Docker image bundling minimap2 + the app, so "local install"
  for a non-bioinformatician is `docker run`. Publish to PyPI as `shm-tools`
  for the CLI crowd.
- **Plots:** compute spectra/density in the engine; render client-side with a
  light charting lib, or reuse the lab's existing matplotlib figure scripts for
  publication-quality exports.

---

## 4. Parameter and preset model

Every UI control maps 1:1 to a field in `shmtrace.config.Params`. A **preset** is
just that dataclass serialised to YAML with assay-appropriate defaults; the UI
loads a preset, then overlays user overrides. This is what lets the tool be
"B-cell SHM by default, but fully tunable."

Shipped presets:
- `jh4_intron` — polyclonal, neutral; germline unknown ⇒ <50% high-frequency
  removal; 529-bp (or full-length) window.
- `n114_vregion` — monoclonal; exact parental reference ⇒ precise germline;
  20-bp end-trim; 615-bp window.

Adding a new assay (e.g. human V(D)J, a different locus) = one new YAML file plus
a reference FASTA. No code change. This is the extensibility story for the
community: the platform is not hard-wired to JH4/N114.

Exposed controls (sliders / text boxes): downsample depth, max expected-error
rate, max variant p-value, germline frequency cutoff, min supporting reads,
end-trim bp, analysis window (start–end), indels on/off, minimap2 preset.

---

## 5. Phased plan with effort estimates

Estimates assume one developer comfortable in Python; halve calendar time with
two. "Done" markers reflect the delivered prototype.

**Phase 0 — Prototype (DONE).**
Working engine + CLI + minimal web UI + two presets + example data.
Demonstrates the full path FASTQ→variants and reproduces WT≫KO. *Effort: done.*

**Phase 1 — Scientific validation & error-model calibration (DONE, July–October 2026).**
Delivered: the quality-calibrated default chosen on J꜀4 and N114 benchmarks and
spike-ins (`docs/calibration/ERROR_MODELS.md`, `GENEIOUS.md`), aggregate
estimation of single-read mutations (`SINGLE_READ.md`), opt-in indel calling
(`INDELS.md`), and validation on data the model had never seen: a public
nanopore deposit (`EXTERNAL_VALIDATION.md`), an Illumina J꜀4 study whose
positive and negative findings both reproduce (`ILLUMINA.md`), an ectopic-AID
versus empty-vector set (`ECTOPIC_AID.md`), and the lab's own matched Sanger
clones (`SANGER.md`). The original Phase 1 description follows.

*(original plan)*
The gating milestone. Run the engine on every sample in the paper and compare
variant counts, densities, and spectra against the Geneious/published values.
Tune the per-column variant test (quality-aware Poisson today; candidates:
empirical per-context error table learned from Aicda-/- controls, strand-bias
filter, min-quality gate) until the Aicda-/- background and WT/KO fold-change
match the paper within tolerance. Deliverable: a validation report + a locked,
documented error model. *This is where the real research effort is.*

**Phase 2 — Local application hardening (mostly done).**
Delivered: the desktop app, multi-sample projects, BAM/Sanger/Illumina/PacBio
input, publication figures and exports, and user docs. Outstanding: a
`docker run` one-liner, and signing the Mac and Windows builds.

*(original plan)*
Robust upload handling (gz, multi-file, BAM input), input validation and helpful
errors, progress reporting, batch mode (many samples at once), export of
publication figures, a `docker run` one-liner, and a proper results view
(per-sample + comparison). Write user docs. *Ship this to your lab.*

**Phase 3 — Hosted multi-user service (3–5 weeks).**
Job queue + workers, object storage, per-job isolation and resource limits,
upload size/rate limits, auth (even lightweight), a job dashboard, and basic
abuse/security hardening (see §6). Deploy behind HTTPS. *This is the public
tool.*

**Phase 4 — Community release (2–3 weeks + ongoing).**
Preset library for common loci, contributor guide for new assays, tutorial +
worked examples, citable release (Zenodo DOI), and a short methods/application
note. Optional: link out to the lab's figure scripts for downstream plotting.

Rough total to a published community tool: **~2.5–4 months** of focused work
(one dev), Phase 1 being the part that must not be rushed.

---

## 6. Risks and open problems (read before committing)

- **Variant-calling specificity (highest risk).** Per-column error models
  over-call on raw Nanopore because some errors are systematic and
  high-quality. Mitigations to test in Phase 1: min-base-quality gate,
  strand-bias filter, an empirical per-context error table trained on negative
  controls, and/or requiring reproducibility across strands. Budget real
  experiment time here; this is the difference between a demo and a tool people
  trust.
- **No UMIs.** Like the paper, per-read frequencies are depth-standardised
  apparent rates, not counts of independent molecules; PCR duplication/jackpot
  cannot be removed without molecular barcodes. State this in the UI.
- **Reference correctness & framing.** Results are only as good as the reference
  the user uploads; window coordinates must match their reference frame. Add
  validation + a coordinate helper.
- **Large uploads (hosted mode).** FASTQs are hundreds of MB–GB. Needs streamed
  uploads, size caps, resumable transfer, and disk quotas — a real infra cost.
- **Security/abuse (hosted, public).** Untrusted file uploads, compute as a
  service, and DoS via huge jobs. Requires sandboxing, per-user limits, and
  scanning. This is the main reason to keep local mode as the default and gate
  the public service.
- **Reproducibility.** Pin minimap2/engine versions and record them in every
  `summary.json` (already emitted) so a result can always be regenerated.

---

## 7. What was needed to proceed to Phase 1 *(all supplied; Phase 1 is done)*

1. The full sample→FASTQ→reference mapping for the paper's samples (so the
   validation harness can run the whole set), plus the Geneious per-sample
   variant counts to calibrate against.
2. Confirmation of the exact reference sequences and window coordinates per assay
   (JH4 1230-bp frame, N114 655-bp frame) — the prototype uses the embedded
   references from `nanopore_pipeline_v12.py`.
3. A decision on the calibration target: match Geneious counts as closely as
   possible, or define an independent, defensible error model and simply report
   concordance with Geneious.

The prototype in this folder is the concrete starting point for all of the above.
