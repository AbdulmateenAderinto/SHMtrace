# SHMtrace audit corrections — October 19, 2026

Baseline: `c3a633ee78a5bb4b024a2ffc8752162749ecfbf7`.

This is a bounded corrective pass following the computational methods review.
It changes implementation consistency and guards unsupported analyses. It does
not introduce a new caller, change presets, enable the singleton quality floor,
change nominal alpha, or regenerate any manuscript figure or result table.

## Changes and scientific rationale

### 1. Sample/decoy candidate parity — audit F01

`variants.py` now uses one candidate-exclusion helper for sample calls and the
quality-calibrated held-out-control calculation. Both respect minimum support,
the optional singleton quality floor, per-allele frequency exclusions and the
opt-in ancestral-reference policy. The Sanger-specific clone-depth safeguards
remain in the sample route; this helper's calibration route is for qcal only.
Recalibrated qualities are used consistently when recalibration is enabled.

Previously a Q10 singleton excluded by a Q30 sample floor remained a decoy
candidate, while a fixed alternate allele allowed by ancestral mode was omitted
from decoys. Both discrepancies are corrected. Margin selection now uses that
same candidate family. With the quality floor off and ancestral mode off, this
refactoring is intended to preserve candidate selection.

`crm/spikein.py` keys its in-memory null cache by all Params fields. Two profiles
with different candidate policies therefore receive separate fits, rather than
silently sharing the first profile's calibration. A toy harness integration test
verifies that the Q30 profile receives its own fitted null for each held-out
control. Old quality-floor sweep outputs must be rerun before interpretation.

### 2. Saved-panel fidelity and visible fallback — audit F02

Null-panel format version 3 preserves named decoy p values and their exposure,
as well as the existing rates, multiplier, recalibration and metadata. A real
JSON-file round trip is tested for identical call sets, calling mode, threshold
and calibration state on identical inputs.

Version 1/2 rate-only panels remain readable. Requested control thresholding
without decoys retains the existing BH fallback but now names that change in
calling notes and summary warnings. Project summaries deduplicate this warning.
No claim of BH validity for data-selected candidates is added.

**Remaining boundary:** restoring decoys is necessary for faithful reuse, but
is not evidence for cross-run transfer. Reference/window/filter compatibility
and provenance are not enforced by this serializer. Reuse for a different run
or analysis configuration still needs independent justification. Existing
panels cannot recover missing decoys; refit matched controls when needed.

### 3. Empty-library and invalid-depth guard — audit F09

Batch preparation raises an informative error naming any sample or control with
no usable reads after preparation. It no longer allows that library to set the
entire batch's common depth to zero. Correcting or explicitly excluding that
library is required before rerunning. Nonpositive explicit or configured batch
depths are rejected. Valid positive-depth behavior is unchanged.

### 4. Restrict UMI indels — audit F03 mitigation

UMI consensus currently omits insertion voting and fragments multi-base
deletions. This pass does not attempt an unvalidated reconstruction algorithm.
Consensus construction rejects `exclude_indels=False`, and the indel caller
also rejects UMI mode to protect direct callers. The error explains how to run
substitution analysis. The restriction is documented in Params and the
consensus docstring. UMI substitution analysis and non-UMI indel analysis remain
available. This is an explicit feature restriction, not a repaired indel model.

## Validation

Validation results are recorded in the adjacent audit implementation package.
The new regression module covers:

- Candidate p-value/count parity at two fixed multipliers across singleton
  floors, ancestral mode, recalibration, minimum support, partial coverage,
  strand orientation, analysis windows, quality gates and end trimming.
- Exact panel round trips, including recalibration and an empty decoy list;
  version 1/2 compatibility and visible fallback warnings.
- Empty cases/controls at automatic and explicit depth, and invalid depths.
- Explicit rejection of unsupported UMI-indel calls.
- Actual spike-in harness refitting for different singleton-floor profiles.

These tests establish implementation invariants, not empirical FDR performance.
The manuscript benchmarks have not been regenerated. The pre-existing external
rRNA UMI benchmark wrapper is handled separately from this regression suite;
its earlier isolated run reproduced the technical error reduction.

## Deferred work and manuscript implications

The different zero-count conventions for replicate tests, sparse-count
intervals, incomplete-coverage denominators, indel class exposure, mixed-run
warnings, standalone/legacy CLI differences and clone/parallel issues remain
open. The empirical-decoy procedure is not promoted to formal FDR control.

Do not replace paper numbers based on this patch alone. In particular, rerun
quality-floor and ancestral-reference analyses affected by changed calibration.
Default physiological candidate selection is intended to remain unchanged;
passing regression tests does not prove all experimental outputs are identical.

## Nondestructive application and reversal

The changes were staged in a separate complete copy of the tracked checkout.
The implementation package holds the baseline manifest, a unified patch, tests,
validation logs, and original copies of every overwritten file. The guarded
application script checks that destination files match the recorded baseline
before writing. Its reversal operation checks that edited files still match
this pass before restoring originals; it refuses to overwrite later edits.
No commit, reset, build, manuscript replacement or figure regeneration is part
of this pass. The precise commands and application status are in the package's
`IMPLEMENTATION_NOTES.md`.

## Second pass: statistical conventions, sparse intervals, legacy route, exposures

This pass addresses audit findings F04, F05, F06 and F08. It also reports the
indel fallback from control mode to Benjamini–Hochberg (BH). Three choices were
made explicitly before implementation: one event as the primary zero-count
convention, profile-likelihood intervals, and labeling (rather than rerouting)
the legacy comparison command. Nominal alpha, presets and calling are unchanged.

### 5. One zero-count convention (F05)

`replicates.ZERO_EVENTS` is the single rule for every log-scale test, fold and
dot plot. A library with no events is placed at one event on its own exposure.
`project.group_table` (all metrics), `animal_tests` (`p_welch_log`,
`fold_geometric`, `fold_ci`, `p_permutation`), the power simulation and the
project dot plots now use it. Previously, project frequency and other
non-density rows used one third of the smallest positive value, and
`animal_tests` used half an event. Both survive only as labeled sensitivity
analyses: `*_sensitivity_p_welch_log_half_event_vs_*` and
`*_sensitivity_p_welch_log_third_of_smallest_vs_*` in the group table, and
`p_welch_log_half_event` and `p_welch_log_third_of_smallest` in `animal_tests`.
On the audit's example (controls 0/0/1, cases 2/3/4, 1,000 bases each), all
four primary outputs now give p = 0.0342 for the comparison where the reviewed
version gave 0.034, 0.022 and 0.008.

`fold_is_lower_bound` is renamed `fold_control_zero_one_event`. It marks a
zero-count convention, not a confidence bound. The web screen marks such folds
with † instead of ≥. The 18 October referee-response decision to keep
`_welch_log`'s third-of-smallest rule is superseded by this convention.

### 6. Sparse-count single-read interval (F04)

`singletons.estimate` inverts a profile likelihood for sample ~
Poisson(signal + background) and controls ~ Poisson(τ × background), scaled by
the control overdispersion. With no events anywhere the interval no longer
collapses to 0 [0, 0]. Coverage results and limitations are in
`docs/calibration/SINGLE_READ_INTERVAL.md`:

- 0.925–0.988 with Poisson background; the old interval fell to 0.70 and was
  zero-width in up to 46% of experiments.
- 0.87–0.90 when three heterogeneous controls carry several background events.
  This shortfall is not corrected.

`combined_density_ci_scope` states that the combined interval holds called
mutations fixed.

### 7. Legacy comparison command (F06)

`shmtrace-compare` is labeled a legacy, pooled, descriptive route in its
output, CLI help, stderr, plots and README. Each report lists the error model
every file was actually called with and surfaces any requested-model fallback;
on the JH4 example, all four files fell back from `qcal` to `read_quality`.
P values are headed "pooled p". The calibrated route is `shmtrace-sheet run` or
an app project.

### 8. Exposures (F08)

`call_variants` reports `analysed_bases`, the exact depth summed over callable
positions (window ∩ end trims). Every density uses that value through
`metrics.density_exposure`, instead of positions × rounded mean depth.
End-trimmed positions no longer enter `analysed_positions`, `mean_depth` or the
density denominator. The control-threshold base total is unchanged so that it
stays matched to the decoys' total. Project rows add
`mutation_frequency_observed_x1e3` and `observed_bases` beside the historical
fixed-window `mutation_frequency_x1e3`.

### 9. Indel multiple testing

Indel calls in control mode use BH, because indel decoys do not exist. This is
now reported as `indel_multiple_testing: "bh"` with a sample warning.

### Effect on the real JH4 and N114 projects

The JH4 and N114 projects (recommended mode, 500 reads, delivered reads) were
run before and after this pass. The JH4 and N114 windows lie inside their end
trims, so densities moved only by mean-depth rounding (≤ 0.01%). Default call
sets are unchanged. Projects on `generic_amplicon` (window `null`, 20-bp trims)
will see larger density changes. These have not been rerun.

| Output | JH4 before → after | N114 before → after |
|---|---|---|
| unique density, project log-Welch p | 0.00267 → 0.00267 | 0.00707 → 0.00707 |
| mutation frequency, project log-Welch p | 0.0052 → **0.0143** | 0.00064 → **0.00161** |
| `animal_tests` unique density p | 0.00117 → **0.00267** | 0.00096 → **0.00707** |
| `animal_tests` unique density geometric fold | 225 → **141** | 136 → **85** |
| `animal_tests` mutation frequency p | 0.0225 → **0.0143** | 0.00067 → **0.00161** |

Every primary p value remains below 0.05. Before replacing any manuscript
number, identify which of these fields supplied it. A value taken from the
project unique-density row (Figures 2 and 4 use one-event censoring) is
unchanged. A value taken from `animal_tests` or a project frequency row has
moved. Single-read intervals widened slightly. In JH4, WT4's single-read
interval now includes zero (upper limit −0.01 → 0.37). WT3's interval remains
entirely below zero, so the control background exceeds what WT3 shows. That
may reflect model misfit at very high mutation loads and deserves a separate
check.

The first-pass reversal script refuses to revert files that this pass edited.
`git diff` against `c3a633e` holds both passes. The first-pass state is
preserved in the audit implementation package's `revision/` directory.

## Third pass: preparation, quality gates, indel exposure, runs, clones, pickling

This pass addresses the remaining P2 findings F07 and F10–F15. Presets, nominal
alpha, the singleton floor default and the one-event zero-count convention are
unchanged. Regression tests are in `tests/test_audit_consistency.py`. Each new
test was also run against the pre-third-pass code, and each failed there for
the reason the audit gave.

### 10. Standalone null preparation (F07)

`null_model.fit` prepares each control file with `pipeline.prepare_reads`, the
sample route. The paralogue-allele filter, UMI consensus, `ignore_qualities` and
Sanger clean-window trimming now apply to controls. Controls keep their
distinct phasing policy (`control_params`: drop every block), and their blocks
are still handed back for the samples. The fitter takes the first reference
record, as `prepare_reads` does. It raises an error if asked for any other
record; no caller asked for one.

The test fits the same control FASTQs both ways: through `fit`, and through
`fit_aligned` on `prepare_reads` output. Rates, quality tails and site qualities
must agree, with residual paralogue reads, UMI-tagged reads and
`ignore_qualities`. Before the change, `fit` kept 48 reads where preparation
kept 40 (the paralogue reads) or 12 (UMI molecules).

**Remaining boundary:** the single-sample CLI still pools all `--null` controls
into one model. It has no per-case run check and no noisy-control handling.
README now says it is not identical to a project.

### 11. Quality gate in the null (F10)

`variants.gate_quality(p)` is the pileup gate: `min_base_quality`, or Q20 under
the Sanger threshold model. `fit_aligned` applies it to every count, depth, site
quality and reference-observation total. Every fitted probability, including
each cumulative quality tail above the gate, is therefore conditional on a base
passing the gate. `window_bases`, the decoys' base total, applies the same gate,
so it matches the sample's control-threshold total.

The test asserts that fitting with a Q20 gate equals fitting on the Q≥20 bases
alone, for rates, strand rates, tails and site depths. It also checks that the
Q20 tail equals the observed error share among eligible bases. In the probe,
low-quality bases make up 90% of the library and err at 20%; eligible bases err
at 2%. The ungated tail is more than five times smaller, and the decoy base
total matches `analysed_bases`. Recalibration is still fitted on all control
bases. Each quality bin is fitted separately, and bins below the gate are never
used.

**Effect:** none at the read-library default Q0. For `min_base_quality > 0`,
fitted tails rise and calls change. Rerun any non-default configuration.

### 12. Indel exposure and ancestral indels (F11)

Class opportunities are now summed over each pileup's window. At each position
and for each class (kind, length 1–4), the count is the reads that span an event
there (`Pileup.depth`). Previously it was every reference position times every
read. Every position of a repeat still counts, while events are left-aligned
onto one of them, so the normalisation at repeats is the same as before. Padding
the reference outside an unchanged window no longer alters class rates or
p values; a test checks this, and also checks the rates against a brute-force
count with unequal coverage.

`ancestral_reference` now applies to indels. The sample frequency filter is
skipped. A fixed event in controls is kept as that site's own rate, and it is
kept out of the class error count. Calibration tests fixed events as well. A
fixed insertion that the controls lack is called. If the controls carry it, the
fitted multiplier lifts its site rate to one, as for substitutions, and it is
not called.

**Effect:** indel class rates rise when the window is shorter than the reference
or coverage is incomplete. The indel spike-in in `docs/calibration/INDELS.md`
predates this and has not been rerun.

### 13. Per-case run check (F12)

`project.run_compatibility` judges each case on its own runs. It reads the
explicit `run` value (a string or a list), or else the run ID in each file's
first read header. A sample with several files is judged on every file, not the
first. Statuses are `same run`, `different run` and `unverified`. A case is
`unverified` when any run is unknown, and is never reported as matched.
`project.run_notes` writes a prominent "CALIBRATION NOT SUPPORTED for …" note,
adds a warning to each mismatched case, names controls that span several runs,
and lists unverified cases when any run is known. The project output carries
`run_check`.

**Deferred:** the project still fits one pooled null from all controls. It does
not select run-matched controls per case and does not refuse to run mismatched
cases. That would change outputs of existing projects and needs a decision
about the user interface. The note says to treat such cases as exploratory.

### 14. Clone summaries and the quality gate (F13)

`clones.call_clones` counts coverage and carriage only from A/C/G/T bases at or
above `gate_quality(p)`, inside the window. Sharing (Jaccard) uses the same
gate. The audit's two-read lineage now gives no clonal C at Q30. Without a gate
it gives the historical count. The majority rule, half of the covering reads
including ties, and equal weighting of each clone's majority profile are
unchanged and documented. That weighting is not the mean cellular burden.

### 15. Bounded clone pair sampling (F14)

When there are more than `max_pairs` pairs, `clones.sharing` draws `max_pairs`
pair indices uniformly without replacement. The draw uses
`random.Random(0).sample(range(total), k)` and decodes each index to its pair,
so memory is bounded by `max_pairs`. Sampling is across all pairs, not stratified
within and between clones. At or below the cap every pair is used, exactly as
before. Above it, the fixed-seed draw differs from the old shuffle, so capped
`sharing` values change slightly. The test brute-forces a 60-read set exactly.
On 3,000 reads (4.5 million pairs), traced peak memory stays under 8 MB; the old
code peaked at 292 MB.

### 16. Pickled reads keep added attributes (F15)

`AlignedRead.__reduce__` returns attributes set after construction as pickle
state, in both the packed and fallback forms. One of them is UMI consensus
`members`. The test checks a pickle round trip and `parallel.pmap` at one and two
processes. UMI truth scoring (`lineage_sim.score_umi`) therefore sees the same
reads at any job count. No full UMI benchmark was rerun at several job counts.

### What was run

`.venv/bin/python -m pytest -q --deselect
tests/test_umi_karst.py::test_consensus_removes_substitution_errors_on_real_umi_reads`
gave 379 passed: the 367 earlier tests plus 12 new ones. `ruff check src tests`
passes. No manuscript benchmark, figure or real-data project was
rerun in this pass. Default read-library calls (Q0, no ancestral reference,
substitutions only, single job) are expected to be unchanged, except where
standalone `--null` controls contain paralogue reads, UMIs or use
`ignore_qualities`. Clone `sharing` values above the pair cap will also change.
Before replacing any number, rerun the outputs that depend on a non-default
quality gate, indels, standalone `--null` fits or clone sharing.

## Fourth pass: run-matched controls, single-read background, full reruns

*20 October 2026.*

### 17. Run-matched controls per case (F12, second half)

`project.run_matched_peers` selects each sample's controls before any null is
fitted. A case whose runs are all known is judged only against controls from
those runs, provided that set has at least `MIN_RUN_CONTROLS` (2) controls. A
control needs one same-run peer. Otherwise the pooled controls are used, as
before, and the existing "CALIBRATION NOT SUPPORTED" note still applies.

The selection controls everything keyed by `peers_of`: the substitution null,
the decoys, the single-read controls and the indel null. One-run projects get
identical controls, so JH4 and N114 do not change. Each sample row records
`null_matching` (`same run`, `run-matched` or `pooled`). Setting
`"match_controls_by_run": false` in the project restores pooling.

### 18. Single-read background at unavailable pairs (new finding F16)

This came from investigating WT3's single-read interval, which lay below zero.
The sample count left out pairs the per-site test had called, and pairs that the
sample carries in two or more high-quality reads, because no single-read event
can be seen there. The control background still counted errors at those pairs.
In heavily mutated samples this overstated the expected background. In WT3,
8 of the 23 control events fell on such pairs.

`singletons.unavailable_pairs` now excludes those pairs from both counts. On
JH4, WT3's expected background fell from 7.5 to 4.9 events, and WT2's from 7.7
to 7.0. The controls barely changed.

**Still unexplained:** WT3 observes 0 events against 4.9 expected, so its
interval stays below zero. WT4 is also low, at 3 against 7.6. Efficiency is the
same in every library (0.83–0.84). The likeliest cause is library-to-library
variation in the high-quality error rate, which three controls cannot measure.
This is the limit documented in `SINGLE_READ_INTERVAL.md`. Treat single-read
estimates for heavily mutated samples as approximate.

### Full reruns

All 23 analyses in `ARCHIVE.md` that use SHMtrace were rerun with this code.
The rerun scripts and logs are in `~/Nanoshm/audit_rerun_2026-10-20/`. Not
rerun:

- the FDR stress tests behind Figure S15, which are not in the archive list and
  whose default call sets are unchanged;
- runs of the external callers (ClairS, DeepSomatic, DeepVariant, Clair3);
- the PacBio switch analysis (`switch.py` is unchanged).

Four regeneration problems turned up during the reruns:

- **Old outputs:** `crm/out/` is not tracked by git, so the 25 outputs
  rewritten before a snapshot was taken have no saved previous version. The
  manuscript-level comparison below is unaffected, because every number was
  saved first.
- **Changed figures:** Figures 3, 4, 7, S3 and S5 changed visibly. The other
  matplotlib figures changed only in file metadata and were restored. No
  SciCanvas figure's data changed.
- **iCloud:** several Sanger chromatogram folders had been offloaded from the
  Mac and stalled the reruns until they were downloaded.
- **Word limit:** the article build reports 7,647 words against the 7,000
  limit. That overage predates this work: the committed text counts 7,605 words
  by the build's own rule. It appeared when the count started to include figure
  legends (18 October). The wording changes in this pass keep the word count
  unchanged. Cutting the remaining words is an editorial task for the authors.

### Manuscript numbers

All 704 values the documents draw from were recomputed. Sixteen changed; twelve
of those appear in the text:

| Value | Before | After | Cause |
|---|---|---|---|
| J_H4 spike-in false calls per library / FDR | 2.4 / 0.048 | 2.6 / 0.050 | F07: standalone controls now prepared like samples |
| J_H4 one-read recall | 30% | 31% | F07 |
| Maea cohort, Wt vs KO frequency p | 0.0011 | 3.0 × 10⁻⁴ | F05 one-event rule |
| Maea cohort, Maea-KO vs KO frequency p | 0.0046 | 0.0039 | F05 |
| J_H4 indel recall, ≥3 reads | 78% | 72% | F11 indel exposure |
| J_H4 homopolymer 1-bp deletion recall | 31% | 19% | F11 |
| UMI simulation, read-mode null false calls J_H4 / N114 | 1.00 / 0.00 | 0.67 / 2.67 | F07/F15 on the UMI route |
| UMI molecules analysed | 141 | 165 | same |
| UMI simulation, read-mode recall range | 91–97% | 90–96% | same |

Every sentence that uses these values still holds:

- UMI-consensus null false calls stay at zero.
- J_H4 homopolymer deletions remain the least-recovered indel class.
- Every p value keeps its side of 0.05.

The four changed values that are not printed are the TAM guide folds and counts.
The Mkln1 evidence analysis was first rerun on the wrong FASTQ folder, which
has no Mkln1 libraries. Rerun on the correct raw-read folder, its sample sheet
and values match the committed ones. The quality-floor reruns are summarised in
`CALLER_REVIEW.md` Part 4: with sample/decoy parity the floor does not buy
precision, so it stays off.

## Fifth pass: false calls under the default, and two stricter settings

*20 October 2026.* The default threshold targets a false-discovery proportion
of 0.05. A heavily mutated library therefore admits more false calls in
absolute terms: 13.7 per library in simulated polyclonal lineages, against
Geneious's 4.0, in exchange for 83% versus 57% single-read recall.

This pass added two opt-in modes. `strict` sets α to 0.01. `strict_max` keeps
α at 0.05 but counts the noisiest control's decoys instead of their mean.
`crm/strict_sensitivity.py` compares both with the default. Table S7 holds
every value, and a STAR Methods paragraph explains the trade-off; STAR Methods
do not count toward the word limit.

- **`strict` (α 0.01):**
  - Lineage false calls fell from 13.7 to 3.0 per library, and single-read
    recall fell from 83% to 38%.
  - On the planted J_H4 controls, recall collapsed from 81% to 17% at an FDR of
    0.061, which is higher than the default's. N114 behaved as intended (76%
    recall, FDR 0.005).
  - With three controls, the decoy counts are too coarse to resolve a 1%
    target.
- **`strict_max`:** small effect on the lineages (12.0 false calls per
  library); J_H4 spike-in recall fell to 55%.
- **Real projects:** both settings removed calls only from mutation-bearing
  libraries. J_H4 Wt vs *Aicda*⁻/⁻ and N114 pMX-*Aicda* vs empty vector stayed
  71–141-fold apart, with p ≤ 0.008 throughout.

**Decision:** the default is kept. The trade-off is stated in the STAR Methods,
and the stricter modes remain available to users.

Two more zero-count conventions were found in the manuscript builder and
aligned with the one-event rule (F05):

- The *Ldha* effect interval dropped mice with no calls. Its label already said
  n = 7 vs 7, but it had used 6. The estimate moves from 1.23-fold (0.46–3.30)
  to 0.55-fold (0.07–4.49). The interval still spans 1, so "no detectable
  difference" stands.
- The HMCES deletion test used log10(k + 1). Its p value moves from
  8.4 × 10⁻⁵ to 1.0 × 10⁻⁴.
