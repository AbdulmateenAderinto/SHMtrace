# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project uses [semantic versioning](https://semver.org/).

## [Unreleased]

### Added
- **Renamed nanoshm → SHMtrace.** The package is `shmtrace`, commands are `shmtrace`, `shmtrace-sheet`, `shmtrace-power`, ... and the app is SHMtrace.app. `import nanoshm` still works (an alias with a DeprecationWarning); the old `nanoshm-*` commands are gone.
- Pooled runs and screens (`demux.py`, MULTIPLEX.md): inline single or dual barcodes split before analysis from a sample sheet's `barcode`/`barcode2` columns (no misassigned reads in simulation at 3% or 8% read error); `shmtrace-screen` ranks conditions against the control with animal-level tests and Benjamini-Hochberg q-values, with a volcano plot.
- Class-switch junction analysis (`switch.py`, `shmtrace-switch`, SWITCH.md): breakpoints, microhomology, insertions, mutations near the join, ambiguity in repetitive acceptors, internal Sμ deletions, unique junctions, animal-level group comparison. On simulated junction reads at 2% error: all found, breakpoints within MH + 3 bp for 95%, ≥4-nt microhomology classed correctly 97%.
- UMI consensus (`umi.py`; `umi = flank | header`): UMIs read next to an anchor at either or both read ends, or from read names; grouped with error tolerance (directional); one quality-weighted consensus per molecule. Simulated UMI libraries from real reads (`lineage_sim.build_umi_library`, `docs/calibration/umi_sim.py`, UMI.md): 0–0.7 false calls and 0.99–1.00 recall of single-molecule mutations, against 1–7 and 0.77–0.94 read by read.
- Group tests with the animal as the unit (`replicates.py`, `animal_tests` in every result, `animal_tests.csv`, app table): negative binomial with shared dispersion (bootstrap p), permutation, Welch on log10, fold with 95% CI; pooled-read p shown only as contrast. Calibrated by simulation: 5.0–5.7% false positives at α = 0.05, against 91–94% for pooled reads.
- `shmtrace-power` (`power.py`): power by simulation from a pilot project's rate, between-animal variation and within-animal clumping, for read-weighted, unique-density or clone-weighted measures.
- Repair-route indices and spectrum decomposition against reference groups (`repair.py`, `shmtrace-repair`, REPAIR.md).
- `docs/calibration/quant_accuracy.py` and QUANT_ACCURACY.md: the accuracy of the reported SHM numbers against truth on J꜀4 and N114, over 168 simulated libraries at 125–1,000 reads. `lineage_sim` now scores unique and combined density, takes a depth, and has an N114-like `subclonal` scenario.
- Sanger pipeline for unseen data (docs/calibration/SANGER_PIPELINE.md):
  - `shmtrace.samplesheet` / `shmtrace-sheet`: explicit file→sample assignment with order, well, label and exclusions. Every problem is reported by line; copied traces and re-run wells are warned about. Read lists name each clone by its order and well.
  - `shmtrace.qc`, run on every Sanger project: duplicate files and clones, mutations shared across animals beyond chance, control clones with more mutations than their qualities predict, label mismatches, mutations lost to the germline rule, load outliers, and every clone not analysed with its reason. Flags only; nothing is changed. The flags appear in the app (Tables), the HTML report, `qc_flags.csv` and `excluded.csv`.
  - `shmtrace.concordance` / `shmtrace-concordance`: the manuscript's two-platform statistics as a tested function; `manuscript_stats.py` now uses it, with unchanged numbers.
  - `shmtrace.sanger_sim` / `shmtrace-validate-sanger`: simulated Sanger projects with known truth and planted bookkeeping errors. Over 40 projects: recall 99.4%, no germline calls, every planted copy, contaminant and labelled misfiling flagged, and no false flags on 6,018 clean clones.
  - The app loads a sample sheet by path. `io_utils.write_abif`. The germline rule now records which positions it removed (`germline_sites`).
  - Tests: `test_samplesheet.py`, `test_concordance.py`, `test_sanger_validation.py`, and `test_real_jh4_sanger.py`, which reproduces the manuscript's Sanger counts from the lab's traces when they are present.
- `organize_jh4_ab1.py --original`: organise the .ab1 traces exactly as filed in the GenBank folders (no evidence-based moves); `manuscript_corrected/export_sanger_calls.py` exports the Sanger calls the manuscript figures read. Ver17 now uses the original classification, drops the TOST analysis and the inclusion figure; `--original` still drops the two Wt2 GenBank records that are renamed MKLN K4 traces; `sanger_window_figure.py --name`; `sanger_vs_nanopore.py` also writes the two-panel `fig3_agreement_position_hotspot`.
- **Sanger against nanopore, same mice** (`docs/calibration/sanger_vs_nanopore.py`,
  `SANGER_VS_NANOPORE.md`). Redraws the manuscript's Figure 2A/B, 3A–E and
  S2/S3 from SHMtrace calls on both platforms, plus position, hotspot and
  per-mouse load agreement. Spectrum r = 0.92 (0.94 under the Geneious rule),
  no platform × type interaction (p = 0.39), positional r = 0.90. It also
  found that the manuscript's Geneious Sanger assemblies mix mice (W2 and W3
  are each two mice) and that the Aicda-/- Sanger denominators are 54% JH6
  clones. Checked from sequence alone (mutation load and shared lineages),
  the misfiled wild-type clones belong to the mouse their names give, but
  three clones filed as WT2 and WT4 are MKLN clones (lineages of 8 and 7
  identical mutations shared with K4 and K2) and are removed.
- **The J꜀4 Sanger traces organised by mouse** (`docs/calibration/organize_jh4_ab1.py`).
  288 .ab1 traces copied into one folder per mouse, renamed
  `<mouse>_<order>-<well>_<label>.ab1`, with a manifest. Labels the sequence
  contradicts are corrected: two WT3 clones labelled w2 (18 and 13 mutations,
  in the w3 wells of their plate), w4-9 to MKLN K2, and an AID-type clone in an
  Aicda-/- sample left unassigned. From the traces, WT2 is nearly unmutated
  (30 of 34 clones), matching its nanopore library's dominant unmutated clone,
  and Sanger alone no longer separates wild type from Aicda-/- (P = 0.27);
  nanopore does (P = 0.0017). The comparison now reads these folders: spectrum
  r = 0.90, positional r = 0.89.
  With the January 2024 orders added (429 traces), 214 of the 220 wild-type
  and Aicda-/- clones are read from traces; the jhextra traces carry exactly
  the mutations their GenBank records do, and every result is unchanged.
- **Simulated libraries with known clones** (`shmtrace/lineage_sim.py`,
  `shmtrace-simulate-lineages`). Real control reads are dealt into clones with
  their own junctions and J trimming, lineage mutations (trunk, branch,
  sub-branch, private) placed by the 5-mer model, and optionally one expanded
  clone; the library runs end to end and is scored against the truth, with a
  threshold-sweep ROC. Clones are recovered at adjusted Rand 0.95–0.97; with
  one unmutated clone holding 60% of reads, clone-weighting recovers 85% of
  the true per-clone frequency against 46% read-weighted. It also found that
  Geneious's ordering edges the recommended model on recall in heavily mutated
  libraries, while making three times its false calls in an empty one; reported
  as open. Variant calling now records each candidate's position and allele
  (`_candidate_keys`) for such sweeps. See `docs/calibration/LINEAGE_SIM.md`.
- **Clones from the V(D)J junction** (`shmtrace/clones.py`; on by default for
  `jh4_intron`). The bases 5' of J꜀4 are each clone's recombination junction,
  so reads are grouped by it and each clone counts once. Validated as a
  barcode: reads sharing a junction share their mutations 26× (Wt3) and 324×
  (Wt4) more than reads that do not. It found that 85% of Wt2's junction reads
  are one almost unmutated clone, which held Wt2's read-weighted frequency at
  1.74 against 12.00 per clone, and clone-weighting strengthens wild type
  against *Aicda*-/- from P = 0.0086 to P = 0.0036. New columns, `clones.csv`,
  a clone figure, and a warning when one clone holds half a library. Per
  clone, not per molecule. See `docs/calibration/CLONES.md`.
- **AID targeting against a 5-mer model** (`shmtrace/targeting.py`; the project
  key `targeting_model`, or the model drop under Advanced). The targeting index
  is the mean model mutability at called positions over the window mean, with a
  permutation p-value. With SHazaM's HH_S5F, J꜀4 wild-type calls score
  1.19-1.34 (p ≤ 0.04) and knockout sequencing errors 1.00 (p = 0.47); on CH12
  real AID scores 2.55 and the R174E/R178E error cluster 0.72, so the index
  flags that artefact by itself. Each call is annotated with its 5-mer's
  relative mutability and substitution probability. Three figures: a
  dose-response by mutability fifth, observed against expected spectrum, and
  the model as a hedgehog plot. Nothing here changes the calls. SHMtrace does
  not ship SHazaM's AGPL models; `scripts/export_shazam_model.R` exports one
  from your own install. See `docs/calibration/TARGETING.md`.
- **Coverage along the amplicon.** Every sample's per-position read depth is
  kept (`coverage.csv`, a coverage figure), and a stretch of ten or more
  positions below half the sample's median depth is named in its warnings.
  Density divides by analysed positions times mean depth, so it assumes even
  coverage; this is where to see whether that holds.

### Fixed
- Manuscript Sanger set: Wt4 clone w4-9 (434141-A11) moved to MKLN K2 (identical over 520 bp, all 7 mutations, to two K2 clones), applied in `organize_jh4_ab1.py --original` too. Wt4 26 clones, 106 calls; Ver17 rebuilt (~39-fold, P = 0.0008; r = 0.96).
- Paralogue (Agouti-type) J꜀4 reads longer than the bundled paralogue record, and PCR chimeras, passed the contaminant filter and were called as mutations: 1–4% of every raw J꜀4 library, with 35–54 false calls per simulated 500-read library. The documented allele-level check is now implemented (`contaminant_min_alleles`, default 10; germline sites ignored). The manuscript's pre-sorted 500-read files contain none of these reads.
- The single-read estimate counted events in every prepared read, not the analysed ones, so a library deeper than the analysis depth overstated the combined density by 1.3–2×. It now counts at the analysed depth, controls included.
- **A Sanger sample of one to three clones reported germline differences as
  mutations.** Too few clones to apply the germline rule to itself, so it was
  switched off; one Aicda-/- clone reported the strain's five J꜀4 germline
  differences as five mutations. Such a sample now uses the germline changes
  measured across the project's other clones.
- **Clone clustering charged one junction indel as two errors.** The
  junction window's 5' edge is an arbitrary cut, so an indel inside it slides
  a base in at that edge; the distance now leaves that edge free. Found by the
  lineage simulation (adjusted Rand 0.96 → 0.99 at ten times the real indel
  rate). On the real libraries it merges a few splinter clones; `CLONES.md`
  carries the re-run numbers.
- **The CHANGELOG had entries repeated up to five times**, from edits that
  matched every `### Added` heading rather than the first; rebuilt with each
  entry once.
- **A sample warning written before the sample's row existed.** The coverage
  warning added earlier in this release referenced the row too early and would
  have raised an error on the first library with a real coverage dip; the test
  data had none, so it never fired. Now covered by a regression test.
- **Base-composition correction in projects.** `shmtrace-compare` already
  divided each substitution class by the sites of its from-base (SHMTool's
  L = sites × files); the project path reported raw counts only, so the two
  entry points gave differently normalised spectra. Projects now report
  `spectrum_corrected_group` and `_sample`, the window's `base_composition`, a
  `per_site_per_sample` column in `spectrum.csv`, and a corrected-spectrum
  figure. On the J꜀4 window, which holds 82 C against 146 G, raw counts
  understate C→T relative to G→A.

### Added (continued)
- **A read-alignment browser** (`shmtrace/browse.py`, the app's **Reads** tab,
  opt-in with "Keep the read alignment"). Every input read is drawn in its
  aligned place against the reference, banded by what became of it — analysed,
  not sampled, minor haplotype, too short, chimeric, contaminant, off target,
  low quality, too long, unmapped — so a dropped read is still visible where it
  would have sat. Zoom to base level for the reference letters and each
  mismatched base, with called positions marked and indels ticked; click a read
  for its strand, span and every mismatch with its quality. Reads travel as
  their differences from the reference rather than their sequence, which is
  what makes a few thousand of them cheap enough to send to a browser, and the
  API serves one window at a time. Every filtering stage now records why it
  rejected a read.
- **Twelve figures beyond the manuscript's seven** (`shmtrace/viz_diagnostics.py`),
  in three groups: detail the project already computed but never drew (the
  saturation curve, the single-read excess by substitution class, indels,
  per-sample positional maps, the combined density with its intervals);
  diagnostics that make the method auditable rather than taken on trust (the
  base-quality histogram the default model keys on, the read-by-site carriage
  behind each phasing decision, read attrition per stage, the fitted control
  background along the amplicon, and the caller's own decision boundary); and
  biology the tables held but nothing showed (clonal structure, replicate
  agreement). Each is drawn only when the run carries what it needs. The
  saturation curve is opt-in, since it re-calls variants at a ladder of depths.
  To support them the caller now records every tested candidate, the base
  qualities of matching and mismatching bases, and the per-read carriage of
  each phasing block; the project records the fitted control background.

### Fixed
- **A calibration bug introduced and corrected in the same session.** Capping
  the quality-calibrated test for speed was not neutral: the same function
  drives the leave-one-out calibration, where more conservative p-values let it
  settle on a *smaller* inflation factor and so a *more permissive* null. On
  deep data it chose 3.0 where the exact computation gives 20.0, roughly
  doubling the calls. The cap is gone; the exact test instead evaluates its
  Poisson tails in one vectorised call rather than summing each term by hand,
  which also removes a catastrophic cancellation in the deep tail. Measured per
  candidate: 3,073x faster at 8,967 supporting reads and 104x at 1,000, but
  about 50x slower at 10, where scipy's call overhead dominates. In practice
  the nanopore J꜀4 project is unchanged at 26 s, because read preparation and
  alignment dominate there, and the Illumina J꜀4 numbers are restored exactly
  to their pre-cap values.

### Added
- **Single-read mutations, estimated in aggregate** (`shmtrace/singletons.py`).
  A mutation in one read cannot be confirmed site by site, so each sample's
  high-quality single-read events are counted and the controls' rate of such
  events is subtracted, corrected for the share of real bases the quality cut
  removes. Projects with two or more controls now report
  `single_read_estimate` and `combined_density_x1e5` (per-site calls plus the
  estimate) with 95% intervals. On J꜀4 spike-ins it recovers 91–104% of planted
  one-read mutations, and it brings Wt2 from 55.7 to 74.9 × 10⁻⁵ while the
  knockouts stay at 0.4–1.2. See `docs/calibration/SINGLE_READ.md`.
  `single_read_min_quality` (35; 30 on Illumina; off for Sanger) sets the cut.
- **Indel calling** (`shmtrace/indels.py`, opt-in with `exclude_indels: false` or
  the app's Advanced checkbox). Indels are left-aligned, counted at quality
  ≥ Q30, and tested against the controls' rate for that exact indel, calibrated
  by leave-one-control-out and corrected by BH; at least 3 reads are needed.
  Reported in `indels.csv` and per-sample counts, not in the substitution
  density. Spike-ins on J꜀4: recall 0.43 at 3 reads, 0.72 at 5, 0.94 at 10;
  1-bp homopolymer deletions stay mostly undetectable; 0.71 false calls per
  library. Sanger clones call every indel at ≥ Q20. See
  `docs/calibration/INDELS.md`.
- **GenBank input** (`.gb`/`.gbk`), as Benchling and sequencing providers
  export Sanger clones: one read per record, flat Q30, ambiguity codes kept.
- **A J꜀4 paralogue in the bundled contaminant reference.** 18% of the lab's J꜀4
  Sanger clones, and up to 41% of raw nanopore J꜀4 reads, are a J꜀4-like
  sequence about 11% divergent from the intron that the J꜀4 primers also
  amplify. Its consensus is now the second record of
  `mouse_jh4_contaminant.fasta`. Left in, it inflated one raw library's
  mutation frequency fivefold. See `docs/calibration/SANGER.md`, which also
  compares Sanger clones with nanopore in the same animals (spectrum r = 0.91).
- **A positional-clustering flag.** Calls packed into a 10-bp window far more
  tightly than chance allows are counted in `clustered_calls` and named in the
  sample's warnings. Real SHM is spread out, and a genuine co-occurring block
  shares reads and is handled by phasing; a localised PCR or sequencing
  artefact is neither, and a control-fitted null cannot see one that is absent
  from the controls. Found on PRJNA918596 (`docs/calibration/ECTOPIC_AID.md`).
- **`pct_CT_GA`**, C→T and G→A as a share of calls: the signature AID leaves
  when the uracil goes unrepaired, and a quick check that calls look like AID's.
- The app's Tables tab shows every column any row has, so a control without
  peers no longer hides the other samples' estimates.

### Fixed (Mac app, found by running the built app end to end)
- SVG and PDF figure exports failed in the packaged app: matplotlib's PDF and
  SVG backends were not bundled. PNG was unaffected.
- Downloads did nothing inside the app window: pywebview blocks them unless
  `ALLOW_DOWNLOADS` is set.
- "Save to a folder" could stay hidden, because the desktop flag was injected
  after the page had read it. The window now opens with `?desktop=1`.
- `shmtrace-app --port N` with N already in use reported success and opened
  whatever else was listening there. It now says the port is taken and exits.
- The packaged app rebuilt matplotlib's font cache on every launch (about 12 s);
  it now keeps it in the user cache folder.

## Unreleased — the standalone app (11 October 2026)

### Added
- **Desktop app.** `shmtrace-app` opens a window with four tabs: set up a project,
  read the results, browse the tables, export. Everything runs locally.
  `packaging/shmtrace.spec` builds a `.app` or `.exe`; see `packaging/BUILD.md`.
- **Projects.** `shmtrace.project.run_project` analyses many samples at once with
  groups and controls, one shared depth, one background calibration, and
  per-sample and per-group tables including Welch and Mann-Whitney tests.
- **Interactive charts** (Plotly) in the app and static figures for export:
  density, mutation frequency, % of reads mutated, substitution spectrum,
  position along the amplicon, mutations per read, hotspot targeting.
- **Exports.** `shmtrace.report` writes CSV, Excel, PNG/SVG/PDF figures, a
  self-contained HTML report, and a ZIP of all of it.
- **More platforms.** Sanger traces (.ab1), Illumina (including paired-end
  merging), PacBio and nanopore, selected with `--platform`.
- **More regions.** Presets for the Sμ switch region, any amplicon
  (`generic_amplicon`) and Sanger clones, plus bundled reference sequences.
- **External validation** on BioProject PRJNA1481264 and the GHVPDN run:
  `docs/calibration/EXTERNAL_VALIDATION.md`.

### Changed (11 October 2026)

- `--mode recommended` is now the quality-calibrated model (`qcal`), chosen on the JH4 and N114 benchmark and spike-ins (`docs/calibration/ERROR_MODELS.md`). The previous default is `--mode v0.1`.
- `shmtrace-benchmark` fits each control set's null once and runs the qcal calibration only for qcal modes.

### Added (11 October 2026, second round)
- **`--mode qcal`, quality-calibrated null.** Geneious trusts a high-quality
  variant base; SHMtrace now measures how far to trust it. From the control
  reads it learns, per base change and strand, how often a base is miscalled
  with at least a given quality, plus each site's own miscalls, and inflates
  those rates by a factor fitted by testing each control against the others.
  See `docs/calibration/ERROR_MODELS.md`.
- **`--mode strand`, strand-aware beta-binomial null.** Each strand is tested
  against its own control error rate with between-library overdispersion;
  strand p-values are combined by Fisher's method (`strand_combine="both"`
  requires each strand).
- **`--strand-filter bias|both`** (or `+filter` on any mode) and per-call strand
  counts and strand-bias p in `variants.csv`.
- **Homopolymer flag** (`homopolymer` column; `--mask-homopolymers`) and the
  reference 5-mer (`context`) for every call.
- **`error_model="pbinom"`**: exact Poisson-binomial over read qualities
  (LoFreq's test).
- **Saved background panels**: `--save-null-panel` / `--null-panel` store and
  reuse a fitted null (rates, strand rates, quality tails, calibration, provenance).
- **`shmtrace-validate`**: plants AID-like mutations at known read counts in
  real control reads and reports recall by support level and false calls for
  each mode.
- **BAM/SAM input** (`pip install SHMtrace[bam]`): reads straight from a Geneious
  assembly, restored to sequenced orientation.

### Fixed (11 October 2026, second round)
- `shmtrace-compare` pooled unique mutations across files into one set.
  SHMTool counts a mutation once per file, with L = sites x files; counts are
  now summed per file. Added non-unique (every read) tests with L = sites x
  reads and a per-class 2x2 test for all 12 substitutions.

### Added (October 2026)
- **`--mode geneious`**: Geneious's approximate variant p-value (binomial tail at
  the mean quality of the variant bases), uncorrected, no phasing. Reproduces
  the manuscript's Table S1 call counts within 0-11 calls per library. See
  `docs/calibration/GENEIOUS.md`.
- **`shmtrace-benchmark`**: runs a sample sheet in every mode at one common
  depth, with leave-one-out empirical nulls from the same run, and compares
  against Geneious call counts or (with exports) site by site.
- **`--rarefaction`**: unique variants and density at several depths with
  repeated draws, since density saturates with depth.
- `shmtrace.batch`: shared haplotype evidence and a common depth for multi-sample
  runs; used by `shmtrace-compare` and the benchmark.
- Calling mode and null-control upload in the web UI.
- `--run-id` / `--null-run-id`; run identity is also read from the MinKNOW
  `runid=` header tag.

### Changed (October 2026)
- **Phasing no longer drops every co-occurring block.** An expanded clone looks
  like a haplotype block; on the bundled wild-type data the old behaviour
  removed 19% of Wt3 reads on a mouse-specific block. New default
  `phase_mode="auto"` drops a block only if it is also seen in mutation-free
  controls or recurs in another sample, and flags the rest.
  `--phase-mode drop` restores the old behaviour. Bundled WT3: 281 -> 310 calls.
- `shmtrace-compare` analyses every file at the same depth.
- A sample below its target depth is now reported in `summary.json` warnings.
- `summary.json` records `shmtrace_version`.

### Fixed (October 2026)
- `minimap2_preset` was ignored (hardcoded `map-ont`).

### Added
- **Minor-haplotype phasing** (`phasing.py`, on by default). A second haplotype,
  paralogue or surviving contaminant appears as a block of positions carried by
  the same subset of reads. Detected by co-occurrence rather than frequency,
  because the block sat at 2% of reads in one JH4 mouse and 17% in another.
  Refuses to act if the population exceeds `phase_max_drop`. Disable with
  `--no-phasing`.
- **Benjamini-Hochberg correction** across all candidates, on by default. The
  caller runs ~1,587 tests per sample; uncorrected p <= 0.05 expected ~79 false
  calls, comparable to the entire observed background. `--fdr` sets alpha,
  `--no-correction` restores the old behaviour.
- **`--null`**: measure the per-position, per-alt error rate from mutation-free
  control reads and test candidates against it. The quality-derived null assumes
  independence across reads and understated real error ~3-fold.

- **`null_rate_floor`** (default on): the empirical error rate is never allowed
  below the read-quality expectation for that column. Fitting a null on
  exceptionally clean controls otherwise drives the rate below the sequencer's
  own implied error and makes the test *more* permissive, which is how a
  cross-run null added calls instead of removing them.
- **Run provenance check on `--null`.** Read IDs encode the sequencing run;
  if the controls come from a different run the tool warns and falls back to
  the read-quality model unless `--force-null` is passed. Null provenance is
  recorded in `summary.json`.

### Changed
- Depth downsampling now happens *after* alignment and phasing, so
  `downsample_to` is the number of reads actually analysed rather than the
  number submitted to the aligner.
- `min_variant_reads` stays at 1. Once the statistics are right, raising it buys
  no background and costs ~15% of signal.

### Fixed
- Held-out *Aicda*-/- background: 30.9 -> 1.0 x 10-5 per analysed base on
  unphased input, with wild-type signal retained. Separation 2.9x -> 76x.
  Manuscript reference value is 1.2 x 10-5. See `docs/calibration/FINDINGS.md`.

## [0.1.0] - 2026-07-31

First packaged release. Extracted from the analysis scripts written alongside
the manuscript and turned into an installable project.

### Added
- `shmtrace` CLI: FASTQ + reference (+ contaminant) to variant table, unique
  variant density and substitution spectrum.
- `shmtrace-compare` CLI: case-versus-control comparison in the spirit of
  SHMTool (MacCarthy et al., 2009), starting from raw FASTQ.
- Competitive minimap2 alignment via `mappy`, expected-error read filter,
  amplicon-length gate, deterministic depth downsampling.
- Quality-aware SNV caller reproducing the manuscript's filtering rules.
- Bundled presets for the JH4 intron and the N114 V region.
- FastAPI web UI with single-sample and comparison tabs.
- Publication-style figures: spectrum, spatial map, mutation-load histogram.
- Test suite and continuous integration.
