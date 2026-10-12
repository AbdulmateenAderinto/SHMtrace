# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project uses [semantic versioning](https://semver.org/).

## [Unreleased]

### Added
- **Single-read mutations, estimated in aggregate** (`nanoshm/singletons.py`).
  A mutation in one read cannot be confirmed site by site, so each sample's
  high-quality single-read events are counted and the controls' rate of such
  events is subtracted, corrected for the share of real bases the quality cut
  removes. Projects with two or more controls now report
  `single_read_estimate` and `combined_density_x1e5` (per-site calls plus the
  estimate) with 95% intervals. On J꜀4 spike-ins it recovers 91–104% of planted
  one-read mutations, and it brings Wt2 from 55.7 to 74.9 × 10⁻⁵ while the
  knockouts stay at 0.4–1.2. See `docs/calibration/SINGLE_READ.md`.
  `single_read_min_quality` (35; 30 on Illumina; off for Sanger) sets the cut.
- **Indel calling** (`nanoshm/indels.py`, opt-in with `exclude_indels: false` or
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
- `nanoshm-app --port N` with N already in use reported success and opened
  whatever else was listening there. It now says the port is taken and exits.
- The packaged app rebuilt matplotlib's font cache on every launch (about 12 s);
  it now keeps it in the user cache folder.

## Unreleased — the standalone app (11 October 2026)

### Added
- **Desktop app.** `nanoshm-app` opens a window with four tabs: set up a project,
  read the results, browse the tables, export. Everything runs locally.
  `packaging/nanoshm.spec` builds a `.app` or `.exe`; see `packaging/BUILD.md`.
- **Projects.** `nanoshm.project.run_project` analyses many samples at once with
  groups and controls, one shared depth, one background calibration, and
  per-sample and per-group tables including Welch and Mann-Whitney tests.
- **Interactive charts** (Plotly) in the app and static figures for export:
  density, mutation frequency, % of reads mutated, substitution spectrum,
  position along the amplicon, mutations per read, hotspot targeting.
- **Exports.** `nanoshm.report` writes CSV, Excel, PNG/SVG/PDF figures, a
  self-contained HTML report, and a ZIP of all of it.
- **More platforms.** Sanger traces (.ab1), Illumina (including paired-end
  merging), PacBio and nanopore, selected with `--platform`.
- **More regions.** Presets for the Sμ switch region, any amplicon
  (`generic_amplicon`) and Sanger clones, plus bundled reference sequences.
- **External validation** on BioProject PRJNA1481264 and the GHVPDN run:
  `docs/calibration/EXTERNAL_VALIDATION.md`.

### Changed (11 October 2026)

- `--mode recommended` is now the quality-calibrated model (`qcal`), chosen on the JH4 and N114 benchmark and spike-ins (`docs/calibration/ERROR_MODELS.md`). The previous default is `--mode v0.1`.
- `nanoshm-benchmark` fits each control set's null once and runs the qcal calibration only for qcal modes.

### Added (11 October 2026, second round)
- **`--mode qcal`, quality-calibrated null.** Geneious trusts a high-quality
  variant base; nanoshm now measures how far to trust it. From the control
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
- **`nanoshm-validate`**: plants AID-like mutations at known read counts in
  real control reads and reports recall by support level and false calls for
  each mode.
- **BAM/SAM input** (`pip install nanoshm[bam]`): reads straight from a Geneious
  assembly, restored to sequenced orientation.

### Fixed (11 October 2026, second round)
- `nanoshm-compare` pooled unique mutations across files into one set.
  SHMTool counts a mutation once per file, with L = sites x files; counts are
  now summed per file. Added non-unique (every read) tests with L = sites x
  reads and a per-class 2x2 test for all 12 substitutions.

### Added (October 2026)
- **`--mode geneious`**: Geneious's approximate variant p-value (binomial tail at
  the mean quality of the variant bases), uncorrected, no phasing. Reproduces
  the manuscript's Table S1 call counts within 0-11 calls per library. See
  `docs/calibration/GENEIOUS.md`.
- **`nanoshm-benchmark`**: runs a sample sheet in every mode at one common
  depth, with leave-one-out empirical nulls from the same run, and compares
  against Geneious call counts or (with exports) site by site.
- **`--rarefaction`**: unique variants and density at several depths with
  repeated draws, since density saturates with depth.
- `nanoshm.batch`: shared haplotype evidence and a common depth for multi-sample
  runs; used by `nanoshm-compare` and the benchmark.
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
- `nanoshm-compare` analyses every file at the same depth.
- A sample below its target depth is now reported in `summary.json` warnings.
- `summary.json` records `nanoshm_version`.

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
- `nanoshm` CLI: FASTQ + reference (+ contaminant) to variant table, unique
  variant density and substitution spectrum.
- `nanoshm-compare` CLI: case-versus-control comparison in the spirit of
  SHMTool (MacCarthy et al., 2009), starting from raw FASTQ.
- Competitive minimap2 alignment via `mappy`, expected-error read filter,
  amplicon-length gate, deterministic depth downsampling.
- Quality-aware SNV caller reproducing the manuscript's filtering rules.
- Bundled presets for the JH4 intron and the N114 V region.
- FastAPI web UI with single-sample and comparison tabs.
- Publication-style figures: spectrum, spatial map, mutation-load histogram.
- Test suite and continuous integration.
