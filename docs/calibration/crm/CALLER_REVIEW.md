# How SHMtrace calls a mutation, and what the Figure 3D–E comparison actually measures

Written 2026-10-19 for a second opinion. Part 1 describes the caller as it stands,
with file and line references, so a reader who has not seen the code can judge the
rest. Part 2 states the claim under review. Part 3 is what was measured and how.
Part 4 is the one change made to the software. Part 5 lists what I am least sure
about.

Nothing in Parts 1–3 changes any published number. The change in Part 4 is off by
default and no manuscript number uses it.

---

## Part 1. How SHMtrace works today

### 1.1 The problem it is built for

Somatic hypermutation is measured by counting distinct substitutions in an
amplicon. On nanopore reads the obstacle is that the per-base error rate and the
per-base mutation rate are the same order of magnitude. A mutation carried by
many reads is easy; a mutation carried by one read is, at a single position,
indistinguishable from a miscall. The software's whole design is about where to
draw that line and how to justify where it was drawn.

Two things make the problem tractable:

1. **Matched AID-null controls sequenced on the same flow cell.** An
   *Aicda*<sup>−/−</sup> library lacks AID-dependent SHM, so after the stated reference
   and preparation filters its mismatches measure technical and AID-independent
   background. That gives a measured null rather than an assumed one.
2. **The unit of interest is the distinct substitution, not the read.** A
   mutation is counted once however many reads carry it.

### 1.2 Pipeline order

`pipeline.prepare_reads` (`src/shmtrace/pipeline.py:38`) then `pipeline.analyse`
(`:100`):

| Step | Module | What it does |
|---|---|---|
| 1 | `reads_filter.prefilter` | expected-error rate and length gate; optional Mott trim for Sanger |
| 2 | `align.build_index` / `align_reads` | minimap2 (via `mappy`) indexes target **and** contaminant references together; each read goes to its best-scoring hit; reads whose best hit is a contaminant, that align too short, or that are supplementary/chimeric are dropped |
| 3 | `align.drop_contaminant_alleles` | allele-level paralogue check: a read carrying ≥ `contaminant_min_alleles` of a contaminant's diagnostic bases is dropped even if it still maps best to the target |
| 4 | `umi.consensus_reads` | optional; one consensus read per molecule |
| 5 | `phasing.partition` | removes a co-occurring block of differences (second haplotype, paralogue) — see 1.3 |
| 6 | `pipeline.downsample` | all samples in a project drawn to one common depth, because distinct-variant density saturates with depth |
| 7 | `variants.call_variants` | the pileup and the test — see 1.4 onward |
| 8 | `clones.cluster` | groups reads by V(D)J junction so a clone counts once |
| 9 | `singletons` | the aggregate single-read estimate — see 1.8 |

Step 6 matters more than it looks. Whether a mutation is seen in "exactly one
read" depends on how many reads there are, so samples and controls must be at the
same depth for any single-read reasoning to mean anything.

### 1.3 Phasing (`phasing.py`)

Competitive alignment removes reads whose *best* hit is elsewhere, but not a
population that still maps best to the target — a second haplotype or a
paralogous locus. Those reads carry a coherent block of differences, and pooling
them makes every position in the block look like a low-frequency variant.

The discriminator is **co-occurrence, not frequency**, with a caveat: shared
mutations can mark either an expanded lineage or a non-target sequence population,
and other evidence is needed to tell them apart. A haplotype block carries every
site on the same reads. On the J<sub>H</sub>4 data, 14% of *Aicda*<sup>−/−</sup>
reads carried all sixteen worst-offending positions while 83% carried none.
Removing that population accounts for roughly a sevenfold drop in apparent
background.

### 1.4 The pileup

`variants.call_variants` (`src/shmtrace/variants.py:246` onward) walks every
aligned base once. Bases below `min_base_quality` (default 0, i.e. all of them)
are excluded entirely. For each reference position it accumulates per-allele
counts, per-strand counts and depths, and — for the quality-calibrated model —
the list of base qualities of every mismatching base, per strand
(`alt_quals[(pos, alt, strand)]`).

Then, for each (position, alt):

- skip if supporting reads < `min_variant_reads` (default 1);
- skip and record as **germline** if the alt reaches `germline_freq_threshold`
  (default 0.50) — the manuscript's "< 50%" rule;
- otherwise compute a p-value under the chosen error model.

Candidates within `end_trim_bp` of the reference ends, and outside
`analysis_window`, never enter.

### 1.5 The error models

Five, selected by `error_model`:

| name | test |
|---|---|
| `qcal` | **the default.** Quality-calibrated against control-measured miscall rates — see below |
| `geneious` | Geneious's approximate p-value: binomial tail at the *mean quality of the variant bases*, `P(X ≥ k), X ~ Bin(depth, 10^(−Q̄/10))`. The error is **not** divided by 3. This is the manuscript's comparator (`variants.py:104`) |
| `strand` | per-strand beta-binomial against control rates, strands combined |
| `pbinom` | Poisson-binomial over every read's quality in the column (LoFreq-like) |
| `sanger` | one read = one clone: every substitution at base Q ≥ `min_base_quality` is called, no test |

**The `qcal` test** (`variants.py:574`). Geneious's insight is that a
high-quality variant base is more believable than a low-quality one; its weakness
is that it trusts the stated quality. Here the trust is measured. For supporting
reads ranked by quality q(1) ≥ q(2) ≥ …, the j-th test asks how likely j or more
miscalls *with quality ≥ q(j)* are among the column's reads, under
`X ~ Poisson(Σ_strands depth × rate)`. The smallest of these p-values, times the
number of tests, is the call's p-value.

The rate comes from `NullModel.qrate` (`null_model.py:185`):

```
g    = qtail[(ref, alt, strand)][q]      # pooled over all sites: share of control
                                         # bases at `ref` miscalled as `alt` on this
                                         # strand with quality >= q
n    = site_depth[(pos, strand)]         # control bases observed at this site
k    = control miscalls at this site with quality >= q
site = (k + prior_trials*g) / (n + prior_trials)
rate = min(1, qscale * max(g, site))
```

Two things to note. First, the site estimate is **floored at the pooled rate**:
site-level evidence can raise the error rate but never lower it. Second,
everything is multiplied by `qscale`.

**`qscale`** (`null_model.calibrate_qscale`, `:308`) is a global inflation factor
fitted by leave-one-control-out: each control is tested against a null fitted
from the other controls, and the smallest factor in `(1, 2, 3, 5, 10, 20, 50)` at
which no more than 1% of its candidates reach p < 0.01 is used. On the
J<sub>H</sub>4 *Aicda*<sup>−/−</sup> libraries this lands at **10**. The purpose
is to make p-values honest across libraries, because nanopore error varies more
between libraries than within one.

### 1.6 Multiple testing, and where the threshold comes from

`multiple_testing` is one of `none`, `bh`, `control`. The default `recommended`
mode uses **`control`** (`pipeline.MODES`, `:285`):

```python
"recommended": {"error_model": "qcal", "multiple_testing": "control"},
```

With two or more controls, `_control_fdr` (`variants.py:123`) runs a
target–decoy procedure:

- **Decoys.** Each control is called against the *other* controls (never itself),
  giving a sorted list of its candidate p-values and its window bases. The decoys
  a sample is judged with are its peers' (`project.decoys_for`, `:412`). A
  control far noisier than the rest can be set aside (`project.noisy_decoys`).
- **The rule.** For every observed p-value `t`, the expected false calls are the
  mean over decoys of their candidates at p ≤ t, scaled **up** (never down) to
  this sample's bases. Keep the **largest** `t` at which
  `expected_false ≤ fdr_alpha × n_called(t)`.
- `fdr_alpha` defaults to **0.05**.

So the threshold is per-library and adapts to how much real signal is present:
a library with many true calls has a large `n_called(t)`, which permits a
proportionally larger expected false count, which permits a looser `t`. **This is
empirical control-based threshold selection targeting a nominal false-discovery
proportion, not a proven FDR guarantee (a counterexample is in the October 2026
audit), and it is the single most consequential design decision in the software.** Its consequence is discussed in Part 3.

Without decoys, `control` falls back to Benjamini–Hochberg over
`analysed_positions × 3` hypotheses.

### 1.7 After the test

- optional strand filter (`off` / `both` / `bias`);
- optional homopolymer masking;
- germline differences supplied by the project (for Sanger samples with too few
  clones to find their own).

### 1.8 Single-read mutations (`singletons.py`)

A mutation carried by one read cannot be told from a miscall by any per-site
test. Rather than call them one by one, the default reports an **aggregate
estimate**:

```
estimate = (observed − expected) / efficiency
```

where *observed* is substitutions carried by exactly one read at quality ≥
`single_read_min_quality` (default **35**) and not already called; *expected* is
the same count in controls, per base, scaled to this sample; *efficiency* is the
share of this sample's correctly read bases at quality ≥ the cut, correcting for
real mutations the cut removes.

The justification for the quality cut is that **errors cluster at low quality and
real bases do not**. The check that the efficiency correction is right is that
the estimate is stable across cuts (Q30 / Q35 / Q40 give 50.3 / 50.2 / 51.3 for
Wt2). See `docs/calibration/SINGLE_READ.md`.

**Note the asymmetry this creates**, because it is the subject of Part 4: the
aggregate estimate applies a Q35 floor to one-read events; the *caller* applied
no quality floor at all. A one-read candidate was admitted on its p-value alone.

### 1.9 Clones (`clones.py`)

Reads are grouped by the V(D)J junction read from the unaligned 5′ sequence,
within `junction_max_edits`. A clone's shared mutations are taken as the clone's
and the clone counts once. This gives counts per clone, not per molecule.

### 1.10 Named modes

```python
"recommended": qcal + control-FDR            (alpha 0.05)
"sensitive":   qcal + control-FDR, alpha 0.20
"geneious":    Geneious p-value, uncorrected at max_pvalue, no phasing
"sanger":      every substitution at Q >= 20, no test
"v0.1", "qcal", "strand", "pbinom": variations with BH
```

---

## Part 2. The claim under review

The Results paragraph for Figure 3D–E, as it currently reads:

> We next simulated approximately 40 B-cell clones per J<sub>H</sub>4 library,
> with clone sizes following a power-law distribution. Mutations were introduced
> at positions sampled according to S5F mutability, including mutations shared
> across an entire clone, restricted to individual branches, or unique to
> individual sequences (Figure 3D). SHMtrace recalled 83% of single-read
> mutations versus 57% for Geneious, with FDR 0.037 versus 0.012 and 13.7 versus
> 4.0 false calls per library. In null libraries, false calls averaged 0.7 versus
> 4.3. …

The question put to me was whether the 0.037-vs-0.012 and 13.7-vs-4.0 gaps are a
fixable defect in SHMtrace.

Elsewhere in the manuscript the same tension appears in the J<sub>H</sub>4
spike-in: SHMtrace recalls **30%** of mutations planted in a single read against
Geneious's **60%**, at 2.4 versus 4.9 false calls per library.

---

## Part 3. What was measured

### 3.1 Method

Three J<sub>H</sub>4 *Aicda*<sup>−/−</sup> libraries, each held out in turn; its
reads become the simulated sample and the other two are its controls. Mutations
are planted by `lineage_sim.design` / `build_library`. Every run below uses
identical planted truth for every caller, so the comparisons are paired.

Instrumentation: `variants._control_fdr` was wrapped to record the p-value list,
the decoy lists and the chosen threshold; `_candidate_keys` already records every
candidate the test saw as `(position, alt, p)`. Joining those to the planted
truth gives, for any threshold, the exact true and false calls and the recall
split by how many reads carry the mutation. No re-run is needed to sweep a
threshold.

Scripts are in the scratch directory, not the repository; the reproduction
commands in §6 use only repository code.

**These re-runs use seeds 1..3 and give 12.33 false calls per library for the
published settings, against 13.67 in the archived
`out/lineage_sim_jh4/lineage_sim.json`.** The difference is draw-to-draw. All
comparisons below are internally consistent (same harness, same seeds); the
archived number is the one the manuscript quotes.

### 3.2 Finding 1 — observed false-discovery proportions were below the target in these simulations

The estimated-FDR curve was traced against the true false calls at every
threshold, per library:

| library | chosen t | calls | true FP | estimated FDR at t | **true FDR at t** |
|---|---:|---:|---:|---:|---:|
| KO1 | 0.076 | 351 | 11 | 0.046 | 0.031 |
| KO2 | 0.146 | 380 | 8 | 0.050 | 0.021 |
| KO3 | 0.332 | 375 | 18 | 0.048 | 0.048 |

The curve is monotone in all three, the rule stops where its own estimate says to
stop, and the estimate tracks the truth. **The realized false-discovery proportion
of 0.037 was below the nominal 0.05 in these evaluated simulations.** That is
evidence for this setting, not general FDR control.

Geneious has no target. Its fixed p-threshold happens to land at ~0.012 on this
benchmark with lower recall.

### 3.3 Finding 2 — the chosen threshold varies fourfold between identical libraries

Note the `chosen t` column: 0.076, 0.146, 0.332. These three libraries differ
only in which control's reads they were built from. And the library given the
*loosest* threshold (KO3) is the *cleanest* at a common threshold:

| | KO1 | KO2 | KO3 |
|---|---:|---:|---:|
| false calls at a common p ≤ 0.05 | 8 | 4 | **1** |
| false calls at its own chosen t | 11 | 8 | **18** |

This is not inconsistent with FDR control — the three libraries have different
candidate distributions — but it means the effective sensitivity of the caller is
not a property of the caller.

The variability is driven by `qscale`. Forcing it to lower values (the fitted
value on these controls is 10):

| qscale | per-library thresholds | false calls/lib | one-read recall |
|---:|---|---:|---:|
| 1 | 0.0104, 0.0113, 0.0091 | 11.67 | 0.804 |
| 2 | 0.0267, 0.0236, 0.0300 | 13.00 | 0.828 |
| 3 | 0.0461, 0.0423, 0.0417 | 12.67 | 0.846 |
| 5 | 0.0757, 0.0745, 0.0825 | 12.33 | 0.846 |
| **10 (fitted)** | **0.0757, 0.1462, 0.3323** | **12.33** | **0.852** |

**Performance is flat across the whole range** — the FDR rule compensates for the
inflation — but the threshold stability is much better below 10. This is a
robustness observation, not a performance gain, and nothing was changed on
account of it.

### 3.4 Finding 3 — the two callers are close to the same ROC curve here

Sweeping both callers' own statistics over the same planted truth:

| false calls/lib | SHMtrace one-read recall | Geneious one-read recall |
|---:|---:|---:|
| ~1.3–1.7 | 0.118 | 0.184 |
| ~2.3 | 0.263 | 0.338 |
| ~3.0 | 0.444 | 0.408 |
| ~4.0–4.3 | 0.538 | 0.695 |
| ~5.0–5.3 | 0.710 | 0.755 |
| ~8.0 | 0.761 | ~0.84 |

At matched false calls, Geneious's statistic is as good or better at most points
on **this** benchmark. The published 83%-vs-57% is produced by the operating
point, not by a better-ranked statistic.

**The counterweight is the null scenario**: with no mutations planted, SHMtrace
makes **0.7** false calls per library and Geneious **4.3**, from the same reads.
That is the direct consequence of §1.6 — FDR control tightens when there is
nothing to find — and it is what produces the 38-fold to 229-fold wild-type /
*Aicda*<sup>−/−</sup> separation the manuscript reports. A fixed per-test
threshold cannot do this.

### 3.5 Finding 4 — base quality separates true from false one-read candidates almost completely

For candidates supported by exactly one read, pooled over the three polyclonal
libraries:

| base quality | candidates | true | false |
|---|---:|---:|---:|
| 0–15 | 322 | 5 | **317** |
| 15–25 | 67 | 9 | 58 |
| 25–35 | 37 | 23 | 14 |
| 35–40 | 41 | **39** | 2 |
| 40–50 | 153 | **151** | 2 |
| ≥ 50 | 39 | **37** | 2 |

Mean quality of true one-read candidates: 40–42. Of false ones: 10–11.

Quality enters the candidate statistic; the optional floor adds a separate eligibility rule on top of it. A true Q45 one-read candidate gets
a median p-value of 0.027, so it only passes at a threshold loose enough to also
admit the Q<15 tail.

For comparison, a feature that does **not** work: the carrying read's own
mutation load (how many confirmed multi-read mutations that read carries) is
6.33 for true one-read candidates and 6.53 for false ones. Every read belongs to
a clone with trunk mutations, so load is uninformative. Removing the one-sided
`max(g, site)` floor in `qrate` (§1.5) was also tested and is roughly a wash.

### 3.6 The caveat that limits Finding 4

`lineage_sim.build_library` and `simulate` plant a mutation by substituting the
base and **leaving the quality string alone**. Planted mutations therefore inherit
the original read qualities, without modelling the altered sequence's effect on
basecalling; the sequencing errors they compete against tend to sit on bases the
basecaller got wrong, which is why those are low-quality.

The direction is physically right — that is what basecaller quality means — but
the magnitude is optimistic. A real mutation falling in a low-quality stretch
(homopolymer, poor signal) would be lost by a quality floor, and **no benchmark
here can see that cost**. `docs/calibration/SINGLE_READ.md` already carried this
caveat for the aggregate estimate ("The planted bases are idealised").

This is the single biggest reason to want a second opinion.

---

## Part 4. The change

> **Superseded, 20 October 2026.** Every quality-floor result in Parts 4 and 5 was
> measured while the decoy controls skipped the floor that the samples applied
> (audit finding F01). The decoys therefore kept low-quality singletons that the
> samples had dropped. That inflated the estimated false-discovery proportion
> and tightened the threshold, which made the floor look like it cut false
> calls. With sample/decoy parity restored, the 20 October rerun
> (`out/spikein_<ds>_qfloor.json`, same command as Part 6) gives:
>
> - **At α 0.05**, the floor does not improve precision. On J<sub>H</sub>4 it
>   raises false calls from 2.56 to 4.11 per library (FDR 0.050 → 0.080). It
>   trades two- and three-read recall for one-read recall. On N114 it is
>   neutral.
> - **At α 0.20**, the floor no longer halves the cost of the sensitive mode:
>   J<sub>H</sub>4 13.44 → 13.33, N114 11.00 → 12.56 false calls per library.
>   It raises one-read recall (J<sub>H</sub>4 0.744 → 0.933, N114 0.867 →
>   0.956).
> - The lineage-simulation rows below were not rerun with the floor.
>
> Options (b) and (c) in Part 5 are therefore **not supported**. Keep the floor
> off by default (option a). The tables immediately below are kept only as a
> record of the earlier, asymmetric measurement.
>
> Corrected spike-in tables (9 libraries per dataset, 60 planted sites each):
>
> J<sub>H</sub>4:
>
> | caller | FP/lib | FDR | recall | 1 rd | 2 rd | 3 rd | 5 rd | 10 rd | 20 rd |
> |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
> | recommended (α 0.05) | 2.56 | 0.050 | 0.806 | 0.311 | 0.722 | 0.878 | 0.944 | 0.978 | 1.000 |
> | geneious | 4.89 | 0.083 | 0.896 | 0.600 | 0.867 | 0.933 | 0.989 | 0.989 | 1.000 |
> | recommended + Q30 | 4.11 | 0.080 | 0.787 | 0.433 | 0.678 | 0.733 | 0.922 | 0.978 | 0.978 |
> | α 0.20 | 13.44 | 0.195 | 0.922 | 0.744 | 0.900 | 0.944 | 0.967 | 0.978 | 1.000 |
> | α 0.10 + Q30 | 6.89 | 0.117 | 0.863 | 0.667 | 0.822 | 0.811 | 0.922 | 0.978 | 0.978 |
> | α 0.20 + Q30 | 13.33 | 0.189 | 0.956 | 0.933 | 0.956 | 0.933 | 0.956 | 0.978 | 0.978 |
> | α 0.20 + Q25 | 15.44 | 0.216 | 0.935 | 0.911 | 0.889 | 0.922 | 0.933 | 0.978 | 0.978 |
>
> N114:
>
> | caller | FP/lib | FDR | recall | 1 rd | 2 rd | 3 rd | 5 rd |
> |---|---:|---:|---:|---:|---:|---:|---:|
> | recommended (α 0.05) | 2.56 | 0.045 | 0.900 | 0.578 | 0.911 | 0.933 | 0.978 |
> | geneious | 2.78 | 0.050 | 0.883 | 0.489 | 0.889 | 0.922 | 1.000 |
> | recommended + Q30 | 2.44 | 0.044 | 0.891 | 0.600 | 0.856 | 0.922 | 0.978 |
> | α 0.20 | 11.00 | 0.160 | 0.963 | 0.867 | 0.956 | 0.967 | 0.989 |
> | α 0.10 + Q30 | 5.56 | 0.089 | 0.948 | 0.833 | 0.933 | 0.956 | 0.978 |
> | α 0.20 + Q30 | 12.56 | 0.176 | 0.981 | 0.956 | 0.978 | 0.989 | 0.978 |
> | α 0.20 + Q25 | 13.56 | 0.187 | 0.981 | 0.967 | 0.978 | 0.978 | 0.978 |

`single_read_call_min_quality` (`src/shmtrace/config.py:245`,
`src/shmtrace/variants.py:351`). A candidate supported by exactly one read whose
best supporting base is below the floor is dropped **before** the test, in the
sample and in the decoy controls alike, so the control-based threshold is
re-derived on the survivors. The count removed is reported as
`single_read_low_quality_removed`.

It is **off by default**. Every published number was produced without it. 161
tests pass unchanged. Commit `b37b41d`.

Only `qcal` and `sanger` read base qualities, so it has no effect on the other
models — which means it cannot be applied symmetrically to the Geneious
comparator as implemented.

### 4.1 At the default rate (α 0.05)

| benchmark | setting | false calls/lib | FDR | one-read recall |
|---|---|---:|---:|---:|
| J<sub>H</sub>4 spike-in | published | 2.44 | 0.048 | 0.300 |
| | + Q30 floor | 2.44 | 0.048 | 0.300 |
| N114 spike-in | published | 2.56 | 0.045 | 0.578 |
| | + Q30 floor | **2.33** | **0.041** | 0.578 |
| Lineage sim (re-run) | published | 12.33 | 0.034 | 0.852 |
| | + Q30 floor | **8.67** | **0.024** | 0.825 |

A no-op on J<sub>H</sub>4; free on N114; on the lineage simulation it removes
30% of false calls for three points of one-read recall.

### 4.2 Under a looser rate, where one-read calls are actually won

| benchmark | α 0.20 | α 0.20 + Q30 floor |
|---|---|---|
| J<sub>H</sub>4 spike-in | 11.89 false/lib, FDR 0.177, one-read 0.722 | **7.00**, 0.114, 0.667 |
| N114 spike-in | 11.00 false/lib, FDR 0.160, one-read 0.867 | **6.00**, 0.094, 0.844 |
| Lineage sim | 50.67 false/lib, FDR 0.121, one-read 0.961 | **15.67**, 0.041, 0.897 |

The floor roughly halves what the sensitive mode costs while keeping most of what
it gains. On the lineage simulation it cuts false calls by a factor of three.

Full J<sub>H</sub>4 spike-in table, recall by how many reads carry the planted
mutation (9 libraries, 3 controls × 3 reps, 60 sites each):

| caller | FP/lib | FDR | recall | 1 rd | 2 rd | 3 rd | 5 rd | 10 rd | 20 rd |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| recommended (published) | 2.44 | 0.048 | 0.809 | 0.300 | 0.756 | 0.878 | 0.944 | 0.978 | 1.000 |
| geneious | 4.89 | 0.083 | 0.896 | 0.600 | 0.867 | 0.933 | 0.989 | 0.989 | 1.000 |
| recommended + Q30 | 2.44 | 0.048 | 0.809 | 0.300 | 0.756 | 0.878 | 0.944 | 0.978 | 1.000 |
| α 0.20 | 11.89 | 0.177 | 0.919 | 0.722 | 0.889 | 0.933 | 0.978 | 0.989 | 1.000 |
| α 0.10 + Q30 | 4.22 | 0.077 | 0.844 | 0.444 | 0.789 | 0.889 | 0.967 | 0.978 | 1.000 |
| α 0.20 + Q30 | 7.00 | 0.114 | 0.906 | 0.667 | 0.889 | 0.933 | 0.967 | 0.978 | 1.000 |
| α 0.20 + Q25 | 9.56 | 0.149 | 0.911 | 0.700 | 0.889 | 0.933 | 0.967 | 0.978 | 1.000 |

N114, same layout:

| caller | FP/lib | FDR | recall | 1 rd | 2 rd | 3 rd | 5 rd |
|---|---:|---:|---:|---:|---:|---:|---:|
| recommended (published) | 2.56 | 0.045 | 0.900 | 0.578 | 0.911 | 0.933 | 0.978 |
| geneious | 2.78 | 0.050 | 0.883 | 0.489 | 0.889 | 0.922 | 1.000 |
| recommended + Q30 | 2.33 | 0.041 | 0.900 | 0.578 | 0.911 | 0.933 | 0.978 |
| α 0.20 | 11.00 | 0.160 | 0.963 | 0.867 | 0.956 | 0.967 | 0.989 |
| α 0.10 + Q30 | 4.22 | 0.071 | 0.920 | 0.667 | 0.944 | 0.933 | 0.978 |
| α 0.20 + Q30 | 6.00 | 0.094 | 0.959 | 0.844 | 0.956 | 0.967 | 0.989 |

### 4.3 On the real J<sub>H</sub>4 libraries

Not a benchmark — there is no truth — but it shows what the floor touches:

| setting | Wt2 | Wt3 | Wt4 | KO1 | KO2 | KO3 |
|---|---:|---:|---:|---:|---:|---:|
| published (α 0.05) | 221 | 344 | 73 | 0 | 2 | 0 |
| + Q30 floor | 217 | 341 | 73 | 0 | 2 | 0 |
| α 0.20 | 249 | 387 | 86 | 0 | 3 | 0 |
| α 0.20 + Q30 floor | 230 | 346 | 82 | 0 | 3 | 0 |

At α 0.05 the floor costs 4, 3 and 0 wild-type calls and no control calls. At
α 0.20 it removes 19, 41 and 4 wild-type calls while the controls are unchanged
at 0, 3, 0 — so on real data the floor is removing calls the controls say nothing
about, and whether those were real mutations cannot be determined.

---

## Part 5. What I want a second opinion on

1. **Is the quality floor legitimate given §3.6?** The benchmarks that support it
   cannot measure its cost, because every planted mutation sits on a
   confidently-read base by construction. Is there a design that could measure
   it — planting mutations with a quality drawn from the empirical distribution
   of *correctly read* bases at that position, say?

2. **Is the asymmetry in §3.4 the right framing for the paper?** SHMtrace loses
   to Geneious on the recall/false-call frontier in a heavily mutated simulated
   library and beats it sixfold in a null library, because one uses FDR control
   and the other a fixed threshold. The manuscript currently quotes only the
   first half. The honest statement seems to be that the comparison is between a
   method with a declared error rate and one without, but I may be rationalising.

3. **Should the per-library threshold be stabilised (§3.3)?** A fourfold spread
   in effective sensitivity across three libraries of the same design is
   uncomfortable for a tool whose output is compared between animals. Lowering
   `qscale` stabilises it at no cost in this benchmark. But `qscale` is fitted
   for a reason (cross-library honesty of the p-values), and overriding it on the
   strength of one simulation seems wrong.

4. **Which operating point should the paper report?** Three options:
   (a) leave the default at α 0.05 with the floor off and correct the prose;
   (b) turn the floor on at α 0.05 — strictly better precision, three points of
   one-read recall;
   (c) turn the floor on and loosen to α 0.10–0.20 — fixes the 30%-vs-60%
   J<sub>H</sub>4 one-read gap at 7.0 false calls per library. Each later option
   means re-running every benchmark in the paper.

5. **Is there a discriminator I have missed?** Read mutation load does not work
   (§3.5). AID targeting context (S5F) would work but is circular here, since the
   simulation plants by S5F. Clone membership is already used for grouping, not
   for calling.

---

## Part 6. Reproduce

```bash
# spike-ins with the new variants (writes out/spikein_<ds>_qfloor.json)
python docs/calibration/crm/spikein.py jh4 n114 --no-lofreq --tag _qfloor \
  --profiles recommended,geneious,recommended_q30,recommended_fdr20,recommended_q30_fdr10,recommended_q30_fdr20,recommended_q25_fdr20
```

```bash
# the published lineage simulation
shmtrace-simulate-lineages --region jh4_intron \
  --control KO1.fastq --control KO2.fastq --control KO3.fastq \
  --reference mouse_jh4_intron.fasta --contaminant mouse_jh4_contaminant.fasta \
  --targeting-model HH_S5F.csv --out lineage_sim
```

The threshold sweeps in Part 3 need the instrumented harness, which is not in the
repository; `lineage_sim.roc` (`src/shmtrace/lineage_sim.py:322`) does the same
sweep for mutations carried by ≥ 2 reads and is already written to
`out/lineage_sim_jh4/lineage_sim.json`.

Relevant reading already in the repository: `docs/calibration/SINGLE_READ.md`
(the aggregate estimate and the new floor), `docs/calibration/ERROR_MODELS.md`
(why `qcal` over the alternatives), `docs/calibration/FINDINGS.md` (the
background measurement), `docs/calibration/GENEIOUS.md` (how the comparator is
reproduced), `docs/calibration/CLONES.md`, `docs/calibration/RECAL.md`.
