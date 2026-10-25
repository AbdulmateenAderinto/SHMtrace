# What was done about the problems ACCURACY_BREADTH found

*25 October 2026. Simulated measurements: `accuracy_breadth.py run fixes`, with every row in `accuracy_breadth/report.md` under "The opt-in caller fixes, against their own base". Real-data measurements: `public_illumina.py hmces_het ldha_margin`, against the published `hmces` and `ldha` runs.*

[`ACCURACY_BREADTH.md`](ACCURACY_BREADTH.md) ran the truth-known design on every public set with three or more mutation-free controls. The nanopore numbers carried over to four further loci; three problems turned up on the Illumina sets. Chasing them turned up a fourth, larger one.

Everything here is **off by default**, so no published number moves.

| | Parameter | Default | Worth using? |
|---|---|---|---|
| A control's second allele is tested as if it were somatic | `mask_control_het` | `False` | **yes** — inert on 8 of 10 sets, and it corrects the published Ramos frequencies |
| A mutation near a read's end is clipped away, so the frequency there reads low | `read_end_margin_bp` | `0` | **no** — the bias is real but removing it buys nothing and costs false calls |
| Base qualities that say nothing about accuracy | *(always reported)* | — | diagnostic, costs nothing |
| A calibration that could not be satisfied was reported as if it had been | `qscale_max` | `50.0` | the warning costs nothing; **raising the cap is a trade, not an improvement** |

## A. A control's second allele: `mask_control_het`

**The problem.** The frequency rule removes a difference from the reference carried by at least `germline_freq_threshold` (50%) of a library's reads. That catches a homozygous difference. A heterozygous one sits *near* 50% and can fall either side. The Ramos IGH amplicon has 29 such sites between positions 132 and 428, and all three AID-KO clones carry the alternate allele at nearly the same fraction — measured on the published run's own reads, at its depth:

| AID-KO clone | Median allele fraction over the 29 sites | `phase_max_drop` (0.45) |
|---|---:|---|
| AIDKO_1 | 0.449 | below it, so phasing drops the block |
| AIDKO_2 | 0.445 | below it, so phasing drops the block |
| AIDKO_3 | **0.470** | **above it, so phasing refuses** |

Every fraction is under the 50% cut, so the frequency rule never sees them, and in the paper's own Ramos run 28 of clone 3's 36 "substitutions" are these sites, with seven of the twelve Wt and HMCES-KO clones gaining 3-7 each.

**The fix.** A (position, allele) that any control library reads in `control_het_min_frequency` (0.25) to `germline_freq_threshold` of its reads, over at least `control_het_min_depth` (20) reads, is a second allele there and not somatic. `null_model.fit_aligned` collects them while fitting the background; `variants._candidate_exclusion` drops them from the samples' candidates and from the controls' decoys alike, so both sides lose exactly the same sites. Reported as `control_het_removed` and `control_het_alleles`.

One control is enough to establish a site, rather than a majority, so a variant private to one library is caught too. Sequencing error does not reach a quarter of the reads at a covered site, so the band's lower edge is well clear of it; what the rule can lose is a real mutation that a contaminated control also carries at that level. Anything at 50% or more is still germline, counted as before, and the rule is off under a declared `ancestral_reference`, where a change driven high is the experiment's subject.

**What it does to the published Ramos numbers.** `hmces_het` is the published run with the mask on and nothing else changed:

| Clone | Group | Unique mutations | Mutation frequency ×10³ | G:C % |
|---|---|---|---|---|
| HMCES_1 | HMCES_KO | 81 → 74 | 8.01 → **0.51** | 93.8 → 97.3 |
| HMCES_2 | HMCES_KO | 52 → 38 | 7.60 → **1.26** | 88.5 → 92.1 |
| HMCES_3 | HMCES_KO | 68 → 61 | 7.77 → **0.34** | 92.6 → 96.7 |
| HMCES_4 | HMCES_KO | 55 → 44 | 6.73 → **0.21** | 87.3 → 93.2 |
| HMCES_6 | HMCES_KO | 65 → 58 | 8.14 → **0.25** | 86.2 → 89.7 |
| WT_1 | WT | 86 → 79 | 8.17 → **0.51** | 88.4 → 89.9 |
| WT_4 | WT | 83 → 81 | 3.82 → **1.66** | 94.0 → 95.1 |
| HMCES_5, WT_2, WT_3, WT_5, WT_6 | — | unchanged | 0.45-0.82, unchanged | unchanged |
| **AIDKO_3** | AID_KO | **36 → 36** | **31.19 → 31.19** | 61.1 |

The published frequencies are bimodal: five clones at 0.4-0.8 per 10³ and seven at 3.8-8.2. Masking leaves all twelve at 0.21-1.66. A handful of allelic calls per clone is enough to do that, because each one is carried by about 46% of 20,000 reads.

**Why AIDKO_3 is untouched: `phase_max_drop`, not the allele ratio.** The threshold is 0.45 and the clones' fractions straddle it, so phasing removes the second allele from clones 1 and 2 and refuses for clone 3. Clone 3's peers are clones 1 and 2, both phased clean, so *its* null carries no second allele to mask and the mask does nothing for it. The twelve samples are called against all three controls, clone 3 among them, so their sites are masked — which is the whole of the table above.

Raising `phase_max_drop` to 0.55 takes clone 3 from 31 calls at 31.20 per 10³ to **1 call at 0.009**, confirming the second allele is all of it. It is not a free change: it also discards about 45% of the reads from every library, and takes WT_1 from 14 calls to 176 at the same depth.

So AIDKO_3 is not a noisy control. It is the one clone whose second allele survived a threshold the others fell below. Nothing on this page fixes that on its own.

**What the mask does not fix, and a correction.** `ACCURACY_BREADTH.md` reported Ramos recall at 0.13 in simulation, with ten of twelve libraries calling nothing, and attributed it to these sites filling the controls' decoys. **That attribution was wrong.** With the mask on, the same libraries score identically, library for library:

| | base | `mask_control_het` |
|---|---:|---:|
| Candidates tested | 915 | 887 |
| Sites dropped as a second allele | 0 | **28** |
| Control-based threshold | none | none |
| Quality-scale inflation | 50.0, unsatisfied | 50.0, unsatisfied |

The mask removes the 28 sites exactly as intended, but the calibration is the blocker — see the last section. The 28 sites are 3% of the candidates. The mask corrects the published numbers, not the simulated recall.

**What it removes, set by set.** Run over every set's own controls, the mask is inert on eight of the ten, the paper's two among them:

| Set | Controls | Window | Sites masked | (germline, ≥50%) |
|---|---:|---:|---:|---:|
| Illumina Ramos VDJ | 3 | 421 bp | **29** | 0 |
| ONT Sμ, CH12 | 3 | 651 bp | **2** | 9 |
| Illumina J꜀4; Illumina CH12 VDJ | 4 | 70, 229 bp | 0 | 0 |
| ONT Sμ LPS blasts; *Pim1*; OT1; OT2 | 3 | 651-1,008 bp | 0 | 0 |
| ONT J꜀4, N114 (this study) | 3 | 529, 615 bp | 0 | 15, 0 |

A mutation-free control set has nothing at a quarter to a half of its reads, which is what those zeros say. The two CH12 Sμ sites are the private variants in control KO3 that `ACCURACY_BREADTH.md` traced a 1.5× overestimate to.

## B. Mutations at a read's end: `read_end_margin_bp` — not recommended

**The problem.** minimap2 soft-clips a mismatch near a read's end rather than align it: clipping costs less than the mismatch penalty when few aligned bases lie beyond it. A read whose mutation is in its first few bases loses both the mutation and its bases at those columns, while its mutation-free neighbours are still counted there. Those columns are then covered by fewer mutant reads than mutation-free ones and read a lower frequency than the library carries.

It matters only where reads begin *inside* the analysis window. On Illumina J꜀4 the reads begin between positions 55 and 66 in a 53-122 window, and 4.5% of planted copies are lost this way at these studies' loads (19% at the paper's J꜀4 loads). Amplicons whose reads span the reference, and the sets with `end_trim_bp` 20, do not have it: in simulation ONT J꜀4 is identical to three decimals at 4 bp and 8 bp.

**The fix.** `read_end_margin_bp` ignores that many aligned columns at each end of every read, in samples and controls alike. It cannot recover the clipped mutation — that is not in the alignment to recover — but it removes the clipped reads' deficit from the denominator as well as the numerator. `align.pair_window` is the one place that decides which columns a pileup may count, and every pileup feeding a call goes through it: the sample's (`variants`), the controls' background and decoys (`null_model`), the per-read table (`pipeline`), the clone tallies (`clones`) and the single-read estimate (`singletons`). In a controlled construction it works exactly as intended: 88 of 96 planted mutations survive an 8-bp clip, and the frequency goes from 0.917 of truth to within 0.1% (`tests/test_caller_fixes.py`).

**It breaks the fixed-window frequency.** `mutation_frequency_x1e3` credits every read with the whole callable window, covered or not. A margin takes mutations out of its numerator and nothing out of its denominator, so it reads low by roughly the share of each read trimmed — which is the whole of the first simulation pass showing J꜀4 getting *worse* (0.78 of truth at base, 0.71 at 4 bp, 0.65 at 8 bp). Use `mutation_frequency_observed_x1e3`, which divides by the bases each read actually showed, or the densities; `call_variants` emits a note saying so, and the simulator scores the matched pair (`true_observed_x1e3` against `est_observed_x1e3`).

**On the real data it does not help.** `ldha_margin` is the published J꜀4 run with an 8-bp margin: calls move both ways (KO2 8 → 5, WT7 9 → 13), `analysed_bp` stays 70 because the window is fixed, and the controls gain calls (AIDKO_1 2 → 4; AIDKO_2's combined density 0.34 → 0.85). On the matched denominator the simulation agrees: 0.81 of truth at base against 0.80 at 8 bp, with false calls in empty libraries going 0.0 → 0.5. The clipping is real but it is not what makes this set read 20% low. **Leave this at 0** unless a particular library shows the clipping and you have measured the effect on it.

A rescue was considered instead — re-align the clipped bases and keep the mutation — and rejected: for J꜀4 the unaligned 5' sequence is the V(D)J junction the clone caller reads, so extending the alignment into it would invent mismatches and move the junction string.

## C. Qualities that carry no information

**The problem.** The single-read estimate separates real mutations from errors by base quality: a real mutated base is read at ordinary quality, errors cluster low, and the efficiency correction is the share of correctly read bases above the cut ([`SINGLE_READ.md`](SINGLE_READ.md)). That needs the quality to vary with how right the base is. In the two public Illumina mouse sets it does not. On the reads the pipeline keeps:

| Set | Quality levels (≥0.5% of reference-matching bases each) | Commonest |
|---|---|---|
| Illumina J꜀4, PRJNA973112 | **1** (Q60) | 100% at Q60 |
| Illumina CH12 VDJ, PRJNA918596 | **2** (Q30, Q60) | 86% at Q30 |
| Illumina Ramos VDJ, PRJNA820976 | 22 | 42% at Q38 |
| ONT Sμ, LPS blasts | 32 | 7% at Q38 |
| ONT J꜀4 (this study) | 27 | 14% at Q50 |

Q60 is what merging a read pair writes where both mates agree. The level says whether both mates covered the base, not how confident the basecall is, so the quality cut sorts merged reads from unmerged ones instead of correct bases from errors. On CH12 the combined density overshoots the truth 1.3×, and 5× at the lowest load.

**The fix.** A flag, not a correction: with no quality information there is nothing to correct with. `singletons.quality_resolution` counts the distinct qualities carrying at least 0.5% of the reference-matching bases, pooled over the sample and its controls. Fewer than four, or one quality on 99% or more, and the estimate comes back with `quality_flat` and a `warning`; `quality_levels` is reported either way. Both sit in the per-sample `single_read` detail of `project.json`, the warning is copied onto the sample's `warnings` (the column `samples.csv` and the HTML report show), and the row carries `single_read_quality_flat`. A NovaSeq's four bins do vary with accuracy and are left alone.

## D. The larger problem: a calibration that cannot be satisfied

`calibrate_qscale` inflates the quality-calibrated error rates by the smallest factor at which a held-out control stops looking like a sample. Until October 2026 the factors it tried were `(1, 2, 3, 5, 10, 20, 50)` — and 50 was both the last step and the ceiling, so a project needing more was given 50 and reported no differently from one that fitted.

| Set | Factor | Satisfied? | Consequence |
|---|---:|---|---|
| ONT J꜀4, N114 (this study) | 5, 10, 10 | yes | nothing to see |
| Illumina CH12 VDJ | 1.0 | yes | no inflation needed |
| Illumina Ramos VDJ | 20.0 | yes | calls are made (52-99 per clone) |
| Illumina J꜀4 | **50.0** | yes, but only at the last factor tried | the strict calling `ACCURACY_BREADTH.md` left unexplained — 1 of 39 mutations carried by 2-10 reads |
| Ramos, simulated at 5,000 reads | **50.0** | **no** | five of six libraries call nothing |

**What was done.** `QSCALE_GRID` now runs to 1000, and `Params.qscale_max` caps which steps the search may use. The cap defaults to **50.0**, the old grid's whole extent, so every published run calibrates exactly as before. What is new is that exhausting the cap is no longer silent: `qscale_warning` says the quality model does not fit the data, and the project carries it as a note — as does any satisfied factor above `QSCALE_SUSPECT` (50), since by then the reported qualities plainly do not describe the data's errors. `--qscale-max` exposes the cap.

**Raising the cap is a trade, not an improvement.** A larger factor is a stricter test. The J꜀4 ONT case in [`QUANT_ACCURACY.md`](QUANT_ACCURACY.md) — 497 control reads analysed at 100 — is itself at the cap, and lifting it takes that run from 33 true calls to 30, recall 1.00 to 0.87, and the frequency ratio 0.93 to 0.88. That is why the search is capped rather than simply extended.

**A related measurement.** The merged-pair Q60 quality in these files is a fiction. Against the controls' own reference-matching bases, second-allele sites excluded:

| Set | Q60 bases | Observed mismatch rate at Q60 | What Q60 claims |
|---|---:|---:|---:|
| Illumina J꜀4 | 2,932,211 (all of them) | 0.034% | 0.0001% — **338× optimistic** |
| Illumina Ramos VDJ | 2,064,792 (15%) | 0.051% | 0.0001% — **505× optimistic** |

That directly explains section C, whose efficiency correction reads the quality at face value.

**What is not established:** that the Q60 fiction is what forces the inflation. The quality-calibrated model learns its rates from the controls empirically, not from the nominal `10^(-Q/10)`, so an optimistic nominal quality does not by itself inflate it. Nor is it systematic site-level disagreement between clones: outside the 29 second-allele sites, no (position, allele) differs between two Ramos AID-KO clones at p < 10⁻⁶ by Fisher's exact test. A per-site rate fitted on two libraries and applied to a third is noisy, and the target — 99% of 915 candidates above p = 0.01 — is demanding; but that is a hypothesis, not a measurement.

**Still worth trying, not tried here:** clamping merged qualities to a realistic ceiling before the pileup. `indels.calibrate` has the same silent-ceiling shape and has not been touched.

## Limits

- **The sample's own heterozygous sites are untouched.** `mask_control_het` works from the controls. Where a sample carries a second allele its controls do not — the cross-context arm in `ACCURACY_BREADTH.md`, where CH12 samples against primary-cell controls read 2.08× the truth — nothing here helps. Matching controls to samples does.
- **`phase_max_drop` is a knife edge and was not changed.** A minor haplotype near 45% of reads is dropped in one library and kept in another, which is what makes AIDKO_3 look noisy. A project-wide phasing decision, or a threshold clear of the spread, would fix it; neither is implemented or measured.
- **One band, one level threshold.** The het band (0.25 to 0.50) and the four-level quality threshold come from these datasets. Both are parameters.
- **The margin does not reach the optional base-quality recalibration.** `recal.fit` (off by default) still fits its quality map over every aligned column; the two options have not been measured together.
- **Substitutions only.** Indel calling has its own germline rule in `indels.py`, which `mask_control_het` does not reach.
- **No cross-run null.** Everything here was measured on each set's own controls and on simulated libraries; none of it has been tested with a saved panel from another run.

## Reproduce

```
python -m pytest tests/test_caller_fixes.py tests/test_quant_accuracy.py
python docs/calibration/accuracy_breadth.py run fixes --jobs 6
python docs/calibration/accuracy_breadth.py report
python docs/calibration/crm/public_illumina.py hmces_het ldha_margin
```
