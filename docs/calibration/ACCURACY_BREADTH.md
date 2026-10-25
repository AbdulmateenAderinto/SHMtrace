# How accurate is SHMtrace beyond J꜀4 and N114?

*25 October 2026. `docs/calibration/accuracy_breadth.py`. All output is in `accuracy_breadth/`: `report.md` has every row, `runs.jsonl` has every run, and `qmaps.json` has the calibration.*

[`QUANT_ACCURACY.md`](QUANT_ACCURACY.md) measures accuracy against planted truth on two nanopore amplicons. Everywhere else, the paper's evidence is agreement: with published effects, with Geneious and with other tools. Agreement can't reveal an error that both sides make. This page runs the QUANT_ACCURACY design on every public set that has three or more mutation-free controls. It then makes the simulation less ideal, one assumption at a time, to see which idealisations flatter the numbers.

## Design

**Sets.** Eight sets from three platforms and four studies:

| Set | Platform | Locus | Controls |
|---|---|---|---|
| PRJNA973112 | Illumina | J꜀4, 70-bp window | 4 AID-null mice |
| PRJNA918596 | Illumina | CH12F3 VDJ, 229 bp | 4 empty-vector replicates |
| PRJNA820976 | Illumina | Ramos IGH VDJ, 421 bp | 3 AID-KO clones |
| PRJNA1481264 | ONT | Sμ in LPS blasts; Sμ in CH12 | 3 *Aicda*⁻/⁻ each |
| PRJNA1481264 | ONT | *Pim1*; CH12 off-targets OT1 and OT2 | 3 *Aicda*⁻/⁻ each |

The settings are the paper's own (`crm/datasets.py`, `crm/public_illumina.py`).

**Libraries.** Each control is held out in turn, mutations are planted into its reads by `lineage_sim`, and the other controls form the null. The libraries are analysed at two depths (Illumina 1,000 and 5,000 reads; ONT 500 and 2,000), with two draws per held-out control.

**Loads are each study's own.** The main scenario plants the median excess of substitutions in that study's AID-positive libraries over its controls. This is measured from the aligned reads, without calls, so it doesn't depend on SHMtrace. "Low" plants a quarter of that. The off-targets plant all of it, with a floor of 0.2 per 10³ bp.

The real excesses are small: 0.15–2.5 per 10³ bp, against J꜀4's 5.0. So a library holds only a few dozen planted mutations, and ratios to the truth swing widely. `report.md` therefore also gives the estimate minus the truth.

> A first pass planted QUANT_ACCURACY's J꜀4 loads on every locus, which is 5–50× these studies' real loads (17–58 per 10³ bp on the 70-bp Illumina window). Its results described libraries that don't exist and were discarded. J꜀4 and N114 keep QUANT_ACCURACY's scenarios. Those plant 11.6 and 1.8 per 10³ bp, against real excesses of 2.6–22 and 2.9–3.9.

## Breadth: at the depth each study used

Estimate ÷ truth: the mean over 6–8 libraries, with the range in brackets. "Empty" is the null scenario: no mutations planted.

| Set | Scenario | Mutation frequency | Combined density | Recall, ≥2 reads | False calls (empty library) |
|---|---|---:|---:|---:|---:|
| ONT Sμ, LPS blasts, 2,000 | polyclonal | 0.95 (0.91–0.98) | 0.98 (0.75–1.25) | 0.81 | 0.3 |
| ONT Sμ, CH12, 2,000 | subclonal | 0.95 (0.91–1.02) | 0.89 (0.86–0.93) | 0.86 | 0.0 |
| ONT *Pim1*, 2,000 | low | 0.86 (0.73–0.93) | 0.97 (0.79–1.15) | 0.87 | 0.0 |
| ONT OT1, 2,000 | low | 0.85 (0.73–0.93) | 0.89 (0.16–2.14) | 0.70 | 0.0 |
| ONT OT2, 2,000 | low | 0.94 (0.82–1.12) | 1.11 (0.68–2.00) | 0.90 | 0.0 |
| Illumina CH12 VDJ, 5,000 | subclonal | 0.93 (0.47–1.40) | 1.27 (0.99–1.70) | 0.81 | 0.5 |
| Illumina J꜀4, 5,000 | polyclonal | **0.78** (0.47–0.94) | 0.58 (0.18–0.99) | 0.43 | 0.0 |
| Illumina Ramos VDJ, 5,000 | subclonal | **0.15** (0.00–0.88) | 0.73 (0.57–0.87) | **0.13** | 0.0 |

**Nanopore generalises.** On four loci the paper never planted into, two of them AID off-targets, the mutation frequency is within 5–15% of the truth. Combined density is within about 10%, and empty libraries give 0–0.3 false calls at 2,000 reads. Frequency runs slightly low at the lowest loads (0.85–0.88).

**Unique density is not a measure for these sets.** It is 0.02–0.46 of the truth at these depths, because most planted mutations sit in one read. QUANT_ACCURACY found the same, at a milder level, on J꜀4.

**Three Illumina problems:**

1. **Illumina J꜀4 underestimates.** Its reads begin anywhere from position 55 to 66, inside the 53–122 window, and the paper set `end_trim_bp` to 0 for this set. A mutation in a read's first ~8 bases is soft-clipped by the aligner and lost: 4.5% of planted copies at these loads, and 19% at J꜀4 loads. The same happens to real mutations.

   Beyond clipping, at 5,000 reads this set calls nearly nothing carried by 2–10 reads (1 of 39), and 57 of 100 mutations carried by more than 10. I haven't established why its calling is this strict. These libraries are heavily PCR-duplicated ([`ILLUMINA.md`](ILLUMINA.md)), which is one candidate.

2. **Ramos recovers almost nothing (recall 0.13).** Nothing is called in ten of the twelve Ramos libraries; only 2 got any calls (11–12 true calls each), both with AIDKO_1 held out. The control-based threshold is never met, because the leave-one-control-out calibration cannot be satisfied at all here: no inflation up to 50-fold, the largest the search tried, makes one AID-KO clone look like a control against the other two, so the rates are inflated to that ceiling and no candidate survives. **Why the clones are that inconsistent is not established here** — see [`CALLER_FIXES.md`](CALLER_FIXES.md), which also gives the calibration a warning for this case instead of reporting it as a satisfied fit.

   Ramos does also carry a second IGH allele — 29 sites where 44–48% of reads differ from the reference in all three AID-KO clones. `phase_max_drop` is 0.45, so phasing drops the block from the clones below that and refuses for the one above it, and the high-allele-fraction rule removes only sites at 50% or more. These sites are *not* what blocks calling here: masking them (`mask_control_het`, [`CALLER_FIXES.md`](CALLER_FIXES.md)) removes 28 of the 915 candidates and leaves every library's score unchanged. They matter in the paper's own run, below, where the threshold is met and they become calls.

3. **CH12 combined density overshoots (1.3×, and 5× at the low load).** Its base qualities take two values: Q30 where one mate covered the base and Q60 where both did (86% and 12% of the reference-matching bases the pipeline keeps). The level says whether the pair merged there, not how confident the basecall is, so the single-read estimate's quality cut sorts merged reads from unmerged ones rather than correct bases from errors. Illumina J꜀4 is worse: after read filtering it is a single value, Q60 on every base.

## The paper's own Ramos analysis has the same artefact

In the paper's Ramos run (`crm/out/public/hmces/`), AIDKO_3's 36 "substitutions" include **28 second-allele sites called at 45–48%**. Its allele ratio is about 47%, against 45% in the other two controls, which at 20,000 reads is significant at every allelic site. Seven of the twelve Wt and HMCES-KO clones gain 3–7 such calls each. Those calls make up their 7–8-per-10³ mutation frequencies, against 0.4–0.8 in the other five clones.

Removing calls at an allele fraction of 0.3 or more:

| | Published | Without the second allele |
|---|---|---|
| AIDKO_3 | 36 substitutions, 61% at G:C | 8 substitutions, 50% at G:C |
| Calls per Wt or HMCES-KO clone | 52–99 | 45–99 |
| G:C share in those clones | 86–97% | 90–97% |

The manuscript describes AIDKO_3 as a noisy control and the effect as "control assignment materially affects the inferred substitution burden". AIDKO_3's excess is in fact its second IGH allele, kept where the other two clones' was removed: `phase_max_drop` is 0.45, clones 1 and 2 read the allele at 0.449 and 0.445 so phasing drops their block, and clone 3 reads 0.470 so phasing refuses. Masking the allele (`mask_control_het`) takes the other twelve clones' frequencies from a bimodal 0.4-0.8 and 3.8-8.2 per 10³ to 0.21-1.66, and leaves clone 3 untouched, because its two peers were phased clean. [`CALLER_FIXES.md`](CALLER_FIXES.md) has the full tables.

## A null from another cell context

Sμ samples were called against the other cell context's controls (2,000 reads):

| Sample against controls | Frequency ÷ truth | False calls (empty library) |
|---|---|---|
| LPS blasts against CH12 | 0.99 | 1.7 (matched controls: 0.3) |
| CH12 against LPS blasts | **2.08** (1.89–2.29) | 2.0 (matched: 0.0) |

The CH12 line's own heterozygous variants are absent from primary-cell controls, so they are called. SHMtrace's run-mismatch check can't see a mismatch of cell context. The same mechanism appears within one context: CH12 Sμ control KO3 has private heterozygous variants at 45–48%, and with it held out at 500 reads, frequency reads 1.5× the truth.

**A fixed 50% germline cut is fragile for anything heterozygous.**

## Realism: which idealisations flatter the numbers

The same libraries at the main depth, with the planted side changed one step at a time (lineage_sim's docstring has the details):

- **qmap:** each planted base takes the quality that real mutant bases show.
- **pcr3 / pcr10:** 3 or 10 reads per molecule, with log-normal skew, and early-cycle PCR errors at 1×10⁻⁴ per base shared by a molecule's copies. The PCR errors are not part of the truth.
- **deep:** bushy trees, 4 levels.
- **all:** qmap, pcr3 and deep together.

| | qmap | pcr3 | pcr10 | deep |
|---|---|---|---|---|
| J꜀4 ONT, frequency ÷ truth (base 1.00) | 1.00 | 1.02 | 1.02 | 1.02 |
| N114 ONT, frequency (base 0.99) | 0.99 | 1.10 | 1.08 | 1.03 |
| Sμ CH12, frequency (base 0.94) | 0.95 | 1.30 | 1.14 | 1.19 |
| False calls in empty libraries: N114 (base 0.3) | – | 10.0 | 7.3 | – |
| Sμ LPS (base 0.3) | – | 4.3 | 11.0 | – |
| *Pim1* (base 0) | – | 24.3 | 24.7 | – |
| OT2 (base 0) | – | 15.7 | 12.3 | – |
| CH12 VDJ, Illumina (base 0.5) | – | 1.8 | 3.5 | – |

- **Planted-base quality doesn't matter.** Real mutant bases have slightly lower quality than reference bases at the same sites (J꜀4: 9.9% against 9.2% below Q30; Ramos: 10.9% against 4.4%). Mapping planted bases onto that distribution moves no estimate by more than 0.05. QUANT_ACCURACY's first limit can be retired.
- **PCR is what matters.** Errors made in the first PCR cycles appear in several reads, as a real mutation does. When a sample's library carries such errors and its controls don't share them, empty libraries gain 2–25 false calls and low-load frequencies run 7–30% high. J꜀4's load is high enough to absorb this (within 2%).

  The donor reads already carry their own library's PCR errors, so these arms model a sample amplified *more* error-prone than its controls. That's the case to avoid: **controls must be amplified the same way, with the same polymerase and cycles, as the samples.** The magnitude scales with the assumed rate: 1×10⁻⁴ per base is roughly Taq over the first few cycles, and a high-fidelity polymerase would make about 100× fewer such errors.
- **Bushier trees cost little.** Frequency stays within 3% on J꜀4 and N114, but on Sμ in CH12 it reads 1.19× (1.04–1.45).

## What this changes

- For ONT data, the J꜀4 and N114 accuracy holds on four further loci and two cell contexts. This is truth-based evidence where the paper so far has agreement only.
- On the public Illumina sets, the mutation frequencies come with these caveats:
  - Illumina J꜀4 reads 22% low at its main load, and erratically at a quarter of it (read-start clipping, plus calling restricted to well-supported mutations).
  - Ramos substitution calling is unreliable while its second allele is unphased.
  - CH12's combined density is not interpretable on flat qualities.
- All three problems now have an opt-in fix, measured in this same framework: [`CALLER_FIXES.md`](CALLER_FIXES.md) (`mask_control_het`, `read_end_margin_bp`, and a flag on qualities that carry no accuracy information). They are off by default, so the numbers on this page are the unfixed caller.

## Limits

- **Small mutation counts.** Each library holds a few dozen planted mutations, so ratios have wide ranges; use the estimate-minus-truth column in `report.md`.
- **Illumina J꜀4's strict calling: the proximate cause, not the reason.** Its published run satisfies the leave-one-control-out calibration only at 50-fold inflation of the quality-calibrated error rates — the largest factor the search tried — so only well-supported mutations survive. That is the mechanism; *why* this data needs 50-fold where the paper's own nanopore libraries need 10 is not established. The other two public sets settle at 20.0 (Ramos) and 1.0 (CH12). See [`CALLER_FIXES.md`](CALLER_FIXES.md).
- **Context, not run.** Sμ "cross-context" swaps cell context. The run table has no flow-cell IDs, so it isn't known to be a cross-run test.
- **No PCR model in breadth.** The breadth arm plants as QUANT_ACCURACY does: one read per molecule and donor qualities. The realism arm shows which of those choices matter.

## Reproduce

```
python docs/calibration/accuracy_breadth.py cache        # reads to ~/Nanoshm/public_data/accuracy_cache/
python docs/calibration/accuracy_breadth.py qmaps        # loads and quality maps
python docs/calibration/accuracy_breadth.py run breadth --jobs 10     # then crossctx, realism
python docs/calibration/accuracy_breadth.py report
```

566 runs take about 1.5 hours on 10 cores.
