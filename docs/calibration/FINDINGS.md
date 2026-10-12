# Calibrating the variant caller against AID-null controls

Run 31 July 2026 on the J꜀4 data. Scripts in this directory reproduce every number.

> **Read this as a record of the July 2026 investigation, not as a description of
> the caller.** Sections 1–4 below describe the caller *as it was before* this
> work. Everything they recommend is now implemented and on by default: read
> phasing (`phase_mode`), Benjamini-Hochberg correction (`multiple_testing`) and
> the control-fitted null. The default mode changed again in October 2026, to the
> quality-calibrated model; see [`ERROR_MODELS.md`](ERROR_MODELS.md) for what the
> caller does today.

**Result: held-out *Aicda*−/− background falls from 9.8 to 1.8 × 10⁻⁵ with no loss of wild-type signal.** The manuscript's fully filtered value is ~1.2 × 10⁻⁵.

> **Implemented 31 July 2026.** All three corrections are now in the caller and on by default. Measured on the *unphased original* FASTQs, which is the harder and more realistic case:
>
> | | old caller | corrected | + `--null` |
> |---|---|---|---|
> | KO background (held out) | 30.90 | 1.77 | **1.01** |
> | WT signal | 88.8 | 73.4 | 77.3 |
> | separation | 2.9× | 41.4× | **76.3×** |
>
> Phasing finds the same 59-site block in all six samples. `--no-phasing --no-correction` reproduces the old numbers exactly.

| strategy | KO background (held out) | WT signal | ratio |
|---|---|---|---|
| A baseline (current caller) | 9.76 × 10⁻⁵ | 70.4 | 7.2× |
| B + min 3 supporting reads | 6.59 | 63.1 | 9.6× |
| C + Benjamini–Hochberg FDR 5% | 2.28 | 62.6 | 27.4× |
| **D empirical null + BH** | **1.77** | **71.3** | **40.3×** |
| E D + min 3 reads | 1.77 | 60.8 | 34.3× |
| *manuscript target* | *~1.2* | *~59* | *~49×* |

Leave-one-out across the three *Aicda*−/− mice: the error model is fitted on two and evaluated on the third, so the background figure is genuinely held out. Wild-type samples are never used for fitting.

---

## 1. The dominant problem is not per-base error

Comparing the three *Aicda*−/− mice at full depth, the apparent error is almost perfectly reproducible between animals:

| | raw r | rank r |
|---|---|---|
| KO1 vs KO2 | 0.998 | 0.494 |
| KO1 vs KO3 | 0.989 | 0.527 |
| KO2 vs KO3 | 0.988 | 0.458 |

and it is extremely concentrated: **the top 1% of candidate slots carry 37% of all error, the top 5% carry 92%.** Independent random error would not behave like this.

Inspecting the worst slots showed something more specific. Sixteen positions across the window all sat at ~10.1–10.4% alt fraction in every KO mouse, spanning different reference bases and different substitution types. A uniform fraction across unrelated positions is not a base-calling phenomenon.

Read-level check (`haplotype_check.py`), KO1:

| reads carrying… | count | share |
|---|---|---|
| none of the 16 sites | 1467 | 83.2% |
| 1–3 of them | 17 | 1.0% |
| ≥ 12 of them | 280 | 15.9% |
| **all 16** | **249** | **14.1%** |

All-or-nothing. This is a **distinct read population** — a second haplotype, a paralogous locus, or residual contaminant that survives competitive alignment — not per-base error. Roughly 14–16% of KO reads and 17–20% of WT reads.

Running the paper's own phased FASTQs through the current caller confirms the size of the effect:

| | unique variants | density |
|---|---|---|
| KO2 unphased | 75 | 61.3 × 10⁻⁵ |
| KO2 phased | 24 | 9.2 |
| KO3 unphased | 78 | 56.8 |
| KO3 phased | 22 | 8.4 |

**Phasing alone accounts for ~7× of the background.** nanoshm had no phasing step at the time; the manuscript pipeline did. It has one now. Everything below is measured on already-phased reads, so it is additional to this.

## 2. No multiple-testing correction *(fixed: BH is now the default)*

The caller then tested 529 positions × 3 alternative bases = **1,587 hypotheses per sample** at p ≤ 0.05 with no correction. That alone expects ~79 false positives by chance, which is the right order of magnitude for the observed background.

Benjamini–Hochberg at FDR 5%: background 9.76 → 2.28 × 10⁻⁵, a 4.3× reduction, with wild-type signal essentially unchanged (70.4 → 62.6). This is the single cheapest fix and it needs no control sample.

## 3. The quality-derived null is optimistic *(fixed: a control-fitted null is now the default)*

`λ = Σ 10^(−Q/10)/3` assumes errors are independent across reads. Measured against the *Aicda*−/− data, observed error runs **3.3× higher** than that prediction (0.417% vs 0.125% mean per-alt rate). Nanopore error is context-driven and correlated: the same position errs in many reads at once, so an independence assumption understates λ exactly where it matters.

Replacing it with a per-(position, alt) rate estimated from AID-null controls, shrunk toward a weak prior, and testing candidates with a binomial against that rate: background 2.28 → 1.77 × 10⁻⁵, and wild-type signal is *not* reduced (62.6 → 71.3, since at some positions the empirical rate is lower than the quality-derived one).

## 4. A minimum-read floor is the wrong lever *(kept: `min_variant_reads` is still 1)*

Raising `min_variant_reads` to 3 helps the uncalibrated caller (9.76 → 6.59) but once BH and the empirical null are in place it buys nothing on background (1.77 → 1.77) and costs 15% of the signal (71.3 → 60.8). It is a blunt proxy for the statistics being wrong. Leave it at 1.

---

## Recommended implementation order *(all done; dates in [`../../CHANGELOG.md`](../../CHANGELOG.md))*

1. **Read-phasing / haplotype partition before pileup.** Largest single effect (~7×). Cluster reads by shared variant profile; treat a coherent co-occurring block as a separate population rather than averaging it into one pileup. Needs designing — see Phase 1 in `../ROADMAP.md`.
2. **Benjamini–Hochberg across all candidate slots.** ~4× further, a few lines, no control sample required. Should be the default.
3. **`--null` option taking one or more AID-deficient FASTQs**, fitting the per-position empirical rate and using it as the test null. Falls back to the quality model when absent.
4. Leave `min_variant_reads` at 1.

Steps 2 and 3 together take the held-out background to within ~1.5× of the manuscript value while retaining full wild-type signal. Step 1 is what closes the rest, and it is also what makes the tool safe on unphased input from other labs.

## Caveats

Later work supersedes parts of this page: [`ERROR_MODELS.md`](ERROR_MODELS.md) (the quality-calibrated default), [`SINGLE_READ.md`](SINGLE_READ.md), [`INDELS.md`](INDELS.md), [`SANGER.md`](SANGER.md), [`ILLUMINA.md`](ILLUMINA.md) and [`ECTOPIC_AID.md`](ECTOPIC_AID.md).


- One locus, one chemistry, six mice. The empirical rates are J꜀4-specific; the *approach* transfers, the numbers do not.
- The residual gap to 1.2 × 10⁻⁵ is not yet explained. Candidates: the manuscript's phasing was more aggressive than what is baked into the bundled examples, and Geneious applies filters that have not been reproduced exactly.
- Wild-type density varies substantially between animals here (34, 53, 124 × 10⁻⁵), so "WT signal" in the table is a mean over three mice and is noisy.
- The 16-site haplotype block was characterised on KO1 and WT3 only.

## Reproducing

```bash
python docs/calibration/build_pileups.py     # full-depth pileups, all six samples
python docs/calibration/haplotype_check.py   # the read-level all-or-nothing result
python docs/calibration/calibrate_loo.py     # the leave-one-out table above
```

Paths at the top of each script point at the deposit FASTQs and will need editing if the data moves.
