# How accurately does SHMtrace quantify SHM? J꜀4 and N114

*12 October 2026. `docs/calibration/quant_accuracy.py`; full output in `quant_accuracy/` (`report.md`, `runs.json`).*

## Design

**Libraries.** Libraries with known mutations were built by `lineage_sim` from each amplicon's own mutation-free controls in the deposit:
- J꜀4: Aicda⁻/⁻ KO1–3, the raw 5,000-read FASTQs.
- N114: empty-vector C3, C4 and D1, with 1,240–1,865 reads each.

The reads therefore carry the sequencer's real errors. Each control was held out in turn, with the other two as the null.

**Depth.** Every library was analysed at 125, 250, 500 and 1,000 reads.

**Scenarios.**
- J꜀4: polyclonal; expanded (one unmutated clone with 60% of reads, the Wt2 situation); low (a light load); null (nothing planted).
- N114: subclonal (a cell line with small subclones, most mutations in one read); low; null.

**Modes.** The recommended (quality-calibrated) model and Geneious's rule. In total, 168 libraries.

For each library, three reported numbers are compared with what was planted:
- **Mutation frequency:** mutations per 10³ bp, read-weighted.
- **Unique density:** distinct per-site calls per base sequenced × 10⁵, the manuscript's measure.
- **Combined density:** unique density plus the estimate of single-read mutations the per-site test can't confirm ([`SINGLE_READ.md`](SINGLE_READ.md)).

## Two bugs found and fixed

Both appear only when the input is deeper than the analysis depth, for example raw FASTQs of thousands of reads analysed at 500. The manuscript's 500-read files were not affected.

**1. Paralogue reads and chimeras got through the contaminant filter.** Full-length paralogue (Agouti-type) J꜀4 reads align better to the 1,230-bp C57BL/6 reference than to the bundled 614-bp paralogue record, so the competitive alignment kept them. They are 1–4% of every raw J꜀4 library, wild type and knockout alike. Each one carries 56–62 paralogue alleles, which were then called as mutations: 35–54 false calls in a simulated 500-read library, against about 2 after the fix.

The documented allele-level check is now implemented (`contaminant_min_alleles`, default 10). A read carrying 10 or more of a contaminant's bases at sites where the contaminant differs from the reference is dropped. Sites the sample itself carries in 20% or more of its reads are germline and are ignored. On the real libraries this separates cleanly:

| | Contaminant alleles per read |
|---|---|
| Real reads | 0–7 |
| Paralogue reads and chimeras | 11–62 |

On the raw libraries at 500 reads, the filter removes 12–68 reads per library and changes the calls as follows:

| Library | Before | After | Why |
|---|---|---|---|
| Wt3 | 371 | 340 | paralogue alleles removed |
| Wt2 | 183 | 205 | the null fitted on the controls no longer includes paralogue alleles, so real calls at those positions pass |
| Knockouts | unchanged | unchanged | |

**2. The single-read estimate counted the whole library.** It counted events in every prepared read, but its density divides by the bases of the reads actually analysed. A 1,500-read library analysed at 500 overstated the combined density by 1.3–2×. It now counts the analysed reads, with the controls drawn at the same depth. `tests/test_quant_accuracy.py` fails without this fix (ratio 1.7–2.2) and passes with it (about 1.0).

## Results after the fixes

Estimate ÷ truth, mean over the three held-out controls. 1.00 is exact.

**Recommended mode**

| Amplicon · scenario | Depth | Mutation frequency | Unique density | Combined density | False calls per library |
|---|---:|---:|---:|---:|---:|
| J꜀4 · polyclonal | 125 / 500 / 1000 | 1.01 / 0.99 / 1.00 | 0.99 / 0.87 / 0.77 | 0.99 / 0.96 / 0.95 | 8.0 / 3.7 / 2.7 |
| J꜀4 · expanded | 125 / 500 / 1000 | 1.01 / 0.98 / 0.98 | 0.98 / 0.83 / 0.70 | 0.99 / 0.95 / 0.93 | 6.7 / 6.3 / 3.0 |
| J꜀4 · low | 125 / 500 / 1000 | 1.00 / 0.99 / 0.96 | 0.97 / 0.74 / 0.42 | 1.01 / 0.94 / 0.97 | 5.3 / 3.0 / 2.3 |
| N114 · subclonal | 125 / 500 / 1000 | 1.08 / 1.02 / 0.99 | 1.11 / 0.93 / 0.80 | 1.11 / 1.00 / 0.97 | 8.0 / 10.0 / 9.3 |
| N114 · low | 125 / 500 / 1000 | 1.05 / 1.01 / 1.00 | 1.11 / 0.93 / 0.78 | 1.11 / 1.00 / 0.94 | 6.3 / 8.3 / 7.7 |

**Empty libraries** (truth 0), the estimate × 10⁵ at 500 reads. Real wild-type J꜀4 is about 28–130 × 10⁵.

| | Unique density | Combined density |
|---|---:|---:|
| J꜀4, recommended | 0.38 | 0.45 |
| J꜀4, Geneious | 1.65 | 1.65 |
| N114, recommended | 0.00 | 0.36 |
| N114, Geneious | 0.54 | 0.69 |

## What this means

- **Mutation frequency** is within 5% of the truth on both amplicons from 500 reads, and within 8% at 125–250. It's the most robust number.
- **Combined density** is within 7% of the truth on both amplicons from 500 reads, and within 11% at 125–250. It is the right measure of distinct mutations when samples differ in depth or load.
- **Unique density**, the manuscript's measure, falls as depth rises: from 0.99 of the truth at 125 reads to 0.42–0.80 at 1,000. The per-site test confirms fewer of the single-read mutations as the number of reads grows. It is a sound comparison only between samples analysed at the same depth, which the manuscript did (500 reads). It should not be compared across depths.
- **Background:** the recommended mode's background in empty libraries is a quarter to a fifth of Geneious's.
- **N114 at low depth:** the recommended mode overshoots by 4–11% at 125–250 reads. That's 6–10 single- or two-read errors per library let through by the 5% false-discovery-rate control, which admits a few errors once a library has real calls; with none, it admits none. Geneious's rule undershoots there instead (0.85–0.94).

## Limits

- Planted bases keep the donor read's quality and have no alignment artefact of their own, so recall here is an upper bound.
- The lineages are simple, and there is no PCR error or amplification bias beyond what the control reads already carry.
- Three backgrounds per amplicon, from one run each.
