# Single-read mutations

*11 October 2026. Data: the six bundled J꜀4 libraries (`examples_compare/`, run DHJGKB, 497 reads each).*

## The problem

At 500 reads, no per-site test can confirm a mutation carried by one read without
also calling errors. In the spike-in of [`ERROR_MODELS.md`](ERROR_MODELS.md) the
recommended mode finds 3% of planted one-read mutations on J꜀4. Geneious's rule
finds 80%, but those calls are where its 3 to 6 false calls per sample come from.
In a polyclonal assay, one-read mutations are much of what "unique mutations"
counts.

## The method

Stop deciding them one by one, and count them instead (`shmtrace/singletons.py`):

| Quantity | Definition |
|---|---|
| observed | (position, alt) pairs carried by exactly one read at base quality ≥ Q35, not already called by the per-site test, outside germline and control-polymorphic sites |
| background | the same count in the controls per base sequenced, times this sample's bases sequenced |
| efficiency | share of this sample's correctly read bases at quality ≥ Q35, which is the share of real one-read mutations the cut keeps |
| estimate | (observed − background) / efficiency, with a 95% interval |
| combined density | (all per-site calls + estimate, floored at 0) per base sequenced × 10⁵ |

Two things matter:

- Nanopore errors sit mostly at Q8–Q11 and real bases do not
  ([`CONFOUNDING.md`](CONFOUNDING.md)). Above Q35 the background is about 4
  events per library, so an excess of a few dozen stands out.
- Low-quality reads that show the same alt are ignored. At about 500 reads a
  column often holds one or two such errors. A first version required the alt
  in exactly one read overall, and it recovered only 85% of planted mutations.

The interval combines Poisson noise in the sample's count, between-library
overdispersion of the background (estimated when there are three or more
controls), and the uncertainty of the control count the background rate comes
from. Samples and controls must be at the same depth. The project already does
this.

## Spike-in

`single_read_spike.py` plants N one-read mutations into one *Aicda*-/- library
at AID-weighted sites, keeping each read's own base quality. It runs the
project with the other two as controls, and compares the estimate with that of
the same library unspiked. Two draws per N, each knockout in turn.

| Planted | Change in estimate (6 runs) | Mean recovery |
|---:|---|---:|
| 0 | 0.0 for all | – |
| 25 | 21.4, 25.0, 22.7, 20.3, 19.0, 27.3 | 91% |
| 50 | 51.2, 50.0, 50.1, 54.9, 48.7, 54.6 | 104% |
| 100 | 102.3, 90.4, 91.9, 97.9, 96.2, 99.8 | 97% |

The unspiked knockouts estimate 1.2, −0.7 and −0.6 one-read mutations against
their peers.

## Real libraries

Q35 cut, each knockout measured against the other two, the wild types against
all three.

| Library | Per-site calls | One-read events observed | Background | Estimate | Combined density ×10⁻⁵ (95% CI) | Recommended density | Geneious density |
|---|---:|---:|---:|---:|---|---:|---:|
| Wt2 | 146 | 49 | 7.4 | 50.2 | 74.9 (68.4–81.4) | 55.7 | 60.9 |
| Wt3 | 334 | 26 | 7.3 | 22.8 | 138.1 (133.2–143.1) | 129.3 | 129.3 |
| Wt4 | 68 | 2 | 7.3 | −6.4 | 26.0 | 26.0 | 25.7 |
| KO1 | 0 | 8 | 7.0 | 1.2 | 0.46 (0–3.5) | 0 | 0.76 |
| KO2 | 3 | 7 | 7.5 | −0.7 | 1.15 (1.15–3.8) | 1.15 | 2.27 |
| KO3 | 1 | 7 | 7.5 | −0.6 | 0.38 (0.38–3.1) | 0.38 | 1.13 |

What it shows:

- Wild-type signal comes back to Geneious's level or above, mostly in Wt2, where
  the recommended mode lost the most against Geneious ([`GENEIOUS.md`](GENEIOUS.md)).
  Knockout background stays below Geneious's.
- Wt4 has almost no one-read mutations. Its mutations sit in multi-read
  clones. Its count of 2 against an expected 7.3 is about a 2% tail, plausible
  once in six libraries. The combined density floors the estimate at zero.
- WT against KO, Welch on log densities: P = 0.0022 for the combined density
  against 0.0017 for the recommended density. Both are n = 3 against 3.

## Does the cut matter?

It should not, once the efficiency correction is applied:

| Cut | Wt2 estimate | Wt3 estimate | Background per library | Knockout estimates |
|---|---:|---:|---:|---|
| Q30 | 50.3 | 20.3 | 13 | −0.5, −5.1, +6.7 |
| **Q35** | **50.2** | **22.8** | **7** | **+1.2, −0.7, −0.6** |
| Q40 | 51.3 | 21.6 | 3 | 0.0, 0.0, 0.0 |

The estimate is stable across cuts, and that stability is the check that the
efficiency correction is right. Q35 is the default for nanopore. Q30 is the
default for Illumina, whose qualities top out near Q40. If a cut keeps under
half of any library's correctly read bases, as with an older basecaller, it is
lowered in steps of 5, but not below Q20. The result records the cut used
(`q_min`) and the one asked for.

## The same floor under one-read *calls*

`single_read_call_min_quality` applies this floor to the one-read mutations the
caller reports, not just to the aggregate estimate. A one-read candidate whose
best supporting base falls below it is dropped before the test, in the sample
and in the decoy controls alike, so the control-based threshold is re-derived on
what is left. It is **off by default**: every published number was produced
without it.

What it buys, measured on three benchmarks (`crm/spikein.py --profiles
recommended_q30,recommended_q30_fdr20`, and the polyclonal lineage simulation):

| Benchmark | Setting | False calls / library | FDR | One-read recall |
|---|---|---:|---:|---:|
| J_H4 spike-in | published | 2.4 | 0.048 | 0.300 |
| | + Q30 floor | 2.4 | 0.048 | 0.300 |
| N114 spike-in | published | 2.6 | 0.045 | 0.578 |
| | + Q30 floor | 2.3 | 0.041 | 0.578 |
| Lineage simulation | published | 13.7 | 0.034 | 0.852 |
| | + Q30 floor | 8.7 | 0.024 | 0.825 |

At the default false-discovery rate the floor is close to free: it is a no-op on
J_H4, removes false calls at no cost on N114, and on the lineage simulation
removes a third of them for three points of one-read recall.

It matters more under a looser rate, which is where one-read calls are won:

| Benchmark | α 0.20 | α 0.20 + Q30 floor |
|---|---|---|
| J_H4 spike-in | 11.9 false / library, one-read 0.722 | 7.0 false, one-read 0.667 |
| N114 spike-in | 11.0 false / library, one-read 0.867 | 6.0 false, one-read 0.844 |
| Lineage simulation | 50.7 false / library, one-read 0.961 | 15.7 false, one-read 0.897 |

The floor roughly halves what the sensitive mode costs while keeping most of
what it gains. On the real J_H4 libraries at α 0.20 it changes wild-type calls
from 249, 387 and 86 to 230, 346 and 82, with the *Aicda*^-/- controls at 0, 3
and 0 either way.

The idealised-quality caveat above applies here with more force: a planted
mutation always sits on a base the basecaller read confidently, so these
benchmarks cannot show what the floor costs on a real mutation that falls in a
low-quality stretch. That cost is real and unmeasured.

## Limits

- **Controls are required.** With no control, or a single control (which has
  no peer to be measured against), there is no estimate, and the project
  reports the per-site calls alone.
- **One run, one locus.** The background is measured per run, which is why the
  controls must share the samples' flow cell. The method has not yet been
  checked on N114 (the deposited reads do not match the Geneious BAMs) or on
  public data.
- **The planted bases are idealised.** They keep the read's own quality, but
  not any alignment artefact a real mutation might bring.
- **It is a count, not a list.** There is no table of which one-read mutations
  are real, so spectra and hotspot targeting still come from the per-site
  calls. `single_read` in the result keeps each library's observed and
  expected counts by substitution class for a background-subtracted spectrum,
  but nothing plots it yet.

## Reproduce

```bash
python docs/calibration/single_read_spike.py 0,25,50,100 2
```
