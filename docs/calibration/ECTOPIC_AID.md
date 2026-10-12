# Ectopic AID against empty vector, in an independent study

*12 October 2026. Data: [PRJNA918596](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA918596) ("Mesoscale DNA feature in antibody-coding sequence facilitates somatic hypermutation"). All 48 FASTQ files checked against ENA's MD5.*

This is the public counterpart of the N114 experiment: an AID-null cell line, one clone, transduced with AID or with empty vector, and the same amplicon sequenced from both. It tests the caller where the only difference between groups is the enzyme.

## The design

CH12F3 *Ung*-/-*Msh2*-/-*Aicda*-/- cells, clone 120, VDJ amplicon F2, day 16 after transduction. Four conditions, four replicates each, Illumina paired-end merged to a 269-bp amplicon (97% of pairs merge):

| Condition | What it is |
|---|---|
| **AID** | wild-type AID, ectopically expressed |
| **EV** | empty vector: the AID-null baseline, used here as the control group |
| **R50E** | an AID point mutant |
| **R174ER178E** | an AID double mutant |

The reference is the majority consensus of the empty-vector reads, as for the Sμ reference. Settings: `--platform illumina`, `generic_amplicon`, window 21–249 (the 269-bp amplicon less 20 bp of primer at each end), recommended mode, the four EV replicates as controls, each tested against the other three.

*Ung*-/-*Msh2*-/- matters: with both uracil-processing pathways gone, an AID lesion is left as an unrepaired U:G pair and read as **C→T or G→A only**. The expected signature is therefore narrower than in a wild-type background, and correspondingly easier to check.

## Result, at a matched 5,000 reads

| Group | Mutations per 10³ bp (mean) | Unique calls per replicate | Welch on log values vs EV |
|---|---:|---|---:|
| **AID** | 0.201 | 5, 6, 6, 9 | **P = 0.00017** |
| R174ER178E | 0.628 | 5, 11, 17, 17 | P = 0.00004 |
| R50E | 0.000 | 0, 0, 0, 0 | P = 0.39 |
| **EV** | 0.005 | 0, 0, 0, 1 | – |

- **AID against empty vector: about 40-fold, P = 0.00017.** Across four EV replicates and 20,000 reads the caller made **one** call in total.
- **R50E makes no calls at all**, and is indistinguishable from empty vector. A catalytically inactive AID should be silent, and it is.
- **The AID spectrum is 100% C→T and G→A** (58% C→T, 42% G→A; 26 calls), with 83–100% of G:C mutations at WRC/GYW hotspots. That is exactly the unrepaired-uracil signature this background predicts, and none of it was built into the caller.
- The one EV call is a G→T, the classic oxidative artefact, not a deamination product.

## R174ER178E: higher frequency, wrong spectrum

R174ER178E comes out three times above AID by mutation frequency. It is not AID activity:

- **The spectrum is wrong.** G→C 24%, G→T 20%, A→T 12%, G→A 12%, T→G 10%, A→G 8%. Only 12% of calls are G→A and none is C→T, in a background where AID can produce nothing else. Hotspot targeting is 29–50% against AID's 83–100%.
- **The calls cluster in 9 bp.** 164–172 holds most of them, recurring across all four replicates at 1.7–2.5% of reads.
- **They are not a haplotype.** Each read carries at most one of them, so read phasing does not see a co-occurring block (`no co-occurring block of >= 4 sites`) and correctly leaves them alone.
- **Those libraries have a locally elevated error rate.** Mismatches per read over 164–173: 0.159 and 0.153 in R174ER178E r1 and r2, against 0.081 in EV r1 and 0.055 in AID r1, at the same median base quality.

So this is a library-specific artefact, not biology and not a general property of the run.

**The lesson for shmtrace.** A control-fitted null measures the background the control libraries share with the sample. It cannot see an artefact present in one library and absent from the controls, and no amount of calibration will make it. Two things catch this one, and both are already reported:

1. **The substitution spectrum.** A real AID signal in this background is C→T and G→A. Anything else is suspect. In a wild-type background the test is weaker but still useful: G:C bias and WRC/GYW enrichment.
2. **Positional clustering.** Calls packed into a few bases, recurring across replicates but not co-occurring within reads, are an error cluster rather than SHM.

**Both are now automatic.** This dataset is what prompted them:

- `clustered_calls` and a per-sample warning name any window of 10 bases holding far more calls than a random scatter would give (Poisson tail, Bonferroni over windows, p < 10⁻³). It fires on R174ER178E r1 and r4 (9 calls in 164–172) and on none of the AID, R50E or EV replicates. On 3,000 random scatters of 25 calls over 529 bp it never fires. Two alternative bases called at one position count as two calls, since that is itself a sign of an artefact.
- `pct_CT_GA` reports C→T and G→A as a share of all calls: 100% in every AID replicate, 0–29% in the mutant ones.

## At a matched 1,400 reads

Repeating with every library cut to 1,400 reads, so the two shallow mutant libraries are no longer the limiting factor, gives the same picture: AID 0.130 per 10³ bp (P = 0.015), R174ER178E 0.664 (P = 0.0038), R50E 0.000, EV 0.000 (no calls in any replicate). The conclusions do not depend on the depth chosen.

## Limits

- **A 229-bp window** of one VDJ amplicon in a cell line, not a polyclonal mouse assay.
- **A repair-deficient background.** The narrow C→T/G→A spectrum is what makes the check sharp, and it is also why the spectrum here cannot be compared with a wild-type one.
- **No deduplication.** As in `ILLUMINA.md`, distinct-mutation counts and the single-read estimate are unreliable on amplicon data this deeply resequenced; the single-read estimates in this run swing with depth and are not interpreted. Mutation frequency is the quantity used.
- **The mutant conditions have one sequencing run each** and 372–1,876 reads, against two runs and 5,280–14,679 reads for AID and EV.

## Still open

- `metrics.spectrum_note` scores calls against an arbitrary set of expected substitution classes, but nothing calls it yet; `pct_CT_GA` covers the common case.

## Reproduce

```bash
python /Users/abdulmateenaderinto/Nanoshm/public_data/run_prjna918596.py 5000
python /Users/abdulmateenaderinto/Nanoshm/public_data/run_prjna918596.py 1400
```

The download, MD5 check, per-replicate pooling and reference-building steps sit with the data, outside the repository.
