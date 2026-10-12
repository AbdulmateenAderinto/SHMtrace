# Calibration and validation

How the defaults were chosen, and how they behave on data they were not chosen on.
Every page names its data, its date and how to reproduce it.

## What the caller does, and why

| | |
|---|---|
| [`ERROR_MODELS.md`](ERROR_MODELS.md) | The error models, and why `recommended` means the quality-calibrated one. Spike-ins and a benchmark on J꜀4 and N114 |
| [`GENEIOUS.md`](GENEIOUS.md) | Reproducing the manuscript's Geneious calls, site by site |
| [`SINGLE_READ.md`](SINGLE_READ.md) | Mutations carried by one read, estimated in aggregate against the control background |
| [`INDELS.md`](INDELS.md) | Opt-in indel calling, and what it can and cannot see |
| [`CLONES.md`](CLONES.md) | Counting B-cell clones rather than reads, with the V(D)J junction as a barcode |
| [`TARGETING.md`](TARGETING.md) | Whether called mutations land where a 5-mer AID targeting model expects; a check independent of the error model |
| [`FINDINGS.md`](FINDINGS.md) | The July 2026 investigation that produced phasing, FDR correction and the control-fitted null. A historical record |
| [`CONFOUNDING.md`](CONFOUNDING.md) | What the J꜀4 result could still be confounded by |

## Validation on data the model was not tuned on

| | |
|---|---|
| [`EXTERNAL_VALIDATION.md`](EXTERNAL_VALIDATION.md) | A public nanopore deposit from another laboratory (Sμ, *Pim1*), and the GHVPDN run |
| [`ILLUMINA.md`](ILLUMINA.md) | An Illumina J꜀4 study: AID dependence recovered, and no LDHA effect, both as published |
| [`ECTOPIC_AID.md`](ECTOPIC_AID.md) | Ectopic AID versus empty vector in an AID-null line — the public counterpart of the N114 design |
| [`HMCES_INDELS.md`](HMCES_INDELS.md) | The indel caller against a published biological effect: HMCES loss raises SHM deletions |
| [`SANGER.md`](SANGER.md) | Matched Sanger clones from the same animals as the nanopore libraries |
| [`PUBLIC_DATA.md`](PUBLIC_DATA.md) | The survey of public SHM data with usable controls, and what is still missing |

## Scripts

`build_pileups.py`, `calibrate_loo.py`, `haplotype_check.py`, `strand_check.py`,
`confounding_check.py` reproduce the July figures. `single_read_spike.py` and
`indel_spike.py` are the spike-in harnesses for the two newer callers. Paths at
the top of each script point at the data and may need editing.
