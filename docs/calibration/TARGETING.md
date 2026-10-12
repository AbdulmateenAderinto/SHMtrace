# AID targeting against a 5-mer model

*12 October 2026. Model: SHazaM's HH_S5F (Yaari et al., Front Immunol 2013), exported from a local SHazaM 1.2.0 install with `scripts/export_shazam_model.R`. Data: the six bundled J꜀4 libraries, and the CH12 ectopic-AID set from [`ECTOPIC_AID.md`](ECTOPIC_AID.md).*

## Why

The hotspot figures put each mutation in a binary motif class (WRC/GYW, WA/TW) and count. A 5-mer targeting model gives every 5-mer a continuous mutability, learned from about 800,000 silent mutations in human repertoires, so every position of the amplicon gets an expected share of AID's mutations. That allows a test that does not depend on the error model at all: **do the called mutations land where AID is expected to land?** Real AID mutations should. Sequencing errors, which know nothing of AID, should not.

A model fitted on the amplicon itself would be useless. Each 5-mer occurs about once in a 529-bp window, so it would reproduce the positional map. The expectation has to come from other data.

## The statistic

The **targeting index** is the mean model mutability at the mutated positions divided by the mean over the window. Values near 1 mean the mutations ignore the model; values above 1 mean they favour where AID acts. The p-value is one-sided, from 2,000 random position sets of the same size. Fewer than 10 mutations give no index: one or four mutations on a hotspot say nothing.

## Does a human heavy-chain model work on a mouse intron?

That was the open question, since SHazaM has no mouse heavy-chain model. It does work.

**Model check.** In HH_S5F, 5-mers centred on the mutated base of a WRCY/RGYW hotspot average 10 times the mutability of SYC/GRS coldspots (0.00222 against 0.00022).

**J꜀4, AID mutations against sequencing error:**

| Mutations | n | Targeting index | p |
|---|---:|---:|---:|
| Wt2 calls | 146 | 1.34 | 0.0005 |
| Wt3 calls | 334 | 1.30 | 0.0005 |
| Wt4 calls | 68 | 1.19 | 0.037 |
| **Every mismatch position in the three knockouts** | 344 | **1.00** | 0.47 |
| Every mismatch event in the knockouts, weighted | 1,063 | 0.98 | 0.78 |

Wild-type calls favour the model's mutable positions; sequencing error in AID-null mice does not, to two decimal places. The dose-response is steeper than the index suggests. Splitting the window's positions into fifths by model mutability, wild-type unique mutations per position rise about fivefold from the least to the most mutable fifth, and the knockouts stay at zero throughout.

The observed wild-type spectrum and the spectrum the model predicts for this window correlate at r = 0.80.

**CH12, real AID against a known artefact** (the ectopic-AID set, where one AID mutant's calls are a library-specific error cluster):

| Calls | n | Targeting index | p |
|---|---:|---:|---:|
| AID | 26 | **2.55** | 0.0005 |
| R174E/R178E (the artefact) | 50 | **0.72** | 0.995 |

The index flags the artefact as non-AID by itself. That makes it a third check independent of the error model, alongside the spectrum (`pct_CT_GA`) and positional clustering (`clustered_calls`). The calls in that mutant are not merely off the model; they sit preferentially on positions AID disfavours.

## What is reported

With `targeting_model` set (the project key, or the model drop under Advanced in the app):

- `s5f_index` and `s5f_p` per sample, and per group under `targeting`;
- each call's `s5f_mutability_rel` (its 5-mer's mutability relative to the window) and `s5f_substitution` (the model's probability of that exact change);
- three figures: the dose-response with the model along the amplicon, observed against expected spectrum, and the model itself drawn as SHazaM's hedgehog plot.

Nothing here changes which mutations are called. The per-call annotation could become a likelihood term in the error model, and the CH12 result suggests it would help, but that needs its own spike-in validation first.

## Licensing

SHazaM is AGPL-3 and nanoshm is MIT, so nanoshm does not ship SHazaM's models. `scripts/export_shazam_model.R` writes one from your own installation:

```bash
Rscript scripts/export_shazam_model.R HH_S5F HH_S5F.csv
```

Cite Yaari et al. 2013 and SHazaM when reporting results that use it.

## Limits

- **A human heavy-chain model on a mouse intron.** It separates AID from error here, but it is not a mouse model, and the absolute index should not be compared across species or loci.
- **One window.** The index depends on how much the window's mutability varies; a flat window cannot produce a large index whatever the biology.
- **The single-read mutations are not scored.** The index uses called mutations only. Scoring the aggregate single-read excess by position is possible but not built.
