# Dataset eligibility and completion ledger

Opened 17 October 2026 for the all-usable-datasets pass, on branch `public-validation-2026-10`. Status is one
of **done** (analysed, outputs committed), **eligible** (checks passed, not yet run), **blocked** (a documented
data or access barrier) or **unverified** (metadata not yet resolved). Sizes are compressed archive sizes from
the saved ENA/SRA tables, not biological sample sizes.

No `AGENTS.md` exists in either the workspace or the repository; none was applied.

## Done

| Dataset | Accession | What it tested | Outputs |
|---|---|---|---|
| Vincendeau 2022, PacBio CSR | PRJNA831666, 12 runs / 58 MB | the junction module; XLF-dependent microhomology recovered | `crm/out/switch_pacbio/`, Figure S16 |
| Compagno 2017, Illumina Sμ | GSE77788 / PRJNA311552, 18 runs / 540 MB | control calibration on an independent lab and platform | `crm/out/public_smu/`, Figure S17 |
| Germline-policy demonstration | synthetic, no download | what the 50% rule does to selected variants | `crm/out/germline_policy/`, below |

Both are written up in `PUBLIC_VALIDATION.md`, including why the GSE77788 drug effect is not claimed.

## Eligible, not yet run — engineered AID diversification

Execution order follows the prompt. Each needs a download, which has not been requested yet.

| # | Dataset | Accession / size | Planned endpoint | Resolved before running |
|---|---|---|---|---|
| E1 | Bassan 2021, MAGE-A1 TCR, PacBio | figshare 12803501, 12 mapped SAM / ~448 MB | linked variants across a ~2-kb engineered construct | the germline policy below; the exact construct reference; whether SAMs are filtered consensus; whether primer UMIs survive |
| E2 | Liu 2018, diversifying base editors | SRP159221, 193 paired runs / ~1.02 GB | AID against other mutators on comparable substrate | which labels are true biological replicates; the "NC" construct from the paper; amplicon overlap |
| E3 | Chen 2020, TRACE | PRJNA555784, 55 runs / ~2.36 GB | mutation accumulation, barcodes, MEK1 selection | editor/library mapping; the filtered mapped-BAM submission; orphan reads in 3-FASTQ exports |
| E4 | Base-editor roadmap 2021 | PRJNA649296, 84 runs / ~35.4 GB | CRISPR-X targeting controls, B2AR selection | which libraries sequence the target and which count guides |
| E5 | Devilder 2019, anti-HLA-A2 | supplement DOCX only | clone-level analysis against a parental construct | the supplement itself (not yet downloaded) |
| E6 | Ma 2016, TAM | SRP075898, 28 runs / ~19.3 GB | guide-dependent localisation, repair context | the mock export's mixed file roles; BCR-ABL baseline vs selected |

**E1–E5 must not be run before the germline policy is settled** (prompt requirement, and the demonstration
below shows why).

## Eligible, not yet run — SHM/CSR candidates

| Dataset | Accession / size | Note |
|---|---|---|
| Lescale 2025, PacBio CH12 CSR | PRJNA1235282, 29 runs / 3.83 GB | same group as PRJNA831666; resolve the nonproductive CH12 allele and repeat runs |
| Marton 2025, PacBio CSR | PRJNA1304899, 39 runs / 6.42 GB | same group again; published analysis code at Zenodo 17266277 |
| Álvarez-Prado 2018, capture | GSE102944, 9 libraries / ~22 GB | one pooled AID-null library; needs a capture-aware preparation path, which does not exist yet |
| Xu 2022, RAD52/IgD | PRJNA789895, 23 runs | Illumina junctions despite a PacBio-sounding description; verify the sequenced controls |
| SWIBRID 2025 | PRJNA1190672, 109 ONT runs | a junction comparator; re-analyses Vincendeau, so not an independent cohort |
| Wu 2026, inducible AID | PRJNA1337487, 659 AMPLICON runs / ~50 GB | run titles do not identify a replicated AID-off series; needs BioSample resolution first |

## Blocked, with evidence

| Dataset | Why |
|---|---|
| Singh 2019, MMTV/AID | Sanger sequences are request-only. MN126120 is a reference clone; GSE134189 is not the same collection. No author contact was made. |
| Bassan 2019, Pmel-1 TCR | The stated figshare deposit holds six TIFFs, not reads (`Pmel1_figshare.json`). |
| Hess 2016, original CRISPR-X | No raw archive verified. An accession will not be guessed. |
| Kling-EVOLVE 2025 | Raw data by request; only selected-antibody supplements are public. |
| ASNA1/Retro-2 | Processed Figure 4 source data only; raw reads unverified. Guide counts are not mutation reads. |
| Yeast yDBE 2023 | No raw accession verified; PDF supplement only. |
| Controlled AID/DT40 2025 | PQ885560–PQ885571 are selected antibody sequences; the full positive/control panel is not established. Gene conversion also confounds an SHM-only reading. |
| Ramos mUb-PCNA 2022, human AID-deficiency Sanger | Clone/trace availability not verified. |

## Out of scope by mechanism

TRIDENT (PRJNA701053) uses PmCDA1, not mammalian AID; AHEAD/OrthoRep likewise. These are optional orthogonal
technical comparisons, not AID validation, and are not counted toward the AID dataset set.

## The germline policy: demonstrated behaviour

`crm/germline_policy_demo.py` builds a synthetic library on a known ancestral reference, with acquired
variants at controlled frequencies and three clean control libraries, and runs the production project
workflow. Results (`crm/out/germline_policy/demo.json`):

| Acquired variant | default (0.50) | `germline_freq_threshold = 1.01` |
|---|---|---|
| 20% of reads | recovered | recovered |
| 50% | **removed as germline** | recovered |
| 90% | **removed as germline** | recovered |
| 95%, 98%, 99.5% | (not run) | recovered |
| 100% (fixed) | **not called** | **not called** |

Two separate behaviours, not one:

1. **The frequency rule** removes acquired variants from 50% up. This is at least *visible*: the removed sites
   are reported per sample as `high_fraction_sites`, with position, allele, reads and depth. For a selected
   population this erases exactly the variant selection enriched.
2. **Exact fixation is invisible.** A variant carried by every read is never called, under either setting, and
   does **not** appear in the removed list. Raising the threshold does not recover it. This is the more
   serious of the two, because nothing in the output says a position was dropped.

Raising the threshold is also not a safe general interface: `null_model.py:416` uses the same parameter to
exclude fixed differences from the error model, so a value above 1 would feed genuine fixed differences into
the background estimate. The threshold is read in nine places (`config.py:78`, `variants.py:346` and `:516`,
`indels.py:177` and `:218`, `project.py:223` and `:247`, `null_model.py:416`, `phasing.py:84`,
`lineage_sim.py:372`), so it is not a single-site switch.

`germline_vcf.py` does not solve this: it *subtracts* known germline sites supplied as a VCF. There is no
existing configuration that says "this reference is the known ancestor, so do not infer germline from
frequency". The only existing opt-out, `germline_on` in `variants.py:319`, is reachable only by the Sanger
model with too few clones.

**No code change has been made.** The smallest justified change, for the authors to approve, is an opt-in
`ancestral_reference` policy that (a) turns the frequency rule off for calling and phasing, (b) leaves the
null model excluding fixed differences as it does now, and (c) makes fixation callable and reported rather
than silent. It needs regression coverage first: the demonstration above is the starting point, plus a
baseline reference mismatch, partial coverage and reverse-complement orientation, and proof that the
physiological presets are unchanged.
