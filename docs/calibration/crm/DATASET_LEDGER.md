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

| Acquired variant, case only | default (0.50) | `germline_freq_threshold = 1.01` |
|---|---|---|
| 20% of reads | recovered | recovered |
| 50%, 90%, 95%, 98%, 99.5% | **removed as germline** | recovered |
| 100% (fixed) | **removed as germline** | recovered |

So the rule does exactly what it says, and it is **visible**: every removed site is reported per sample as
`high_fraction_sites`, with position, allele, reads, depth and frequency. For a selected population this
removes exactly the variant selection enriched, which is the problem for the engineered datasets.

A **baseline mismatch** — a fixed difference from the deposited reference carried by the controls as well as
the case — is correctly not called at either threshold, because the control-calibrated test sees the same
thing on both sides. The control calibration already handles that case; the germline rule is not what
protects against it.

*(An earlier version of this section reported that fixation was silently invisible. That was an artefact of
the demonstration script, which planted an alternate base that happened to equal the reference at that
position, so no variant existed there. The script now derives each alternate from the reference base. The
corrected behaviour is the table above: fixation is removed by the rule, and reported.)*

Raising the threshold above 1 does recover the acquired variants, and in the baseline-mismatch test it did not
inflate the control calls. It is still not the right interface: `germline_freq_threshold` is read in nine
places (`config.py:78`, `variants.py:346` and `:516`, `indels.py:177` and `:218`, `project.py:223` and `:247`,
`null_model.py:416`, `phasing.py:84`, `lineage_sim.py:372`), including the null model's exclusion of fixed
differences from the error estimate, so a value above 1 changes more than calling and says nothing about
intent. Whether that pooled-background effect is material has **not** been demonstrated.

`germline_vcf.py` does not solve this either: it *subtracts* known germline sites supplied as a VCF. There is
no existing configuration that says "this reference is the known ancestor, so do not infer germline from
frequency". The only existing opt-out, `germline_on` in `variants.py:319`, is reachable only by the Sanger
model with too few clones.

## The patch: `ancestral_reference`

Added at the authors' direction, tests first (`tests/test_ancestral_reference.py`, six tests). Default
`False`, so every preset and existing project is untouched.

| File | Original | Change |
|---|---|---|
| `config.py:79` | — | `ancestral_reference: bool = False`, with the reasoning in a comment |
| `variants.py:319` | `germline_on = True`, only the Sanger path could clear it | `germline_on = not p.ancestral_reference`, and a note is added to the sample when it is off |
| `phasing.py:84` | candidate fraction bounded above by `germline_freq_threshold` | bound becomes 2.0 under the policy, so a selected variant at fixation can still be phased |
| `null_model.py:416` | a control allele at or above the threshold is excluded from the error model | under the policy it is kept, as a position-specific rate |

The last one was not in the original design and was forced by a failing test. Leaving the null model alone
made a **baseline mismatch** — a fixed difference from the deposited reference that the controls carry too —
into a false positive: with the frequency rule off, nothing stopped the case being called for carrying what
its own controls carry. Keeping those control alleles in the error model restores the right answer, because
the controls are the un-diversified baseline. This is why
`test_baseline_mismatch_is_not_called_either_way` exists.

Verified behaviour (`crm/germline_policy_demo.py`, and the tests):

| | default (0.50) | `ancestral_reference` |
|---|---|---|
| acquired at 20% | recovered | recovered |
| acquired at 50–99.5% | removed, and reported | recovered |
| acquired at 100% | removed, and reported | recovered |
| baseline difference shared with controls | not called | not called |
| reverse-complement reads, partial coverage | — | recovered |
| physiological presets | unchanged | unchanged (option defaults off) |

**Scope deliberately left out.** `indels.py` keeps the frequency rule under the policy: the datasets this was
built for (AID-diversified receptors, base-editor substrates) are substitution experiments, and widening the
change to indels would need its own evidence. `project.py`'s Sanger germline and `lineage_sim.py` are
likewise untouched. Anyone enabling the policy for an indel-bearing experiment should read this first.
