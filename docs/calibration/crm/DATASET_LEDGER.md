# Dataset eligibility and completion ledger

Opened 17 October 2026 for the all-usable-datasets pass, on branch `public-validation-2026-10`. Status is one
of **done** (analysed, outputs committed), **eligible** (checks passed, not yet run), **blocked** (a documented
data or access barrier) or **unverified** (metadata not yet resolved). Sizes are compressed archive sizes from
the saved ENA/SRA tables, not biological sample sizes.

No `AGENTS.md` exists in either the workspace or the repository; none was applied.

## Done

| Dataset | Accession | What it tested | Outputs |
|---|---|---|---|
| Vincendeau 2022, PacBio CSR | PRJNA831666, 12 runs / 58 MB | the junction module; XLF-dependent microhomology recovered | `crm/out/switch_pacbio/`, Figure S19 |
| Compagno 2017, Illumina Sμ | GSE77788 / PRJNA311552, 18 runs / 540 MB | control calibration on an independent lab and platform | `crm/out/public_smu/`, Figure S18 |
| Germline-policy demonstration | synthetic, no download | what the 50% rule does to selected variants | `crm/out/germline_policy/`, below |
| E1 Bassan 2021, AID-diversified TCR | figshare 12803501, 12 SAM / 228 MB | variant recovery against a known ancestral construct | `crm/out/tcr/`, below |
| E5 Devilder 2019, anti-HLA-A2 clones | supplement DOCX, 53 KB | calls against an independent published annotation | `crm/out/devilder/`, below |
| E2 Liu 2018, diversifying base editors | SRP159221, 386 files / 1.02 GB | whether the spectrum separates AID from an APOBEC | `crm/out/liu/`, below |
| E6 Ma 2016, TAM reporter | SRP075898 amplicon arms, 3.5 GB | guide-dependent localisation | `crm/out/tam/`, below |

Both are written up in `PUBLIC_VALIDATION.md`, including why the GSE77788 drug effect is not claimed.

## E1, Bassan 2021: what it showed

Reference `hT27` reconstructed exactly from the deposit's CIGAR+MD tags: 1,785 bp, matching the header, with
zero disagreements between reads and no ambiguous base. The headers also answer two of the prompt's
questions: the reads are CCS aligned with `bwa mem -x pacbio`, so they are not subreads, and read names are
`movie/zmw/ccs`, which identify a ZMW and not a molecular UMI.

Eleven libraries were usable at 6,000 reads each. `11_h12_HA` holds 19 reads and is excluded by count.

| sample | unique positions | mutation frequency /10³ | C→T or G→A |
|---|---|---|---|
| **h5_wt (the available baseline)** | **182** | **0.079** | 61% |
| h5_MHA / h5_HA / h5_HAHA | 39 / 9 / 11 | 0.56 / 0.59 / 0.59 | 36 / 22 / 27% |
| h7_HA / h7_HAHA | 28 / 47 | 1.44 / 1.77 | 32 / 23% |
| h8_MHA / h8_HA / h8_HAHE | 66 / 54 / 57 | 1.10 / 0.81 / 0.94 | 35 / 33 / 26% |
| h12_MHA / h12_HAHE | 47 / 49 | 0.69 / 0.79 | 36 / 53% |

The baseline carries the **most** positions and by far the **lowest** frequency: many scattered rare
differences, against 7–23× higher mutation frequency at fewer positions in the diversified populations. That
is what selection on a diversified pool looks like, and it is the comparison this deposit can support.

What it cannot support: there is one baseline library and it belongs to founder h5, so no false-discovery
rate is claimed and `1_h5_wt` is not treated as a matched AID-off control for h7, h8 or h12. Selected
endpoint frequencies are not intrinsic mutation rates.

**Policy sensitivity, measured rather than assumed.** Re-running three libraries with
`ancestral_reference` on and off changed one of them: h7_HAHA gave 29 called positions with the policy and 27
without, the two differences being variants the frequency rule removed as germline. The other libraries were
identical. So on this deposit the policy recovers real variants but the effect is small, because most
selected variants here sit below 50%. It remains the right default for this class of experiment, and the
claim is kept to what was measured.

## E5, Devilder 2019: concordance with a published annotation

Additional file 1 publishes the parental heavy variable region (423 bp) and, for each selected cluster of
rounds 1 and 2, its nucleotide sequence *and* the nucleotide mutations the authors annotated. That annotation
is an independent answer to check calls against, on a defined reference, with no read processing in between.

| round | clusters | with a published annotation | reproduced exactly | published but not recovered | recovered but not published |
|---|---|---|---|---|---|
| R1 | 17 | 17 | **17** | 0 | 0 |
| R2 | 26 | 25 | **25** | 0 | 0 |

All 42 annotated clusters reproduce the authors' mutation lists exactly, in both directions. Through
SHMtrace's own clone path the sample gives 32 unique mutations with `ancestral_reference` and 31 without:
the one difference is a mutation carried by at least half the selected clusters, which the frequency rule
removes as germline. That is the smallest clean illustration of why the policy exists.

What this does not test: it is a curated collection of *selected* clones, so it carries no population
frequencies, no unselected distribution, no replication and no technical error rate, and none are computed.
It exercises reference handling and substitution identification, not the error model or the false-discovery
calibration.

## E2, Liu 2018: the spectrum separates the mutators, the transition share does not

The deposit edits one 238-bp EGFP amplicon with a panel of deaminases delivered the same way, three
replicates each, plus a no-editor arm. Pairs merge cleanly (97.8%), so the ordinary Illumina path applies.
The contrast was fixed before looking at any result: DBE-AIDmono against DBE-humanA3a, in wild-type and in
UNG/MSH2-null 293T, with NC as the control. Motif shares are descriptive and were never used to choose what
to call.

| background | measure | AID (3 reps) | A3A (3 reps) | Welch p |
|---|---|---|---|---|
| wild-type 293T | **WRC/GYW share** | **56.5%** (46.9–65.2) | **14.4%** (7.1–21.6) | **0.004** |
| wild-type 293T | C→T or G→A | 44.7% (39.8–52.1) | 17.9% (13.3–20.3) | 0.007 |
| wild-type 293T | frequency /10³ | 0.66 | 0.56 | 0.32 |
| UNG/MSH2-null | **WRC/GYW share** | **37.1%** (35.2–39.6) | **20.0%** (17.3–22.0) | **0.0009** |
| UNG/MSH2-null | C→T or G→A | 78.1% (71.3–83.5) | 75.6% (64.0–84.1) | 0.74 |
| UNG/MSH2-null | frequency /10³ | 3.8 | 3.6 | 0.77 |

The no-editor arm gave **0, 0 and 1** calls in 20,000 reads each, so the background control holds on an
independent laboratory's construct.

Two things worth taking from this. The **motif preference separates the two deaminases in both backgrounds**,
with non-overlapping replicate ranges. The **transition share does not**: it separates them in wild-type
cells but collapses in the repair-null background, where both sit near 78%, because losing uracil excision
and mismatch repair drives C→T regardless of which enzyme made the uracil. The same loss raises the mutation
frequency about five-fold for both mutators, so frequency identifies the repair background, not the enzyme.
A spectrum check that rested on the transition share alone would therefore have mistaken an APOBEC for AID in
repair-deficient cells; the motif check would not.

Depth caveat: three of the twelve editor libraries supplied fewer than 20,000 usable reads (4,792–12,518),
which is recorded per sample in `liu_mutators.json`.

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

## E6, Ma 2016 (TAM): guide-dependent localisation, and a limitation confirmed

The reporter arm has four libraries - no editor (mock), AID alone, dCas9-AIDX with one sgRNA against EGFP,
and dCas9-AIDX with a pool - one per condition, so nothing is replicated and no group statistic is computed.

Two things had to be resolved before any number meant anything:

- **The mock is not what its label says.** ENA calls `SRR3608504` SINGLE, but it exports three FASTQ files:
  spots 1-809,813 are a paired block and 809,814-1,292,558 are single-end, together matching the run's
  stated 1,292,558 reads. The other three libraries are single-end only, so the single-end block is used and
  the paired block excluded, rather than mixing two preparations.
- **The libraries are multiplexed.** Five amplicons share each file, each with its own start. A single
  position-wise consensus over a whole library averages different sequences together - it gave a worst-case
  consensus of 27.8%, which is what first exposed the problem. Built per amplicon instead, four of the five
  reach 86-93% and one reaches 62%, which is flagged rather than used confidently.

Mismatches per 1,000 bases at Q ≥ 30, against each amplicon's own reference:

| library | ACGTCTATAT | CCCCTTGAAC | **CCTGACCTAC** | GCCCTATGTC | TTTTCTCCCA |
|---|---|---|---|---|---|
| mock (no editor) | 9.5 | 2.6 | **5.5** | 1.5 | 3.0 |
| AID alone | 10.5 | 3.1 | **4.2** | 2.1 | 17.3 |
| dCas9-AIDX + sgRNA | 20.4 | 2.1 | **22.6** | 2.0 | 17.1 |
| dCas9-AIDX + pool | 11.1 | 2.8 | **20.6** | 1.3 | 30.3 |

At `CCTGACCTAC` the two **tethered** arms raise the rate about four-fold over mock while **untethered AID
does not** (4.2 against mock's 5.5). That is guide dependence, not merely editor dependence, and two
amplicons (`CCCCTTGAAC`, `GCCCTATGTC`) stay flat in every arm as internal negatives. The `TTTTCTCCCA`
column is raised by every editor including untethered AID, and is also the amplicon whose consensus is
weakest, so it is not relied on.

**Why these are raw rates and not calls.** The arm has a single no-editor library, and SHMtrace fits its
error model from the controls, so with one control there is nothing to calibrate against. Run anyway on the
targeted amplicon at 15,000 reads apiece, the calibrated path does not merely lose power, it inverts:

| library | raw rate /10³ | calibrated unique calls |
|---|---|---|
| mock (the control) | 5.5 | **117** |
| AID alone | 4.2 | 6 |
| dCas9-AIDX + sgRNA | **22.6** | **0** |
| dCas9-AIDX + pool | 20.6 | 46 |

The arm with the highest raw excess returns no calls, and the control itself accumulates the most. That is
the manuscript's own requirement for two or more AID-null controls, confirmed on an independent dataset
rather than asserted, and it is why only the raw rates above are interpreted here. Those rates still contain
sequencing error and are descriptive, not mutation calls.

### The ABL exon-6 footprint arm: the position follows the guide

Five libraries target the same locus with five different sgRNAs (#4, #6, #7, #9, #15), which makes the arm
self-controlled for localisation: each library is the negative control for the others at their own peaks, so
no no-editor library is needed and none exists here. This arm is multiplexed too - three amplicons share
each file - and analysing only the commonest one showed a footprint for a single guide. All three are
analysed below; the mistake is recorded because it is easy to repeat.

Peak interior position per library, with the median of the other four libraries at that same position
(mismatches per 1,000 bases at Q ≥ 30, 20,000 merged pairs per library):

| guide | amplicon | peak position | rate | others there | fold |
|---|---|---|---|---|---|
| sgRNA_9 | amp1 | 154 | 111.1 | 0.30 | **369×** |
| sgRNA_15 | amp1 | 94 | 136.4 | 2.26 | **61×** |
| sgRNA_6 | amp1 | 95 | 141.5 | 4.61 | **31×** |
| sgRNA_4 | amp1 | 49 | 5.3 | 0.61 | **8.6×** |
| sgRNA_7 | amp1 | 98 | 43.6 | 34.8 | 1.3× |
| sgRNA_9 | amp0 | 38 | 13.9 | 0.36 | 39× |
| any guide | amp2 | 228-234 | 2.4-3.4 | ~3.0 | 0.8-1.2× |

**Three of the five** - #6, #9 and #15 - give an unambiguous localised peak at a distinct position, each
reaching 116-137 mismatches per 1,000 against a background near 1, and 31-369 times the other libraries at
that same base. The third amplicon shows nothing in any library, as an internal negative.

The other two need stating carefully, and the fold column alone would mislead about both:

- **sgRNA_4** passes the fold test at 8.6x, but only because the other libraries happen to be very low at
  position 49 (0.61). Its own peak reaches just 6 per 1,000 on a broad, noisy profile, so it is weak
  evidence of a footprint, not a clean one. Figure S18B shows this directly.
- **sgRNA_7** reaches 43 per 1,000 but scores only 1.3x, because its peak at 98 sits inside the region where
  #6 (95) and #15 (94) are also high, so the comparison term is itself elevated. The fold statistic assumes
  guides target different places and understates for overlapping targets.

So the honest count is three clean footprints, one overlapping and one weak - which is still a position that
follows the guide, but not five for five.

### The saturation arm: more guides, more localised sites

Four libraries carry 3, 5, 7 and 10 sgRNAs against the same locus. Counting positions where a library sits at
least five-fold above the lowest of the **other three** at that same position (and at least 2 per 1,000):

| guides carried | localised sites |
|---|---|
| 3 | 0 |
| 5 | 32 |
| 7 | 37 |
| 10 | 66 |

The count rises with the number of guides. The leave-one-out baseline matters and was checked: scored
against the lowest of all four, the library with the fewest guides is usually the minimum and would be
compared with itself. Recomputed leaving each library out of its own baseline, the totals are unchanged, so
the zero at three guides is not an artifact of that choice.

Two caveats. A zero is not proof of no editing - those three guides may simply target outside the three
amplicons recovered here. And the reference is the pooled consensus of the same libraries being measured, so
heavy editing degrades it: the amplicon carrying most of the signal has a weakest-position consensus of
84.8% here, against 96.9% in the less heavily edited footprint arm. The counts are therefore a trend, not a
calibrated yield.

### The reporter-footprint pair (single sgRNA ± UGI): not usable for the repair contrast

This looked like a chance to test, in a second system, the repair-context effect E2 showed - that losing
uracil excision drives the transition share up whichever enzyme made the uracil. It does not survive
inspection. The two libraries differ in three ways at once:

| | no UGI (`SRR3620780`) | + UGI (`SRR3620777`) |
|---|---|---|
| layout | single-end, one file | paired, two files |
| read length | 300 bp | 251 bp |
| dominant amplicons | `ACGTCTATAT`, `CCTGACCTAC`, `CCCCTTGAAC` | `GCCCTATGTC` (96% of reads) |

They share `GCCCTATGTC`, so a comparison could be forced onto it using R1 of the paired library over the
common span. But with one library per condition, and UGI confounded with layout, chemistry and amplicon
composition, any spectrum difference would have three other explanations. No number is reported from it.

The project's off-target arm has two mock replicates and would support a calibrated analysis, but it is
genome-wide fragmented material and would have to be forced through amplicon-shaped gates. It was
deliberately not downloaded (17 of the 19.3 GB); a fragmented-read route is a separate capability.

## E3, Chen 2020 (TRACE): downloaded, then blocked on metadata

All 110 files were fetched (1.1 GB) and inspected before this was called. The blocker is not the reads.

The archive resolves the **assay** but not the **condition**. Every run's `sample_alias` is a figure label
(`FigS5`, `Fig2A_and_2B`, `Fig2`, `FigS3`, `Fig1H`), `sample_title` is "Human sample from Homo sapiens", and
`library_name` is a bare index (1, 3, 6, 7, 10, 13…). Nothing names the editor, the target, the induction,
the drug or the replicate. The paper's data-availability statement gives only the BioProject, and the saved
text carries no table mapping a library index to a condition. Assignments were not guessed.

What the reads themselves do settle, and it is worth recording:

- The groups sequence **different constructs**: the MEK1-screen reads carry the T7 promoter
  (`TAATACGACTCACTATAG`), the Fig2A_and_2B reads are lentiviral backbone, the off-target reads are genomic.
- Runs exporting **three FASTQ files** are not three replicates. In `SRR9722420` the `_1`/`_2` files hold
  147-bp pairs and the unsuffixed file holds 302-bp reads with *different* spot ids, so the third file is
  unpaired singletons.

To unblock: a library-index-to-condition table from the authors or a supplementary file that carries one.
Without it the deposit supports no prespecified contrast, and inventing one would be worse than leaving it.

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
