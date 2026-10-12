# Sanger clones against nanopore, same animals

*11 October 2026. Sanger: the lab's J꜀4 clone set ("JH4 MKLN sanger sequencing data", Benchling GenBank exports, 762 files). Nanopore: the six bundled J꜀4 libraries (run DHJGKB). The Sanger PCRs were separate reactions from the nanopore ones, so clone-level identity is not expected. MKLN animals have no nanopore libraries and are left out.*

## Preparing the clones

- nanoshm now reads GenBank files (`.gb`, `.gbk`) directly. Each record is one read with flat Q30 qualities. IUPAC ambiguity codes are kept, so they are never counted as mutations.
- The 762 files hold 407 distinct clone records. The per-animal folders mix animals (the W2 folder holds w3 clones, and K3 holds JH4K1 clones), and folder `24` is the union of all of them. Clones were therefore assigned by name: `w2`/`wt2` → WT2, `jhextra-w3` → WT3, `JH4K1` → KO1, `k2` → KO2, and so on. A re-sequenced clone keeps its longest read.
- Matched animals and clones: WT2 20, WT3 38, WT4 28, KO1 25, KO2 32, KO3 25, KO4 35.

## A paralogue in the J꜀4 amplicon

46 of these clones (18%) are one sequence, about 70 differences (11%) from the J꜀4 intron over 614 bp, and 78 from the deposit's `Contaminant_JH4`. The clones agree with each other to within 0–11 differences, with or without AID, so this is not hypermutation. It is a J꜀4-like sequence that the J꜀4 primers also amplify. Its majority consensus is now the second record of the bundled `mouse_jh4_contaminant.fasta`. Competitive alignment drops it, and heavily mutated genuine J꜀4 clones stay far closer to J꜀4: the most mutated, at 117 differences from J꜀4, is 159 from the paralogue.

The deposited nanopore libraries contain none of it; all 497 reads are kept in each. Raw nanopore reads do contain it. In the FAM72A runs it is 0–41% of a library (CGQZ98 K2 41%, K32XSR W2 38%). Left in, it read as heavily mutated J꜀4: K32XSR W2's mutation frequency came out at 54.6 × 10⁻³ against 10.5 × 10⁻³ once it is removed.

## Comparison

Sanger mode (every substitution in a clone is a mutation, Q20 on flat Q30 input), with both contaminants removed. Nanopore in the recommended mode at 497 reads.

| Animal | Sanger clones | Sanger mutations / 10³ bp | Nanopore mutations / 10³ bp | Nanopore combined density × 10⁻⁵ |
|---|---:|---:|---:|---:|
| WT2 | 20 | 7.8 | 1.7 | 74.9 |
| WT3 | 34 | 29.1 | 21.3 | 138.1 |
| WT4 | 18 | 8.8 | 4.8 | 26.0 |
| KO1 | 15 | 3.4 | 0.0 | 0.5 |
| KO2 | 11 | 4.6 | 0.1 | 1.1 |
| KO3 | 23 | 1.4 | 0.02 | 0.4 |
| KO4 | 34 | 2.3 | – | – |

The clone counts are after the paralogue and contaminant are removed: KO2 drops from 31 to 11.

- **The ranking agrees.** WT3 is the most mutated animal on both platforms, and WT2 and WT4 sit well below it.
- **The substitution spectrum agrees.** Pooled across the three wild types (445 Sanger and 548 nanopore unique mutations), the 12-class spectra correlate with Pearson r = 0.91 (Spearman 0.84). WT3 alone gives 0.92. WT2 and WT4 alone give r ≈ 0.45, with 66 and 77 Sanger mutations, too few for a 12-class spectrum.
- **The Sanger knockout background is about 50 times nanopore's**: 1.4–4.6 against 0.02–0.13 × 10⁻³. These GenBank exports carry no qualities, so poor bases at the end of a read cannot be trimmed. Mismatches in knockout clones are enriched in the last 100 bases of the read (57 of 142), at reference positions 125–175. Some recur at wild-type hotspots (143, 160, 165 and 166, mostly G>A and C>T). That pattern suggests some wild-type material in the knockout clone preparations, not sequencing error alone. Worth checking at the bench.

## Indels

Sanger finds indels in all three wild types. WT3 has 13, including deletions of 22, 60 and 89 bp and a 26-bp insertion. WT2 has 3 and WT4 has 1. The knockouts have 0–2, mostly 1-bp. Nanopore calls 10 indels in WT3 and none in WT2 or WT4. Sanger counts each indel in a clone, while nanopore calls only indels carried by three or more reads (`INDELS.md`), so the two agree that WT3 is the indel-rich animal. They do not share individual events, as expected from separate PCRs and 34 clones against 497 reads.

Indels as a share of mutations: Sanger WT3 4.3%, nanopore WT3 3.0%.

## What this needs next

- Raw `.ab1` traces for these clones would give real qualities. Clone ends could then be trimmed, which would probably remove most of the Sanger knockout background.
- Whether the knockout clone sets carry some wild-type material. The recurrent hotspot mutations in the knockouts are the thing to check.
