# Sanger clones against nanopore, same animals

*11 October 2026. Sanger: the lab's J꜀4 clone set ("JH4 MKLN sanger sequencing data", Benchling GenBank exports, 762 files). Nanopore: the six bundled J꜀4 libraries (run DHJGKB). The Sanger PCRs were separate reactions from the nanopore ones, so clone-level identity is not expected. MKLN knockouts have no nanopore libraries; they are shown for reference only.*

## Preparing the clones

- nanoshm now reads GenBank files (`.gb`, `.gbk`) directly. Each record is one read with flat Q30 qualities. IUPAC ambiguity codes are kept, so they are never counted as mutations.
- The 762 files hold 407 distinct clone records. The per-animal folders mix animals (the W2 folder holds w3 clones, and K3 holds JH4K1 clones), and folder `24` is the union of all of them. Clones were therefore assigned by name: `w2`/`wt2` → WT2, `jhextra-w3` → WT3, `JH4K1` → KO1, `k2` → KO2, and so on. A re-sequenced clone keeps its longest read.

## A paralogue in the J꜀4 amplicon

46 of these clones (18%) are one sequence, about 70 differences (11%) from the J꜀4 intron over 614 bp, and 78 from the deposit's `Contaminant_JH4`. The clones agree with each other to within 0–11 differences, with or without AID, so this is not hypermutation. It is a J꜀4-like sequence that the J꜀4 primers also amplify. Its majority consensus is now the second record of the bundled `mouse_jh4_contaminant.fasta`. Competitive alignment drops it, and heavily mutated genuine J꜀4 clones stay far closer to J꜀4: the most mutated, at 117 differences from J꜀4, is 159 from the paralogue.

The deposited nanopore libraries contain none of it; all 497 reads are kept in each. Raw nanopore reads do contain it. In the FAM72A runs it is 0–41% of a library (CGQZ98 K2 41%, K32XSR W2 38%). Left in, it read as heavily mutated J꜀4: K32XSR W2's mutation frequency came out at 54.6 × 10⁻³ against 10.5 × 10⁻³ once it is removed.

## Which clones are which

Confirmed by Abdulmateen on 11 October 2026. The Sanger W2, W3 and W4 mice are nanopore WT2, WT3 and WT4. A1–A3 are the *Aicda*-/- mice. K1–K4 (also labelled `MKLN_KO`) are MKLN knockouts, which keep AID and have no nanopore libraries. The `jh6-` clones are another amplicon and are left out. Of the three `w1`-labelled clones, none duplicates another clone (each is 1–2 differences from its nearest wild-type clone). `w1-jh4` sits in the WT2 folder and is counted as WT2. `w1.1extra` and `w12extra` sit only in the W3 folder, which also holds clones from other mice, so they are left out. Sanger A1, A2 and A3 are nanopore KO1, KO2 and KO3.

After the paralogue and contaminant are removed: WT2 21 clones, WT3 34, WT4 18, A1 4, A2 5, A3 1, MKLN K2 11, K3 23, K4 34.

## Trimming the messy ends

Sanger reads are messy for their first and last stretch, and these exports carry no qualities to trim by. The manual practice is to align the clones and cut each one back to where it cleans up. nanoshm now does the same for the Sanger platform (`sanger_clean_window`, `sanger_clean_max_mismatch`). After alignment each clone is cut back, from each end, to the first run of 20 reference bases with no ambiguity code, no indel and at most one mismatch. Mutations inside the clean part are kept.

On these 150 clones the clean part starts at reference 50–75 and ends at 650–675, about where the insert meets the primers. That lies outside the 529-bp comparison window (101–629), which is why the manuscript chose that window. Inside the window, trimming changed only a handful of calls.

## Comparison

Sanger mode (every substitution in a clean, trimmed clone is a mutation, Q20 on flat Q30 input), with both contaminants removed. Nanopore in the recommended mode at 497 reads.

| Mouse | Sanger clones | Sanger unique mutations | Sanger mutations / 10³ bp | Nanopore mutations / 10³ bp | Nanopore combined density × 10⁻⁵ |
|---|---:|---:|---:|---:|---:|
| WT2 | 21 | 67 | 7.5 | 1.7 | 74.9 |
| WT3 | 34 | 290 | 28.2 | 21.3 | 138.1 |
| WT4 | 18 | 77 | 8.8 | 4.8 | 26.0 |
| *Aicda*-/- A1 (nanopore KO1) | 4 | 0 | 0 | 0.0 | 0.5 |
| *Aicda*-/- A2 (KO2) | 5 | 1 | 0.4 | 0.1 | 1.1 |
| *Aicda*-/- A3 (KO3) | 1 | 5 | 9.5 | 0.02 | 0.4 |
| MKLN-/- K2, K3, K4 | 11, 23, 34 | 18, 13, 30 | 4.6, 1.4, 2.3 | – | – |

- **AID dependence agrees.** Wild types rank the same on both platforms (WT3 well above WT2 and WT4). The *Aicda*-/- Sanger clones carry 1 mutation in A1 and A2's 9 clones. A3 is a single clone with 5 mutations, too little to read anything from. As the manuscript says, ~10 Sanger clones put the knockout background below the detection floor, while nanopore at ~500 reads measures it.
- **The substitution spectrum agrees.** Pooled across the three wild types (433 Sanger unique mutations, before `w1-jh4` was added to WT2):

  | Nanopore mode | Unique mutations, Pearson r | Recurrence-weighted, Pearson r | Mean absolute difference per class |
  |---|---:|---:|---:|
  | geneious (the manuscript's rule) | 0.947 | 0.962 | 1.3–1.4 points |
  | recommended | 0.918 | 0.963 | 1.4–1.5 points |

  The manuscript reports r = 0.97 and 1.4 points from Geneious with 515 Sanger mutations. The small gap is consistent with a slightly different clone set.
- **The MKLN knockouts are mutated** (1.4–4.6 × 10⁻³, with G:C mutations at WRC/GYW hotspots 40–79% of the time), as expected for an AID-proficient line. An earlier version of this page read them as AID-null controls and took their hotspot mutations for carry-over. That was wrong.

## Indels

Sanger finds indels in all three wild types. WT3 has 13, including deletions of 22, 60 and 89 bp and a 26-bp insertion. WT2 has 3 and WT4 has 1. The *Aicda*-/- clones have none; the MKLN knockouts have 0–2, mostly 1-bp. Nanopore calls 10 indels in WT3 and none in WT2 or WT4. Sanger counts each indel in a clone, while nanopore calls only indels carried by three or more reads (`INDELS.md`), so the two agree that WT3 is the indel-rich animal. They do not share individual events, as expected from separate PCRs and 34 clones against 497 reads.

Indels as a share of mutations: Sanger WT3 4.5% (13 of 290), nanopore WT3 3.0%.

## What this needs next

- Few *Aicda*-/- clones were Sanger-sequenced (10 across three mice). They are kept as the Sanger control, but they cannot calibrate a Sanger background; the nanopore knockouts do that.
- Raw `.ab1` traces would give real qualities in place of the flat Q30. That matters less now that clones are trimmed to their clean part.
