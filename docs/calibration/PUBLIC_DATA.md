# Public data for validation beyond nanopore

*Surveyed 11 October 2026 (NCBI SRA and ENA). Aim: show that the error models hold on other platforms, against AID-null (or otherwise mutation-free) controls.*

## Found, with usable controls

| Study | Platform | Locus | Samples | Why it is useful |
|---|---|---|---|---|
| [PRJNA973112](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA973112), LDHA in germinal centres (GEO GSM7360896–913) | Illumina NovaSeq, 2×151 | **Mouse J꜀4 intron**, germinal-centre B cells | 7 WT, 7 LDHA-KO, **4 AID-KO**; about 12 M read pairs each | Same locus and assay as the paper, on Illumina, with AID-null controls. The direct Illumina counterpart of the nanopore J꜀4 calibration |
| [PRJNA820976](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA820976), HMCES protects Ig genes from deletions | Illumina MiSeq, 2×300 | Human IGHV in Ramos cells, plus a reporter line (RASH-1C) | 6 Ramos WT clones, 6 Ramos HMCES-KO clones, **3 Ramos AID-KO clones**; RASH-1C with UNG, MSH2, MSH6, EXO1, APE1 and APE2 knockouts | HMCES loss raises SHM **deletions**, a known real-indel effect against an AID-null control. The best public check of the indel caller. The reference needs building from the AID-KO consensus |
| [PRJNA918596](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA918596), mesoscale DNA feature | Illumina HiSeq X | VDJ in CH12F3 Ung-/-Msh2-/-Aicda-/- cells, B1-8 | AID versus **empty-vector (no-AID)** clones; 880 runs | V-region SHM on a repair-deficient, AID-null background, with matched no-AID controls |
| [PRJNA1481264](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1481264), MLLT1/MLLT3 | Oxford Nanopore | Sμ, Pim1 | WT and Aicda-/- | Already used (`EXTERNAL_VALIDATION.md`) |

## Found, weaker controls

| Study | Platform | Locus | Note |
|---|---|---|---|
| [PRJNA590097](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA590097), Rev7 | Illumina HiSeq X | J꜀4 and Jκ5 introns, Peyer's patch PNA-high and PNA-low | No AID-KO; PNA-low (non-germinal-centre) cells are a low-mutation comparison, not a null |
| [PRJNA1124255](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1124255), ELOF1 Rep-SHM-seq | Illumina NovaSeq | IGHV | Germinal-centre against non-germinal-centre B cells; 754 GB in all |
| [PRJEB52221](https://www.ebi.ac.uk/ena/browser/view/PRJEB52221), IgH Eμ-MARs | Ion Torrent S5, 125 bp single | Sμ | WT, MAR-deleted, UNG and MSH mutants; no AID-null. A fourth platform, but only a ranking |
| [PRJNA588067](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA588067), dietary fibre | Illumina MiSeq 2×300 | VH186.2 | TCR-KO comparison, not AID-null |

## Not found

- **Sanger.** No public Sanger SHM traces with AID-null controls turned up in SRA, and GenBank holds the J꜀4 reference (J00440, NG_005838), not clone sets. The JoVE J꜀4 protocol (Sible et al. 2021) reports Sanger counts (226 mutations in 105 WT clones, 2 in 113 AID-KO) but deposits no traces. The best Sanger calibration is the lab's own Sanger clones from the same animals as the nanopore runs.
- **PacBio** SHM amplicons with AID-null controls: none found.

## Order of work

1. ~~PRJNA973112 on J꜀4~~ done: `ILLUMINA.md`. Wild type is about 18-fold above *Aicda*-/- (P = 0.015). The amplicon covers 70 bp (J꜀4 segment and start of the intron), not the intron window.
2. PRJNA820976: indels (HMCES-KO against WT, both against AID-KO), and substitutions in a human V region.
3. ~~The lab's matched Sanger clones~~ done: `SANGER.md`.
4. PRJNA918596 if a V-region Illumina set is still needed after 2.
