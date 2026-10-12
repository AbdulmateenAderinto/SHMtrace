# Illumina: J꜀4 in an independent study

*11 October 2026. Data: [PRJNA973112](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA973112) (GEO GSM7360896–913), germinal-centre B cells from NP-CGG-immunised mice, day 14. Illumina NovaSeq 2×151. 7 wild type, 7 LDHA-KO and 4 *Aicda*-/- mice. All 72 FASTQ files checked against ENA's MD5.*

## What this amplicon covers

The study amplified from a primer in the rearranged VDJ to one in the J꜀4 intron (Fwd `CATAACAAAGGTTAAAAATAAAGAC`, Rev `GCCTGACATCTGAGGACTC`) by reverse-complement PCR. On nanoshm's J꜀4 reference, the intron primer ends at position 123. The J꜀4 gene segment is 53–104, and the intron starts at 105.

The usable window is therefore **53–122 (70 bp)**: the J꜀4 coding segment and the first 18 bp of the intron. It is not the 529-bp intron window of the nanopore assay. Both mates of a pair cover all of it, so pairs are merged before alignment.

Settings: `--platform illumina`, `analysis_window: [53, 122]`, `min_align_len: 60` (the profile's 100 would discard nearly every read, since only ~70–97 bp of each read lies on the reference), the bundled J꜀4 contaminants, and 20,000 reads per mouse from the first 300,000 pairs of one run each. Mode: recommended (quality-calibrated); the four *Aicda*-/- mice are the controls, each tested against the other three.

## PCR duplicates

80–100 ng of genomic DNA is about 15,000 genomes, sequenced to 10–20 million read pairs. Each starting molecule is therefore read hundreds of times. The study collapsed duplicates, using the VDJ junction as a cell barcode, and kept clusters of at least 10 reads. nanoshm does not deduplicate.

Mutation frequency (mutations per base sequenced) and the share of mutated reads stay valid under duplication. Counts of distinct mutations, the single-read estimate and per-site support do not, because one molecule's mutation appears in many reads. Those are therefore not interpreted here. The single-read estimates in this run swing between −12 and +20 for exactly this reason.

## Result

| Group | Mice | Mutations per 10³ bp (mean) | Reads mutated (mean) |
|---|---:|---:|---:|
| Wild type | 7 | 1.81 (0.72–3.96) | 10.7% |
| LDHA-KO | 7 | 1.50 (0.00–2.21) | 9.0% |
| *Aicda*-/- | 4 | 0.10 (0.00–0.20) | 0.6% |

- Wild type against *Aicda*-/-: about 18-fold, Welch on log values P = 0.015 (P = 0.013 for reads mutated).
- LDHA-KO against *Aicda*-/-: P = 0.034.
- The *Aicda*-/- mice carry 0–4 calls each in 20,000 reads.
- G:C mutations in wild type fall at WRC/GYW hotspots 47–71% of the time.

The quality-calibrated model, chosen on nanopore J꜀4, separates AID-active from AID-null mice on Illumina reads from another laboratory, another amplicon design and another part of the locus.

**LDHA-KO mouse KO4 has no calls, and that is correct.** Only 2.3% of its reads carry any Q30 mismatch in the window, against 16% for LDHA-KO3 and 27% for wild-type 2. It looks like an AID-null library, so the cause is probably the sample (few germinal-centre B cells, for instance), not the caller.

## Speed

Merging 300,000 read pairs took 95 of the 140 seconds per mouse. Comparing the overlap with an early exit and looking up quality characters from a table brought it to about 38 seconds, with identical output.

## Limits

- **A 70-bp window.** It tests the error model on a short, heavily targeted stretch. It does not repeat the intron window.
- **No deduplication.** A UMI- or junction-based collapse would make distinct-mutation counts and the single-read estimate usable on data like this.
- **One run per mouse, first 300,000 pairs.** The second run of each mouse is downloaded but unused.

## Reproduce

`/Users/abdulmateenaderinto/Nanoshm/public_data/run_prjna973112.py 20000`, after `subsample.sh` for each run. Both scripts sit with the data, outside the repository.
