# Validation on independent public data, and the GHVPDN run

*11 October 2026.*

The calibration in [`ERROR_MODELS.md`](ERROR_MODELS.md) was measured on the
paper's own libraries. Everything below is from data SHMtrace had never seen: a
public nanopore SHM deposit from another laboratory, and the AID-variant /
FAM72A run GHVPDN.

## 1. BioProject PRJNA1481264 (independent laboratory)

ONT PromethION amplicons from the Institut de recherches cliniques de Montréal,
released 18 August 2026, 115 runs. It is the only public ONT amplicon SHM set
with matched *Aicda*-/- controls on the loci SHMtrace targets. There are no JH4
*Aicda*-/- runs, so the Sμ switch region and the AID off-target *Pim1* were used
instead. Each run was subsampled to 6,000 reads and analysed at 2,000.

There is no published reference for these amplicons, so for each locus the
reference is the majority consensus of one *Aicda*-/- library, iterated to
convergence. A consensus of knockout reads carries no somatic mutations by
construction, which is what a reference has to be here.

**Signal against background** (mean unique-mutation density ×10⁻⁵, each control
called against the other controls only):

| Locus | Mode | AID-competent | *Aicda*-/- | Separation |
|---|---|---:|---:|---:|
| Sμ, LPS-activated B cells | geneious | 7.73 | 0.28 | 28-fold |
| | v0.1 | 4.62 | 0.08 | 56-fold |
| | **recommended** | **9.25** | **0.06** | **167-fold** |
| | strand | 3.60 | 0.11 | 33-fold |
| Sμ, CH12 | geneious | 8.99 | 0.44 | 20-fold |
| | v0.1 | 4.32 | 0.11 | 39-fold |
| | **recommended** | **14.57** | **0.00** | **no background** |
| | strand | 3.75 | 0.00 | no background |
| *Pim1* (AID off-target) | geneious | 0.33 | 0.04 | 8.5-fold |
| | v0.1 | 0.02 | 0.04 | none |
| | **recommended** | **0.21** | **0.02** | **11-fold** |

The recommended mode finds more real mutations than Geneious's rule on both Sμ
sets and calls less background, so the separation is five to eight times better.
At *Pim1*, where AID acts weakly and the true signal is near the noise, it is the
only mode that still separates the genotypes; the v0.1 default finds nothing.

**Spike-ins into held-out controls** (60 planted sites × 3 draws, planted into
one *Aicda*-/- library with the null fitted from the other two):

| Locus | Mode | Recall | False calls per sample | FDR |
|---|---|---:|---:|---:|
| Sμ | geneious | 0.53 | 5.0 | 0.14 |
| | v0.1 | 0.21 | 1.7 | 0.12 |
| | **recommended** | **0.57** | **1.7** | **0.05** |
| | strand | 0.21 | 0.7 | 0.05 |
| *Pim1* | geneious | 0.52 | 1.0 | 0.03 |
| | v0.1 | 0.13 | 0.3 | 0.04 |
| | **recommended** | **0.57** | **0.3** | **0.01** |
| | strand | 0.11 | 0.0 | 0.00 |

The recommended mode dominates on both: more recall than Geneious's rule with a
third of the false calls. Recall is lower than on JH4 (0.57 against 0.70)
because these amplicons are shorter and were analysed at 2,000 reads, so a
planted mutation in one read is a smaller fraction of the column.

**What this shows.** The default was calibrated on the paper's own flow cells,
a different laboratory, different libraries and a different basecaller version.
It transfers. That is the main thing an external set can tell us, and it is the
argument for shipping this default rather than one tuned per run.

## 2. Run GHVPDN: AID variants and FAM72A

20 N114 libraries, five in each group. Confirmed by Abdulmateen on 11 October 2026:

| Group | What it is |
|---|---|
| AF | AID + FAM72A |
| AUF | AID up-mutant (Neuberger laboratory) + FAM72A |
| ADF | AID + G1-degron FAM72A (FAM72A degraded outside G1, matching its cell-cycle-restricted expression) |
| AUDF | AID up-mutant + G1-degron FAM72A |

The run has no empty-vector or no-AID library of its own, so the background
comes from the four 75K3B3 empty-vector N114 libraries, a different run.
**A cross-run background is not a valid null**: nanopore error is run-specific and
the basecaller may differ between orders. Treat the table below as a ranking,
not as calibrated mutation frequencies.

Mean unique-mutation density ×10⁻⁵ at 1,000 reads per library (500 reads in
brackets, which agrees):

| Group | Recommended | Geneious mode | Per library (recommended, 1,000 reads) |
|---|---:|---:|---|
| AF | 15.9 (26.1) | 15.7 | 0.0, 22.9, 20.5, 11.1, 25.2 |
| ADF | 16.2 (26.7) | 16.1 | 27.7, 5.7, 13.4, 20.7, 13.7 |
| AUDF | 6.5 (9.9) | 9.1 | 0.0, 13.6, 6.4, 9.8, 2.6 |
| AUF | 5.3 (9.7) | 7.0 | 4.7, 7.7, 4.2, 4.2, 5.4 |
| 75K3B3 AID (reference point) | 16.5 | 17.4 | 14.8, 12.4, 11.2, 12.8, 31.5 |
| 75K3B3 empty vector | 0.0 | 0.7 | 0.0 ×4 |

Three things to take from it:

1. **The up-mutant groups are lower, not higher.** AUF and AUDF come out at about
   a third of AF and ADF in both modes. An AID up-mutant is expected to mutate
   more, so either the construct is not behaving as an up-mutant in these cells,
   the libraries differ in transfection efficiency or input, or something about
   the amplicon differs between groups. This is the result that most needs
   checking at the bench.
2. **AF_1 and AUDF_1 are at background** (0.0 in recommended mode; Geneious mode
   agrees for AF_1 at 1.0), and AUDF_6 is close to it. These look like libraries
   where the construct did not work. Decide whether to exclude them before
   averaging, because with n = 5 one dead library moves the group mean a long way.
3. **Within-group spread is large** (AF ranges from 0 to 25). With this spread and
   n = 5, the AF-versus-AUF difference is the only one worth claiming, and only
   once the dead libraries are resolved.

The Geneious exports for this run cover only about 1,000 to 1,500 of the 3,000
to 5,000 reads per library, and the read set Geneious used is not recorded, so
a site-level comparison against those exports is not meaningful here. The
densities above are from SHMtrace reading the FASTQs directly.

**To make this run publishable**, sequence an empty-vector or catalytically dead
AID library on the same flow cell as the samples. That single addition turns the
ranking above into calibrated frequencies.

## Reproducing

```bash
# public data: one run, then the whole set
fasterq-dump SRR39285857            # Sμ, LPS, Aicda-/-
shmtrace-benchmark --sheet sheet_su_lps.csv --reference Igh_Su_consensus.fasta \
  --preset preset_su_lps.yaml --modes geneious,v0.1,recommended,strand --out pub_su_lps
shmtrace-validate --control Su_LPS_AIDKO_r1.fastq.gz \
  --null Su_LPS_AIDKO_r2.fastq.gz --null Su_LPS_AIDKO_r3.fastq.gz \
  --reference Igh_Su_consensus.fasta --preset preset_su_lps.yaml \
  --profiles geneious,v0.1,recommended,strand --depth 2000 --out pubval_su_lps

# GHVPDN
shmtrace-benchmark --sheet ghvpdn_sheet.csv --reference N114_geneious_655.fasta \
  --preset n114_vregion --modes geneious,recommended --depth 1000 --out ghvpdn1000
```

The sheets, the derived references and the per-locus presets are in
`shmtrace_data/public_PRJNA1481264/`.
