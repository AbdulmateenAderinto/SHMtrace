# Bundled reference sequences

Pick one of these in the app instead of supplying a FASTA, or supply your own.

| File | What it is | Where it came from |
|---|---|---|
| `mouse_jh4_intron.fasta` | Mouse IgH JH4 downstream intron, 1,230 bp | The Aderinto & Martin deposit (`Jh4 - Base`) |
| `mouse_jh4_contaminant.fasta` | Sequences to exclude from JH4 analyses: the deposit's contaminant, and a JH4-like paralogue the JH4 primers also amplify (614 bp, about 89% identical to the JH4 intron) | The deposit (`Contaminant_JH4`); the paralogue is the majority consensus of 46 Sanger clones from the lab's JH4 clone set, October 2026 (`docs/calibration/SANGER.md`) |
| `n114_vregion.fasta` | N114 hybridoma heavy-chain V region, 655 bp | The same deposit (`seq_23353`) |
| `n114_vregion_contaminant.fasta` | A second V-region template the N114 F1/R2 primer pair also amplifies (476 bp, 4.0% divergent from the N114 V region over positions 180-655) | Majority consensus of 37 Sanger clones from the lab's N114 clone set, orders of February-March 2024; the F1/R2 pair was replaced by F1/R1 from April 2024 |
| `mouse_igh_smu_ch12.fasta` | Mouse IgH Sμ switch-region amplicon, 651 bp | Majority consensus of the Aicda-/- CH12 reads of BioProject PRJNA1481264 |

A consensus built from knockout reads is the sequence those reads agree on, so it
carries no somatic mutations, but it is not a curated germline entry. For
publication, check it against your own germline.
