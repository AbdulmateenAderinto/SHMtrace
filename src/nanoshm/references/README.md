# Bundled reference sequences

Pick one of these in the app instead of supplying a FASTA, or supply your own.

| File | What it is | Where it came from |
|---|---|---|
| `mouse_jh4_intron.fasta` | Mouse IgH JH4 downstream intron, 1,230 bp | The Aderinto & Martin deposit (`Jh4 - Base`) |
| `mouse_jh4_contaminant.fasta` | The JH4 contaminant sequence to exclude | The same deposit |
| `n114_vregion.fasta` | N114 hybridoma heavy-chain V region, 655 bp | The same deposit (`seq_23353`) |
| `mouse_igh_smu_ch12.fasta` | Mouse IgH Sμ switch-region amplicon, 651 bp | Majority consensus of the Aicda-/- CH12 reads of BioProject PRJNA1481264 |

A consensus built from knockout reads is the sequence those reads agree on, so it
carries no somatic mutations, but it is not a curated germline entry. For
publication, check it against your own germline.
