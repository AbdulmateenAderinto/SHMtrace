# Changes in Ver19–Ver21

## Ver21 (13 October 2026): SHMtrace end to end on the delivered reads

Up to Ver20 the ONT results were SHMtrace calls made on 500-read subsets produced by an earlier
preprocessing script (15-mer competitive selection, a 900–2,000 bp JH4 / 400–900 bp N114 length
window, a 651-nt N114 selection reference, reservoir sampling with seed 42). Ver21 has SHMtrace do
all of the work from the FASTQ files as delivered by the sequencing provider:

1. raw length ≤ 2,000 bp (JH4) or ≤ 900 bp (N114), no minimum; mean expected error ≤ 0.01 per base
2. minimap2 (map-ont) against the target plus competing references (JH4: the JH4-like paralogue and a
   recurrent contaminant haplotype); reads aligning best elsewhere, chimeric reads and alignments
   < 800 bp (JH4) / < 400 bp (N114) removed
3. JH4 reads with ≥ 10 paralogue-specific bases removed
4. 500 usable reads drawn at random per library (seed 0); usable reads 962–2,905 (mean 1,620)
5. calling unchanged: control error model per substitution, strand and quality; leave-one-control-out
   factor 5; Benjamini–Hochberg FDR 0.05; sites ≥ 50% removed

| Result | Ver20 (subsets) | Ver21 (delivered reads) |
|---|---|---|
| JH4 Wt calls | 155 / 347 / 69 | 205 / 340 / 73 |
| JH4 Aicda−/− calls | 0 / 3 / 1 | 0 / 2 / 1 |
| JH4 fold (ONT) | ≥ 108×, p = 0.0018 | ≥ 137×, p = 0.0025 |
| Background | 0.55 × 10−5 | 0.48 × 10−5 (Q53) |
| N114 pMX-Aicda / empty | 69 / 68 / 150 vs 0 / 0 / 1 | 75 / 59 / 166 vs 0 / 0 / 0 |
| N114 fold | ≥ 89×, p = 0.0034 | ≥ 92×, p = 0.0048 |
| Sanger–ONT spectrum | r = 0.91 | r = 0.93 |
| Depth: share of maximum at 500 reads | 91% | 94% (JH4), 91% (N114) |
| 30-read JH4 controls | uncalibrated (12–29 calls) | calibrated (0–3 calls) |
| Wt3 tree clone | 84 reads, 87% fit | 69 reads, 70% fit |

Text: Methods 2.7 rewritten as SHMtrace's read preparation with per-step losses; 2.8, 2.9, Results
3.1–3.5, Discussion, Conclusions, Abstract, Tables 1, S1, S2, all legends and Data availability
updated. The depth section is rewritten: density now declines with depth as a saturating numerator
should. The earlier Geneious analysis is kept only as a labelled comparison in Section 3.1.

Figures: 2, 3, 4, S1, S3–S7, S9, S10 and the graphical abstract regenerated from the new runs.

Builder: `revise_manuscript_v21.py`; package: `build_submission_v21.py`.

## Ver20: Figure 1 and Supplementary Figure S11 in SciCanvas

Figure 1 redrawn as a high-level overview; the detailed workflow became Supplementary Figure S11.
Both meet Cell Press figure rules (174 mm, Arial 6–8 pt, strokes 0.5–1.5 pt, ≥ 1,000 dpi, no
red–green pairing) and pass SciCanvas's figure check. PCR is drawn as a thermal cycler (traced icon),
and nanopore sequencing as a pore. Builder: `revise_manuscript_v20.py`.

## Ver19: mutation trees

Supplementary Figure S10 (Wt3 clone mutation tree) with Methods and Results paragraphs.
Builder: `revise_manuscript_v19.py`.
