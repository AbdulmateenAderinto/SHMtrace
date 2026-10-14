# Changes in Ver22 (14 October 2026)

Ver22 reports the analyses as rerun on 13 October 2026 with SHMtrace v0.2's recommended mode, the
same analysis as the Cell Reports Methods article:

1. **False-discovery rate estimated from the controls.** Ver21 applied Benjamini–Hochberg at 0.05.
   Ver22 keeps the 0.05 target but estimates it with the held-out controls as decoys (target–decoy;
   Elias and Gygi, 2007, added to the references). The leave-one-control-out margin is now 5 for JH4
   and 10 for N114 and the full-length run.
2. **N114 empty-vector subclone C3 is its second library** (75K3B3_2_C3_Empty_March_2024_2). The
   first library was noisier than the other controls (Cell Reports Methods article); the two shared
   no calls.

| Result | Ver21 | Ver22 |
|---|---|---|
| JH4 Wt calls | 205 / 340 / 73 | 221 / 344 / 73 |
| JH4 Aicda−/− calls | 0 / 2 / 1 | 0 / 2 / 0 |
| JH4 fold (ONT) | ≥ 137×, p = 0.0025 | ≥ 141×, p = 0.0027 |
| Background | 0.48 × 10−5 (Q53), 3 calls | 0.48 × 10−5 (Q53), 2 calls |
| N114 pMX-Aicda / empty | 75 / 59 / 166 vs 0 / 0 / 0 | 70 / 49 / 173 vs 1 / 0 / 0 |
| N114 fold | ≥ 92×, p = 0.0048 | ≥ 85×, p = 0.0071 |
| Sanger–ONT spectrum | r = 0.93, slope 1.00 | r = 0.92, slope 0.97 |
| Depth: share of maximum at 500 reads | 94% (JH4), 91% (N114) | 93% (JH4), 74% (N114) |
| Full-length Wt vs Aicda−/− | 39 vs 0.51 (~78×) | 40 vs 0.18 (~229×) |
| Usable reads per library | mean 1,620 | mean 1,671 (C3's second library) |

Every changed number is computed from the outputs by `revise_manuscript_v22.py`; Ver21's highlights
are cleared, so the yellow text in Ver22 is what changed. The text now states that the empty-vector
controls had a single call between them, and Supplementary Figure S7's title no longer says that
discovery approached the maximum by 500 reads.

Figures 1, 2, 3, 4, S1, S3–S7 and S11 and the graphical abstract were regenerated from the rerun
(`scripts/`, `figure1_scicanvas/`); Figure 1 and S11 describe the control-based FDR.
