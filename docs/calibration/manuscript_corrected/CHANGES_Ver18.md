# Changes in Ver18 (since Ver16)

*12 October 2026.*

Every number in Ver18 now comes from one pipeline:
- **Sanger:** SHMtrace in Sanger mode, run on the original `.ab1` chromatograms.
- **Nanopore:** SHMtrace in its recommended mode, run on the deposited 500-read subsets.

Geneious is no longer used for any result. Its exports stay in the data deposit, and Section 3.1 compares them with the new calls.

## How the nanopore calls are made

Each candidate substitution is tested against the error rate measured in the negative-control libraries:
- the controls are *Aicda*−/− mice for J~H~4 and the empty vector for N114;
- the error rate is specific to the substitution, the strand and the base quality.

The control rates are inflated by a factor chosen by leave-one-control-out; the factor was 5 for both assays. Benjamini–Hochberg correction is then applied at FDR 0.05, and each control is called against the other two controls.

Geneious instead applied an unadjusted p ≤ 0.05 per candidate, with no reference to the controls.

## What changed in the results

| | Ver16 (Geneious ONT; Sanger from GenBank records) | Ver18 |
|---|---|---|
| J~H~4 ONT calls, Wt mice | 161 / 342 / 68 | 155 / 347 / 69 |
| J~H~4 ONT calls, *Aicda*−/− mice | 2 / 6 / 3 | 0 / 3 / 1 |
| J~H~4 ONT, Wt vs *Aicda*−/− | ~47-fold, p = 0.0037 | ≥ 108-fold, p = 0.0018 |
| ONT background | ~1.2 × 10⁻⁵ (Q49) | ~0.55 × 10⁻⁵ (Q53) |
| Sanger clones (Wt) | 36–38 per mouse | 26–40 per mouse (re-assigned; see Ver17) |
| Sanger, Wt vs *Aicda*−/− | ~41-fold, p = 0.0045 | ~38-fold, p = 0.0009 |
| Spectrum r (Sanger vs ONT) | 0.97 | 0.91 [0.69–0.97], slope 1.02 |
| Platform × type interaction | p = 0.875 | p = 0.50; smallest adjusted p 0.46 |
| N114 calls, pMX-*Aicda* | 80 / 71 / 147 | 69 / 68 / 150 |
| N114, empty vector | density 1.69 × 10⁻⁵ | 0 / 0 / 1 calls |
| N114 fold | ~18-fold, p = 0.0147 | ≥ 89-fold, p = 0.0034 |
| Full-length J~H~4 spectrum vs 529 bp | r = 0.989, n = 739 | r = 0.974, n = 680 |
| Calls recovered at 500 reads | 97% (J~H~4), 95% (N114) | 91% for both |

The Wt mice are almost unchanged. The difference is in the controls: calling against them removes most of the background. The J~H~4 and N114 group differences therefore become larger. Because one *Aicda*−/− mouse and two empty-vector subclones have no calls, those folds are lower bounds: each such library is plotted and tested at its one-event density.

## What changed in the depth analysis (Supplementary S6, S7)

The depth curves are no longer a simple saturating curve, for two reasons:

- **Under-calling below ~250 reads.** A mutation seen in only one read rarely passes FDR 0.05 against the controls. So a mouse whose mutations sit mostly in single reads (Wt2) is under-called: it had 3 and 2 calls at 60 and 125 reads, against 110 at 250. Clonal mutations (Wt3) are called at every depth.
- **No calibration at 30 reads.** The background model cannot be calibrated, so the J~H~4 controls themselves get 12–29 calls each. That point is shaded and described as not interpretable.

The Results text (3.5) and the S6/S7 legends now say this. The sentence claiming that density declines with depth "as a mathematical property" was removed, because with this caller density rises with depth for N114.

## Text and figures changed

- **Methods 2.8** is one SHMtrace section covering both assays. The separate N114 paragraph was merged into it, and the Geneious description was replaced, with Geneious cited for the original calls.
- **Text updated with the new results** (Results 3.1–3.5, the Abstract, the Discussion and the Conclusions):
  - every number;
  - the background paragraph, with the Geneious comparison;
  - "no multiplicity correction" removed: the new calls are FDR-controlled.
- **Figures replaced:** Figures 1–4 and S1, S3–S7 and S9. In Figure 1, step 5 now reads "vs. controls, FDR ≤ 0.05". The graphical abstract was regenerated.
- **Figure 2:**
  - The censored *Aicda*−/− mouse is shown as an open circle.
  - In panel C, the shaded Sanger window is corrected to positions 101–629; it was drawn to 540.
- **Tables:**
  - S1: ONT rows now give reads analysed, candidates tested, FDR-passing calls and density, and the column headers are renamed to match.
  - S2: all ONT rows updated.
  - Table 1: the background figure updated.

## Still to check by hand

- **Read-length window.** Figure 1B says 800–2,000 bp for J~H~4, but Methods 2.7 and the Figure 1 legend say 900–2,000 bp. This is unchanged since Ver13; confirm which is right.
- **Highlighting.** Yellow highlighting marks every changed passage; remove it before submission.
- **Gene-name formatting.** It was restored automatically in the changed passages; skim them once in Word.
