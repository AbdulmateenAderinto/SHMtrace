# Manuscript figures, regenerated with the corrected Sanger data

The manuscript's own figure scripts (Ver13, `6_Code/figure_scripts`), copied here and
changed only where the Sanger data enter. Nanopore inputs are read, unchanged, from the
original analysis folder; nothing there is modified.

| Figure | What changed |
|---|---|
| `Figure2_JH4_combined` | Panel B: corrected Sanger densities (unique mutations over the bases each clone read). Panels A and C unchanged. |
| `Figure3_concordance` | Corrected Sanger matrices and mutations per clone. Panels B and D show each platform's pooled proportion with a 95% Wilson interval (mice as points sized by their mutations). The per-type and Ts/Tv P values are computed (paired t, Šídák), not typed in. |
| `FigureS5_individual_Sanger_matrices` | The manuscript's Supplementary Figure S2, from the corrected Sanger calls. |
| `FigureS6`, `FigureS7` | Unchanged inputs (nanopore J꜀4 and N114), regenerated for completeness. |

Corrected Sanger calls: `sanger_calls/` (one Geneious-shaped CSV per mouse, positions on the
529-bp reference, and `corrected_sanger.json`), built by `export_sanger_calls.py` from
the .ab1 traces organised with `../organize_jh4_ab1.py --original` (`~/Nanoshm/sanger_ab1_original`):
every clone stays with the mouse it was filed under in the GenBank folders, as the previous student
assigned it. This is the author's decision: the manuscript measures nanopore performance and
concordance against the Sanger data as collected, not a re-curation of it. The evidence-based
reassignment (`~/Nanoshm/sanger_ab1_by_mouse`) is documented in `../SANGER_VS_NANOPORE.md` only.

| Mouse | Clones | Calls | Density (×10⁻⁵) |
|---|---|---|---|
| Wt2 / Wt3 / Wt4 | 41 / 38 / 27 | 177 / 230 / 111 | 840 / 1183 / 833 |
| Aicda⁻/⁻ 1 / 2 / 3 | 11 / 10 / 10 | 1 / 1 / 2 | 17.3 / 19.0 / 37.8 |

Wt vs Aicda⁻/⁻ ~41-fold (Welch on log10, P = 0.0012); spectrum Pearson r = 0.96, slope 1.06;
interaction P = 0.76; smallest Šídák P = 0.46; Ts/Tv adjusted P = 0.67 / 1.00.

`revise_manuscript_sanger.py` builds Ver17 from Ver16 (Ver16 untouched; edits highlighted
yellow): new numbers and figures, the TOST analysis removed (Methods, Results, abbreviation,
Table S2 row and legend), and two new supplementary figures — S8, the 529-bp window
(`figS_sanger_window_original`), and S9, agreement beyond the spectrum
(`geneious_original/fig3_agreement_beyond_spectrum`). The inclusion/curation figure is not used.

Rebuild: `python ../organize_jh4_ab1.py <ab1 dir> ~/Nanoshm/sanger_ab1_original <genbank dir> --original && python export_sanger_calls.py && cd scripts && python figure2_jh4_combined.py && python figure3_concordance.py && python figureS5_S7_matrices_redesign.py`
