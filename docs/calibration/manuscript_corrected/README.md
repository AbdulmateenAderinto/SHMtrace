# Manuscript figures, regenerated with the corrected Sanger data

The manuscript's own figure scripts (Ver13, `6_Code/figure_scripts`), copied here and
changed only where the Sanger data enter. Nanopore inputs are read, unchanged, from the
original analysis folder; nothing there is modified.

| Figure | What changed |
|---|---|
| `Figure2_JH4_combined` | Panel B: corrected Sanger densities (unique mutations over the bases each clone read). Panels A and C unchanged. |
| `Figure3_concordance` | Corrected Sanger matrices and mutations per clone. Panels B and D show each platform's pooled proportion with a 95% Wilson interval (mice as points sized by their mutations; Sanger Wt2's 50% A→T point is drawn at the axis top). The per-type and Ts/Tv P values are computed (paired t, Šídák), not typed in. |
| `FigureS5_individual_Sanger_matrices` | The manuscript's Supplementary Figure S2, from the corrected Sanger calls. |
| `FigureS6`, `FigureS7` | Unchanged inputs (nanopore J꜀4 and N114), regenerated for completeness. |

Corrected Sanger calls: `sanger_calls/` (one Geneious-shaped CSV per mouse, positions on the
529-bp reference, and `corrected_sanger.json`), built by nanoshm from the organised .ab1
traces (`~/Nanoshm/sanger_ab1_by_mouse`; see `../SANGER_VS_NANOPORE.md`).

Rebuild: `cd scripts && python figure2_jh4_combined.py && python figure3_concordance.py && python figureS5_S7_matrices_redesign.py`
