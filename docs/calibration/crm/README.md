# Cell Reports Methods article and STAR Protocols protocol

The expanded SHMtrace paper (method, comparison with related callers, validation on several datasets)
and its companion protocol. Every number in the text comes from a file in `out/`, every figure from a
script in `figs/` or `scicanvas/`.

## Analyses (`out/`)

| Script | What it does | Output |
|---|---|---|
| `benchmark_real.py jh4 n114 su_lps su_ch12 pim1 pim1_deep` | SHMtrace, the Geneious rule and LoFreq (three settings) on the same analysed reads of real libraries with AID-null controls | `benchmark_real_*.json` |
| `spikein.py jh4 n114 su_lps pim1` | planted mutations in held-out controls, every caller | `spikein_*.json` |
| `lofreq_runner.py` | LoFreq 2.1.5 on exactly SHMtrace's reads (conda env `lofreq`) | – |
| `datasets.py` | the datasets: delivered J_H4/N114 reads, PRJNA1481264 subsets | – |
| `ircm_prepare.py` | PRJNA1481264: 6,000-read subsets, consensus references from AID-null libraries | `~/Nanoshm/public_data/PRJNA1481264/` |
| `shmtrace-simulate-lineages` (see LINEAGE_SIM.md) | simulated lineages on the delivered J_H4 controls | `lineage_sim_jh4/` |
| `../quant_accuracy.py --amplicons N114 --out crm/out/quant_accuracy_n114` | quantification accuracy, N114 on delivered reads | `quant_accuracy_n114/` |
| `single_read_spike.py` | single-read estimate against planted one-read mutations | `single_read_spike.json` |
| `public_illumina.py ldha ch12 hmces` | the three public Illumina datasets with the current SHMtrace | `public/*/` |
| `targeting_jh4.py`, `targeting_ch12.py` | 5-mer targeting index and dose-response | `targeting_*.json` |
| `jh4_extras.py` | J_H4/N114 projects with the targeting model, clones, single-read estimate | `extras_*/` |
| `make_tables.py` | Table S2 | `Table_S2_benchmark.xlsx` |

Inputs outside the repository: the deposit's delivered FASTQs, `~/Nanoshm/public_data/` (public sets, not
redistributed; Karst et al. data are GPL-3 and not committed), `~/Nanoshm/models/HH_S5F.csv` (exported from
SHazaM with `scripts/export_shazam_model.R`; AGPL-3, not committed).

## Figures

- `scicanvas/crm_figure1.js`, `crm_ga.js`, `star_figure1.js`: built with
  `manuscript_corrected/figure1_scicanvas/build.js` (graphical abstract with `GA=1`); Cell Press checks and
  SciCanvas's figure check run on every build.
- `figs/fig2_calibration.py` … `fig6_platforms.py`, `figS5_estimates.py`, `star_figs.py`: matplotlib, Arial,
  174 mm, 600-dpi PNG/TIFF and vector PDF.

## Documents

```bash
python manuscript/resolve_refs.py                      # references from PubMed -> manuscript/refs.json
python manuscript/build_crm.py build                    # article, highlights/eTOC/motivation, Document S1
python manuscript/build_star.py build                   # protocol
python manuscript/build_letters.py build                # cover letters
```

Bracketed items in the text ([to confirm]) are laboratory details only the authors can supply.
