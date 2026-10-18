# Two public datasets added to the validation, 17 October 2026

Branch `public-validation-2026-10`. This records what was run, what it establishes, what it does not, and what
was rejected. Numbers here come from the committed outputs; nothing is quoted from the source papers as if it
were our result.

## What was added

| | PacBio class-switch junctions | Illumina Sμ mutation libraries |
|---|---|---|
| Accession | PRJNA831666 | GSE77788 / PRJNA311552 |
| Source | Vincendeau et al. 2022, Nat Commun 13:3707 | Compagno et al. 2017, Nature 542:489 |
| Design | WT, *Shld1*^−/−, *Xlf*^−/−, double, 3 mice each | WT and *Aicda*^−/−, DMSO / idelalisib / duvelisib, 3 each |
| Downloaded | 12 runs, 58 MB | 18 runs, 30 files, 540 MB |
| Tests | the junction module (`shmtrace.switch`) | the control-calibrated substitution caller |
| Does **not** test | the substitution caller or its FDR calibration (no AID-null arm) | the drug effect (see below) |
| Scripts | `crm/build_switch_reference.py`, `crm/switch_pacbio.py` | `crm/public_smu_pi3k.py` |
| Outputs | `crm/out/switch_pacbio/` | `crm/out/public_smu/` |
| Figure | S16 (`crm/figs/figS16_switch_pacbio.py`) | S17 (`crm/figs/figS17_public_smu.py`) |

Raw reads live outside the repository, in `~/Nanoshm/public_data/PRJNA831666/` and
`~/Nanoshm/public_data/GSE77788/`, as the other public datasets do. They are not redistributed.

## What the new evidence establishes

**The junction module recovers a repair phenotype it was never tuned for.** Mean microhomology was 1.01 nt in
wild type, unchanged by loss of SHLD1 (1.02 nt, p = 0.92), and rose to 1.69 nt without XLF (p = 0.0092) and
1.59 nt in the double mutant (p = 9.3 × 10⁻⁵). Junctions with 4 nt or more of microhomology rose from 7% to 18%
and blunt junctions fell from 21% to 7% in the double mutant. Replicate ranges do not overlap between wild type
and either XLF-deficient genotype. This is the shift toward microhomology-mediated end joining expected when
classical non-homologous end joining is impaired, obtained with no parameter fitted to the dataset.

The source paper's own emphasis is different — it reports SHLD1 restricting resection and being required for
orientation-specific joining — and that analysis was not reproduced here. Recovered deletion boundaries are not
a measurement of pre-ligation resection, and no claim of that kind is made.

**The control calibration transfers to an independent laboratory, platform and perturbation.** Over the Sμ
repeats, using only the twelve libraries sequenced the same way, the nine *Aicda*^−/− libraries produced one
call between them while the three wild-type libraries produced 81, 162 and 172 (p = 0.0021). The wild-type
calls carry the AID signature: 5-mer targeting index 1.52–1.67 and 31–35% C→T or G→A.

## What it does not establish, and why

**The drug effect is not resolved, and we do not claim it.** In the window every library reads, wild-type means
rose from 18.7 unique mutations under DMSO to 36.3 under idelalisib and 54.7 under duvelisib — the direction the
original study reports — but replicates within a treatment differed by more than the treatments did (p = 0.64
and p = 0.38). Two properties of the deposit cause this:

1. **Layout is confounded with genotype.** All nine *Aicda*^−/− libraries are paired; six of the nine wild-type
   libraries are 300-bp single-end.
2. **The paired mates do not overlap.** R1 and the single reads cover reference positions 401–650 and 351–650 of
   the 651-bp amplicon; R2 covers 1–250, with a gap between. `io_utils.merge_pairs` keeps only 9–21% of pairs,
   and those are the shortest molecules — a length-biased subset that would have differed systematically
   between the layouts. Analysed as deposited, this alone could have manufactured a case/control difference.

The harmonisation — one read per library, over the 216 bases all eighteen read — removes the confound but lands
in the half of the amplicon where AID acts least. The repeats, where the signal is, are reached only by R2, and
only one wild-type replicate per treatment is paired. So the dataset supports a genotype contrast and not a
treatment contrast. That is a property of the deposit, not of the method.

**Junction counts are not switching efficiency.** Reads per unique junction varied from 1.1 to 23 between
samples, so junction yield reflects recovery of long-range PCR products. Only junction shape is compared.

## Code-change log

**Confirmed defect, not patched — `switch.py` ambiguity flag.** `summarise()` reports `ambiguous_pct`, which
reads **0.0% in all twelve samples**. That is not evidence of unambiguous breakpoints. The flag fires only when
minimap2 returns a second alignment of the acceptor part, and for these reads it returns exactly two alignments
per junction read and no secondaries, so the test can never fire. Measuring the same property directly — by
re-aligning each junction's first 60 acceptor bases on their own, with secondary alignments kept, and counting
placements within five matched bases of the best — puts **19–47%** of unique junctions at more than one equally
good acceptor position. `tests/test_switch.py` has five tests and none covers ambiguity.

This was left unpatched deliberately: it is core calling code used by the app and CLI, the manuscript never
quotes the built-in figure, and a behaviour change late in a revision would need regression coverage it does not
have. The direct measure lives in `crm/switch_pacbio.py` and is what Figure S16D and the legend report. **A fix
is proposed, not applied** — either compute ambiguity the direct way inside `summarise()`, or return `None`
rather than `False` when no secondary alignment was available, in both cases with tests added first.

**No change to `src/`.** Nothing in the package was modified. The module's `map-ont` default was overridden to
`map-hifi` through the existing `preset` argument, which needed no code change.

**New files only**, all of them under `crm/`: `build_switch_reference.py`, `switch_pacbio.py`,
`public_smu_pi3k.py`, `figs/figS16_switch_pacbio.py`, `figs/figS17_public_smu.py`, this file. Edits to existing
files are confined to `manuscript/` (one Results sentence, two STAR Methods sections, two figure legends, two
key-resources rows, the numbers block) and `manuscript/resolve_refs.py` (two PubMed queries).

## Two things the archives get wrong

**PRJNA831666 files are named `_subreads` but are high-accuracy reads.** Median identity against the bundled Sμ
fragment is 99.85% over alignments of 600 bp or more (10th percentile 99.4%). Raw PacBio subreads sit near
85–90%. `map-hifi` is therefore used. Had the file name been believed, the preset and every downstream
expectation would have been wrong.

**GSE77788 single-end files are not merged pairs.** They are uniform 300-bp single-end reads from the same
forward primer as R1, not pre-merged mates — checked directly, since merged reads would vary in length.

## Word budget

The main text stood at 6,994 of 7,000 words before this work, so there was room for six words. The substance is
therefore in the STAR Methods and the two figure legends, which do not count toward the limit, and the Results
carry one sentence citing both figures. The six words were recovered by reverting a clause added earlier the
same day to the Figure 5E paragraph (the per-mouse G:C transition difference), which the figure shows anyway.
**If the authors want these results stated in the Results, roughly 150 words must be cut elsewhere; that is an
authorial decision and was not made here.**

## Exact manuscript changes

Every change to the article, with its reason. Nothing else in `manuscript/` was touched.

| File | Change | Why |
|---|---|---|
| `crm_text.py` | **Added** to Results, after the mOrange paragraph: *"Two public datasets extended these tests: Sμ amplicons under PI3Kδ blockade[@compagno2017], and PacBio class-switch junctions[@vincendeau2022] (Figures S16 and S17)."* | Every supplementary figure must be cited in the main text. Held to 23 words by the word cap. |
| `crm_text.py` | **Reverted**, in the Figure 5 paragraph: *"…not equivalence, and nanopore gave a {pm_gcts} percentage-point lower share of G:C transitions in each mouse."* → *"…not equivalence."* | Added earlier the same day; Figure 5E shows it. Freed the words for the sentence above. `pm_gcts` is still computed and unused. |
| `crm_methods.py` | **Added** two STAR Methods sections, *Public PacBio class-switch junctions* and *Public Illumina Sμ mutation libraries* | Where the detail lives, since STAR Methods do not count toward the word limit. |
| `crm_methods.py` | **Added** two key-resources rows (PRJNA831666, GSE77788/PRJNA311552) | Deposited data used must be listed. |
| `crm_methods.py` | *"Document S1. Figures S1–S15…"* → *"…S1–S17…"* | Two figures added; the index must match. |
| `build_crm.py` | **Added** the Figure S16 and S17 legends | New figures. |
| `crm_numbers.py` | **Added** two blocks computing `pb_*` and `smu*` placeholders from the committed outputs | Every number in the text is computed, never typed. |
| `resolve_refs.py`, `refs.json` | **Added** PubMed queries `35764636[uid]` and `28199309[uid]` | Two new citations, fetched from PubMed rather than composed. |

## Datasets examined and not used

| Dataset | Why not |
|---|---|
| PRJNA1235282 (Lescale 2025, PacBio CH12 Sμ–Sα) | Same question as PRJNA831666 and the same research group; 3.83 GB and more sample bookkeeping for no independent validation. One PacBio CSR example was the brief's instruction. |
| PRJNA1304899 (Marton 2025, PacBio CSR) | As above, 6.42 GB; richer but a larger project than this side paper warrants. |
| Bassan 2021 (PacBio TCR, figshare) | Selected populations with one extraction per group, not replicated AID-on/off. Known high-frequency selected mutations would be removed by the ≥50% germline heuristic, needing a special-cased reference. |
| GSE102944 (Álvarez-Prado 2018) | One pooled AID-null library. Splitting it would manufacture controls that are not independent. Capture reads also need a different preparation path from full-length amplicons. |
| PRJNA1337487 (Wu 2026) | 907 runs, ~1 TB mixed assays; run titles did not identify a replicated AID-off mutation series. |
| Singh 2019 (mBio), DT40 AID study | Sanger sequences available on request, or deposits that do not establish the full positive/control panel. |

These were assessed from the metadata gathered in the dataset-search folder and from the archives; their reads
were not downloaded.

## Rerun

```bash
PY=~/Nanoshm/nanoshm_v0.2/.venv/bin/python
cd ~/Nanoshm/nanoshm_v0.2/docs/calibration

$PY crm/build_switch_reference.py      # fetches the IgH locus once, verifies it, writes the Sμ/Sγ1 reference
$PY -W ignore crm/switch_pacbio.py     # -> crm/out/switch_pacbio/
$PY -W ignore crm/public_smu_pi3k.py   # -> crm/out/public_smu/  (slow; run it with jobs=1 if memory is tight)

cd crm/figs && $PY figS16_switch_pacbio.py && $PY figS17_public_smu.py
cd ~/Nanoshm/nanoshm_v0.2/docs/calibration && $PY -W ignore crm/manuscript/build_crm.py crm/build
```

Accessions are listed in the two analysis scripts (`SAMPLES` and `RUNS`), each run tagged with its genotype,
treatment, replicate and, for GSE77788, its library layout. `crm/out/switch_ref/reference_provenance.json`
records the Ensembl request, the window coordinates and the two alignment checks that must pass before the
reference is written.
