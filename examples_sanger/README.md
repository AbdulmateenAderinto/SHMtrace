# Sanger example: J<sub>H</sub>4 clones from Maea-deficient B cells

290 Sanger traces of J<sub>H</sub>4 intron clones from 12 mice, the cohort behind the
J<sub>H</sub>4 mutation analysis of

> Barbulescu P, Chana CK, Wong MK, et al. FAM72A degrades UNG2 through the GID/CTLH
> complex to promote mutagenic repair during antibody maturation.
> *Nature Communications* (2024).

| Group | Mice | Clones |
|---|---|---|
| Wild type (`WT`) | 2w1–2w4 | 22, 28, 23, 21 |
| *Maea*<sup>fl/fl</sup> *Mb1*<sup>cre</sup> (`Maea_KO`) | 2k1–2k4 | 58, 30, 50, 42 |
| *Aicda*<sup>−/−</sup> (`Aicda_KO`, the controls) | 2a1–2a4 | 4, 2, 9, 4 |

Shared with the authors' agreement. In the app, **Try it with Sanger sample data** on the Set up tab
loads it; from Python:

```python
from shmtrace import project, samplesheet
sh = samplesheet.load("examples_sanger/samples.csv")
res = project.run_project(sh.project("src/shmtrace/references/mouse_jh4_intron.fasta", "jh4_intron",
                                     contaminants=["src/shmtrace/references/mouse_jh4_contaminant.fasta"]))
```

## What is in the folder

- `traces/berrumar_<order>_chromats/`: the traces as delivered, one folder per sequencing order.
- `traces/filed_twice/`: four traces the original analysis filed under two mice (orders 391226
  and 393491). Each is a byte-for-byte copy of a trace in its order's folder; `samples.csv`
  marks both rows `shared_trace`, so each is analysed under both mice, as it was filed.
- `samples.csv`: one row per clone: file, mouse (`sample`), group, control, order, well, the
  label on the trace and a note naming the analysis folder the clone was filed in.

## How the clones were assigned to mice

From the plates themselves: the well numbers in order 390697 fall in 8-well blocks, one block
per mouse, and the mouse of each clone follows from its plate and well, never from its
mutations. Folder contents were checked against the original analysis project.

## What to expect

With SHMtrace's defaults (JH4 preset, Sanger platform): 290 clones analysed, three set aside
by the clone checks (one contaminant, two that do not map). Mutation frequency per 10³ bp:

| Mouse | WT 2w1 | 2w2 | 2w3 | 2w4 | Maea 2k1 | 2k2 | 2k3 | 2k4 |
|---|---|---|---|---|---|---|---|---|
| per 10³ bp | 20.4 | 17.0 | 11.8 | 10.0 | 3.5 | 2.3 | 4.5 | 2.4 |

Wild type is 4.66-fold above *Maea*-deficient by geometric mean (4.7-fold in the paper), and
both are well above *Aicda*<sup>−/−</sup> (0–0.5 per 10³ bp). The clone checks also flag
clones that are identical across mice, the signature of carry-over between plates; they are
listed in the QC table and left for the person analysing to decide.
