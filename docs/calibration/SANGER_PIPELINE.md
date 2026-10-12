# The Sanger pipeline: sample sheet, QC flags, concordance, validation

*12 October 2026.*

The errors that changed the lab's JH4 Sanger result were not calling errors. They were bookkeeping errors:
- clones filed under the wrong mouse;
- trace copies renamed into another mouse's folder;
- contaminating clones under knockout labels.

The pipeline therefore makes the sample assignment explicit, checks it, and flags anything suspicious for a person to decide. **Nothing is corrected automatically.**

## 1. The sample sheet (`nanoshm.samplesheet`)

A CSV or tab-separated file with one row per trace:

```
file,sample,group,control,order,well,label,exclude,aliases
traces/434141_A7.ab1,WT2,WT,no,434141,A7,w2-1,,w2;wt2
traces/432625_E6.ab1,Aicda_A1,KO,yes,432625,E6,a1-1,AID-type mutations in a control,a1
```

| Column | Required | Meaning |
|---|---|---|
| `file` | yes | A trace, a folder of traces, or a `*` pattern. Paths are relative to the sheet, or absolute, or start with `~`. |
| `sample` | yes | The animal or independent source the trace belongs to. |
| `group` | | The sample's group. Defaults to the sample name. |
| `control` | | `yes` or `no`: whether the sample is a mutation-free control. |
| `clone_id` | | The clone's name in every table. Defaults to `order-well`, or to the file name. |
| `order`, `well`, `label` | | Where the trace came from. Kept with every flag. |
| `exclude` | | Any text here leaves the row out, and the text is recorded as the reason. |
| `aliases` | | How this sample is written on tube labels (e.g. `w3;wt3`). Enables the label check. |

`load()` reports every problem at once, each with its line number, and the sheet cannot be used until all are fixed:
- a missing file, or a file that isn't a read file;
- the same file listed twice;
- the same clone name used twice;
- a blank sample;
- an unreadable `control` value;
- a sample given two different groups or control statuses.

It warns about, but accepts:
- a trace that is byte-for-byte the same file as another;
- two different traces from the same order and well (usually a re-run).

```bash
nanoshm-sheet draft traces_folder/ -o samplesheet.csv   # a first draft from sub-folders; check every row
nanoshm-sheet check samplesheet.csv
nanoshm-sheet run samplesheet.csv --reference ref.fasta --region jh4_intron \
      --contaminant contaminant.fasta -o results/
```

In the app, type the sheet's path under **3 · Samples** and click **Load sheet**. The QC flags and the excluded clones appear under **Tables**, in the HTML report, and as `qc_flags.csv` and `excluded.csv`.

## 2. QC flags (`nanoshm.qc`, run on every Sanger project)

| Check | Severity | Fires when |
|---|---|---|
| `duplicate_file` | problem | Two clones are the same file, byte for byte. |
| `duplicate_clone` | problem | Clones from two animals agree at every base both read (at least 150 bp) and share 3 or more mutations. Two animals can't make the same mutated clone. |
| `shared_mutations` | problem / check | A clone shares more mutations with another animal's clone than independent animals share by chance. The test is Poisson: each mutation is weighted by how many animals carry it, counting each animal once, and the result is Bonferroni-corrected over the pairs. |
| `control_mutations` | problem | A control clone carries more mutations than its own base qualities predict, calibrated on the other control clones. |
| `label_mismatch` | check | The tube label names a different sample (needs `aliases`). |
| `germline_or_lineage` | check | The ≥ 50% germline rule removed differences that no other animal shows. This can mean an expanded lineage, or a hotspot mutated in several lineages, has lost its mutations. |
| `load_outlier` | check | A clone carries far more mutations than the rest of its animal: more than the upper quartile plus 3 × IQR, and at least 6. |
| `short_clone` | info | A clone reads less than a third of the window. |

Every clone that was not analysed is listed with the reason: excluded by the sheet, no quality segment, unmapped, contaminant or paralogue, alignment too short, no clean stretch, or minor haplotype.

## 3. Concordance between two measurements (`nanoshm.concordance`)

This is the manuscript's Sanger-vs-ONT comparison as a tested function, for any two runs of the same animals:
- **Spectrum:** pooled proportions with 95% Wilson intervals; Pearson r with a Fisher interval, and slope; the platform × type repeated-measures interaction; paired t tests with Šídák adjustment per type and for G:C Ts/Tv.
- **Position:** the profile along the window.
- **Hotspots:** WRC/GYW and WA/TW shares.
- **Density:** per-animal densities.

```bash
nanoshm-concordance sanger/project.json ont/project.json --pair WT2=WT2 --pair WT3=WT3 \
      --pair WT4=WT4 --labels Sanger,ONT
```

`docs/calibration/manuscript_corrected/manuscript_stats.py` now uses this module and gives the manuscript's numbers unchanged.

## 4. Validation with known answers

### Simulation (`nanoshm.sanger_sim`)

Each simulated project has:
- **Animals:** three wild-type and three knockout animals, about 30 clones each.
- **Lineages:** four to eight per wild-type animal, with AID-weighted mutations in and around the window.
- **Germline:** strain differences shared by every animal, and occasionally one difference private to an animal.
- **Realistic `.ab1` qualities:** a noisy start, a fading end, and errors drawn from the qualities.
- **Unusable traces:** failed reads, paralogue clones and short reads.
- **Planted bookkeeping errors**, one of each:
  - a misfiled clone with its label intact;
  - a misfiled clone whose label was also changed;
  - a re-typed copy (FASTA) of a clone, filed under another animal;
  - a copied trace file;
  - a mutated clone filed under a control;
  - a whole lineage filed under another animal.

The pipeline then runs exactly as on real data, from the sample sheet. Results from 40 projects (`nanoshm-validate-sanger --reps 40 --seed 101`; full output in `sanger_validation_summary.json`):

| | Result |
|---|---|
| Recall of true mutations (per animal) | mean 99.4% (lowest 83%, a knockout animal with 6 true mutations) |
| Precision in mutated animals (≥ 20 true mutations) | mean 95% (lowest 85%) |
| False calls | 21 per 10⁵ bases read: the single-pass Sanger floor at Q20 (your real Aicda⁻/⁻ clones: 17–38 per 10⁵) |
| Germline differences called as mutations | 0 |
| Failed or paralogue traces analysed | 0 |
| Mutated clone filed under a control | 40/40 flagged |
| Copied trace file | 40/40 |
| Re-typed copy under another animal | 40/40 |
| Misfiled clone, label intact | 40/40 |
| Whole lineage misfiled, labels intact | 125/125 clones |
| Misfiled clone, label changed too | 22/40; see below |
| Mutations lost to the germline rule, flagged | 9/9 |
| Flags on correctly filed clones | 0 of 6,018 |

**The detection limit.** When the label has been changed too, the only evidence is shared mutations. Detection depends on how many mutations the clone shares with its lineage mates still in its true animal:

| Shared with mates | 3 | 4 | 5 | 6 | 7 | 8+ |
|---|---|---|---|---|---|---|
| Detected | 0/3 | 0/6 | 4/10 | 9/10 | 6/7 | 3/4 |

The misfiled clones in your real data shared 10–21 mutations.

`tests/test_sanger_validation.py` runs two of these projects, plus a clean one with nothing planted, on every test run. It fails if calling or QC gets worse.

### Your real data (`tests/test_real_jh4_sanger.py`)

`docs/calibration/samplesheets/jh4_manuscript.csv` is the manuscript's six mice, with every clone kept where it was filed. Run through the pipeline, it gives the manuscript's numbers exactly:

| | WT2 | WT3 | WT4 | Aicda⁻/⁻ 1 | Aicda⁻/⁻ 2 | Aicda⁻/⁻ 3 |
|---|---|---|---|---|---|---|
| Clones | 39 | 38 | 27 | 11 | 10 | 10 |
| Unique calls (101–629) | 153 | 230 | 111 | 1 | 1 | 2 |

The test runs where the traces are and is skipped elsewhere. On this sheet the QC reports, without changing anything:
- **Duplicates (problem):** three clones labelled `wt3-3`, `wt3-4` and `wt3-6`, filed under WT2, are identical to WT3 clones (sharing 10–21 mutations).
- **Label mismatches (check):** 14 clones labelled w3 are in WT2, carrying 1–26 mutations. 13 clones labelled w2 are in WT3: 12 carry no mutation and one carries 3. Two Aicda⁻/⁻ clones carry another knockout mouse's label.
- **Clean controls:** no Aicda⁻/⁻ clone carries more mutations than its qualities explain.
- **Germline:** the same five germline positions are removed in all six mice, so they are true strain differences, and no lineage was lost to the rule.

With the MKLN mice added (`jh4_with_mkln.csv`), the QC also finds that WT4 clone w4-9 (434141-A11) is identical over 520 bp, with 7 shared mutations, to two MKLN K2 clones (434141-E11 and 434355-H3).

## Known limits

- **A whole lineage filed under the wrong animal** can't be seen in the sequence, because no clone of it is left in the true animal to share mutations with. Only the label check catches it, and only when labels and aliases are given.
- **The germline rule:** a lineage that fills half of an animal's clones, or a hotspot mutated in several lineages, loses its mutations. Heterozygous germline differences look the same to the rule. The flag says when this may have happened; deciding is up to you.
- **The Q20 error floor:** single-pass Sanger at Q20 leaves about 2–3 false calls per 10⁴ bases read. Knockout controls sequenced alongside measure this floor for a given run.
