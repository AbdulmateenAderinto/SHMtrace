# Pooled runs and screens

*12 October 2026. `shmtrace/demux.py`, `shmtrace-screen`; tests in `tests/test_demux.py`.*

## Inline barcodes

A screen of SHM regulators (CRISPR knockouts in Ramos or CH12, or a panel of hybridoma subclones) is cheapest as one pooled run. An inline barcode in each PCR primer tells the conditions apart: `<barcode><primer> … <primer′><barcode2>`. ONT native barcodes are split by the basecaller; inline barcodes are split here.

Each read is oriented by finding the left primer anchor near its start, in the read or its reverse complement. The bases just outside the anchor are its barcode. A read is assigned to a sample only if all three of these hold:
- the best barcode is within 2 edits;
- it is at least 2 edits better than the runner-up;
- with dual barcodes, both ends name the same sample.

A doubtful read is left unassigned rather than put in the wrong sample.

In the sample sheet, give every row the same pooled FASTQ in `file`, plus `barcode` (and `barcode2` for dual barcodes):

```
file,sample,group,control,barcode,barcode2
pool.fastq,EV_1,EV,yes,ACGTACGTAC,AGCTTCGACA
pool.fastq,EV_2,EV,yes,TTGCAAGCTG,CTAGGATCGT
pool.fastq,geneX_1,geneX,no,GGATCCTAGA,TCCAGGTTAC
...
```

The run is split into one FASTQ per sample before the analysis. By default the anchors are the first and last 20 reference bases; for primers outside the reference, give `left_anchor` and `right_anchor`.

### Accuracy on simulated pooled runs

Four samples, 250 reads each, 10-nt barcodes, half the reads reverse-complemented, and substitution and indel errors throughout:

| Read error | Barcodes | Assigned correctly | Assigned wrongly | Unassigned |
|---|---|---:|---:|---:|
| 3% | single | 98.3% | 0.00% | 1.7% |
| 3% | dual | 96.8% | 0.00% | 3.2% |
| 8% | single | 84.0% | 0.00% | 16.0% |
| 8% | dual | 72.8% | 0.00% | 27.2% |

Errors cost reads, not accuracy: no read went to the wrong sample. Barcode sets should differ by at least 4 edits pairwise; the margin rule relies on that.

**In the app:** load a sample sheet with `barcode` columns (Set up → Samples). The pooled file is split when you press Analyse, and a note reports the reads assigned to each sample.

## The screen report

`shmtrace-screen project.json -o screen/` compares every condition with the control group using the animal-level tests (`REPLICATES_AND_POWER.md`). Each replicate culture or mouse is one unit. It reports:
- the fold change with its 95% CI;
- the p-value;
- a Benjamini–Hochberg q-value across the conditions.

Conditions are ranked by q-value. It writes `screen.csv` and a volcano plot (`screen.png`: log2 fold against −log10 q, with conditions past the q cut-off labelled).

Plan replicates with `shmtrace-power`. With between-replicate variation like the J꜀4 mice, a 2-fold change needs many replicates per condition. Cell-line replicates usually vary much less than mice, so a pilot of the control condition is the way to find out.

**In the app:** the *Screen* tab, after an analysis. It ranks every condition against the control group, marks folds against a control with no mutations (†: one event per control animal, a zero-count convention, not a confidence bound), draws the volcano and downloads the table.
