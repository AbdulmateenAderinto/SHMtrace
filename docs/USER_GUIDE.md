# nanoshm user guide

## What the app does

nanoshm takes the read files from a sequencing run, a reference sequence, and
your note of which samples are controls, and gives you the mutation tables and
figures for a somatic hypermutation experiment. It replaces the part of the work
that is done by hand in Geneious: mapping the reads to the reference, separating
the contaminant amplicon, trimming and quality-filtering, equalising the number
of reads per sample, calling the mismatches that are more than sequencing noise,
removing germline positions, and then counting densities, spectra and hotspot
targeting in Prism. One run of the app does all of that for every sample at
once, and writes the tables and the figures. It reports the same quantities as
Aderinto & Martin. Indels are ignored throughout; this is a substitution assay.

## Getting started

### The built app

There is a built application for macOS and for Windows. It contains everything
it needs, including Python, so nothing else has to be installed.

- **macOS:** drag `nanoshm.app` to Applications. The first time you open it,
  right-click (or control-click) the icon and choose **Open**, then confirm.
  The app is not signed by Apple, so a double-click on the first launch is
  refused. After that first time it opens normally.
- **Windows:** unzip the folder anywhere and run `nanoshm.exe` inside it. Keep
  the whole folder together; the `.exe` needs the files beside it.

The app opens one window. Inside that window is a page served by the app itself
on `127.0.0.1`, which is your own machine. Nothing listens on the network and no
reads, references or results are sent anywhere. Everything stays on the computer
you run it on, including the files you drop in, which are copied into a temporary
folder for the length of the session.

### Installing with pip instead

If you would rather run it from a Python environment (Python 3.10 or newer):

```bash
pip install -e ".[app]"
nanoshm-app
```

`nanoshm-app` starts the same local server and opens the same window. If no
window toolkit is available, it prints an address such as
`http://127.0.0.1:52431/` and keeps serving; open that address in a browser.
`nanoshm-app --no-window` does that deliberately, which is how to use the app on
a machine you are logged into over SSH.

The `[app]` extra brings in the window toolkit, the Excel writer and the BAM
reader. `pip install -e ".[web]"` gives a working app without those three, so
the `.xlsx` export will not be available.

## A worked experiment, step by step

The example is six JH4 intron libraries from one Plasmidsaurus order: three
wild-type mice and three *Aicda*−/− mice, one FASTQ per mouse.

### Set up

**1 Region and data type.**

- **Region:** `jh4_intron`. The hint under the menu says what the preset is for.
- **Sequencing platform:** `ont`.
- **Calling mode:** `recommended`, which is already selected.
- Leave the **Advanced** box closed. Its three fields (reads per sample,
  analysis window, significance threshold) override the preset, and the preset
  is right for this assay.

**2 Reference.**

Open **Built-in reference** and choose *the JH4 intron entry (1,230 bp)*. The
contaminant line below then fills in by itself and says
`(added for you)`: picking the built-in JH4 intron also selects the bundled JH4
contaminant, which is the sequence that has to be competed away for this
amplicon. If you have your own reference instead, leave the menu on
**Use my own file…** and drop your FASTA on the box below it, and drop your own
contaminant FASTA on the second box. Only the first record of a reference FASTA
is used.

**3 Samples.**

Drop all six read files on the sample box, or click it and select them. Each
file you add becomes one sample. The names are taken from the file names with
the extension removed, so `JH4_WT_1.fastq.gz` becomes `JH4_WT_1`.

Now press the two buttons:

- **Group by name** fills the group box of each sample. It splits each name on
  underscores, dots, hyphens and spaces, throws away the parts that every
  sample shares and the parts that are only digits, and uses what is left. With
  these six names it gives `WT` for three samples and `Aicda KO` for the other
  three. Check the result and type over anything it got wrong; the group box
  accepts any text, and samples with the same text are one group.
- **Controls by name** ticks the **control** box of any sample whose name
  contains `ko`, `knockout`, `empty`, `aicda`, `null`, `control`, `ctrl`,
  `germline`, `untransfected` or `mock`. Here it ticks the three *Aicda*−/−
  samples. It matches anywhere in the name, so read down the list and correct
  the boxes by hand if it has ticked something it should not have.

The count at the right should read `6 samples, 3 controls`.

**4 Run.** Type a project name, for example `JH4 WT vs Aicda KO Oct 2026`, and
press **Analyse**. A progress bar appears and the messages run through
`aligning and filtering reads`, `calibrating the background model`,
`calling <sample>` for each sample in turn, and `summarising`. Six libraries of
about 5,000 reads each take a few minutes. When it finishes, the three other
tabs become usable and the page moves to Results.

### Results

At the top is a summary block: the region, reference name and analysed window;
the platform, mode and the number of reads per sample the run settled on; how
many samples, groups and controls there were; and the mean mutation density of
the non-control samples against the mean of the controls, with the fold
difference. For a healthy JH4 experiment the wild-type mean should be tens of
mutations per 10⁵ bp and the *Aicda*−/− mean should be close to zero.

Below are seven figures. They are interactive: hovering over a point names the
sample, and the camera button in the corner saves the figure as a PNG.

- **Mutation density**, unique mutations per 10⁵ bp, one dot per sample with the
  group mean and standard deviation.
- **Mutation frequency**, mutations per 10³ bases sequenced.
- **Substitution spectrum**, the twelve substitution classes as a percentage of
  each group's unique mutations, with the G:C-origin and A:T-origin halves
  shaded.
- **Mutations per read**, the distribution of mutations per read for each group.
- **Position along the amplicon**, the percentage of reads mutated at each
  position, with WRC/GYW hotspot positions marked along the bottom axis.
- **AID hotspot targeting**, the percentage of each sample's G:C mutations that
  sit in a WRC or GYW motif.

### Tables

The menu at the top left picks between five tables: Samples, Groups and
statistics, Every mutation, Mutations per position, and Substitution spectrum.
The filter box matches text anywhere in a row. Clicking a column heading sorts
by it, and clicking again reverses. **Download CSV** saves the table you are
looking at. The table shows at most 3,000 rows; the CSV has all of them.

### Export

**Everything (.zip)** is the one to take. It contains the six CSV tables, the
seven figures as PNG, a standalone `report.html`, `project.json` with every
number and the settings that produced them, and `tables.xlsx` if the Excel
writer is installed.

Single buttons give each table on its own, and each figure as PNG, SVG or PDF.
In the built app only, a **Save to a folder on this computer** box appears at
the bottom: type a folder path and everything is written there directly instead
of going through the downloads folder.

## Choosing the settings

### Region

The region preset carries everything specific to the assay: how long an
alignment has to be, the analysis window, the read depth, and how much is
trimmed from the ends.

| Region | What it is for | Analysis window | Reads per sample | Minimum alignment |
|---|---|---|---|---|
| `jh4_intron` | Mouse JH4 downstream intron, the neutral polyclonal assay. Germline is unknown, so high-frequency positions are removed by the under-50% rule | 101–629 of the 1,230 bp reference | 500 | 800 bp |
| `n114_vregion` | N114 hybridoma V region; the parental sequence is known, so germline is defined exactly. Excess mutation over empty-vector controls is AID activity | 21–635 of the 655 bp reference | 500 | 400 bp |
| `igh_smu` | Mouse IgH Sμ switch region. Heavily mutated, so mutation loads per read are high | 21–631 | 2,000 | 420 bp |
| `generic_amplicon` | Any other region: JH6, Sγ, a V gene, an AID off-target. Everything except 20 bp at each end is analysed and nothing is assumed about the germline | whole reference minus 20 bp each end | 1,000 | 200 bp |
| `sanger_clones` | Sanger traces, one clone per read, with no statistical test | whole reference minus 20 bp each end | all reads | 200 bp |

For a new region, start with `generic_amplicon`, look at one run, and then
narrow the window with the Advanced box. A preset is a small YAML file, so a
permanent new assay is a new file rather than new code.

### Platform

The platform leaves the assay settings alone and changes the aligner and the
read filters.

| Platform | Aligner | What changes |
|---|---|---|
| `ont` | `map-ont` | Oxford Nanopore. The preset as written. The default |
| `pacbio` | `map-hifi` | HiFi reads. Multiple-testing correction counts every analysed position times every alternative base, because accurate reads show only a few error singletons |
| `illumina` | `sr` | Read pairs are merged on their overlap, at least 20 bp, by quality-aware consensus; pairs that do not overlap are dropped. The quality gate is relaxed to 0.02 expected errors per base and the minimum alignment is lowered to 100 bp. Same correction change as PacBio |
| `sanger` | `map-ont` | One trace is one clone. Modified-Mott quality trimming at Q20, quality gate 0.05 on the trimmed read, minimum alignment 200 bp, no depth equalising and no phasing |

Short reads cannot span a long-read amplicon, so with `illumina` set an analysis
window in the Advanced box that your reads actually cover.

### Calling mode

The mode decides how a mismatch is judged against sequencing error.

| Mode | What it does | When to use it |
|---|---|---|
| `recommended` | Quality-calibrated against your controls. How often a base of a given quality is a miscall is measured in the control reads, per substitution and strand, and the measured rate is inflated by a factor calibrated on the controls themselves | Always, unless you have a reason not to. The default |
| `geneious` | Geneious's own approximate p-value, uncorrected, at p ≤ 0.05, with no phasing | Reproducing Geneious numbers for a comparison or a revision |
| `v0.1` | The earlier default: a per-position error rate from the controls, Benjamini-Hochberg corrected | Comparing with results produced before October 2026 |
| `strand` | Each strand tested against its own control error rate, with overdispersion between libraries, combined by Fisher's method | A run where errors appear on one strand only |
| `pbinom` | An exact test on the base qualities alone | No controls at all on the run |
| `qcal+filter` | `recommended` plus a filter that drops calls whose reads are split between the strands far more unevenly than the rest of the column | A run where you suspect strand artefacts but still want the default's behaviour |
| `sanger` | Every mismatch with base quality 20 or better is a mutation, with no statistical test | Sanger traces. The platform picks it for you, so `recommended` already means this |

`recommended` is the default because of the benchmark in
`docs/calibration/ERROR_MODELS.md` (11 October 2026, on the paper deposit). In
the wild-type JH4 libraries it agrees with Geneious at the level of individual
sites (Jaccard overlap 0.85 to 0.97) and gives essentially the same density,
73.4 against Geneious's 74.9 per 10⁵ bp, while its background in the *Aicda*−/−
libraries is a quarter of Geneious's, 0.38 against 1.52 per 10⁵ bp. The
separation between wild type and knockout is therefore 192-fold instead of
49-fold. On the same document's spike-in test, where known mutations are planted
into control reads, it made 0.0 false calls per JH4 sample and 0.3 per N114
sample, against Geneious's 2.7 and 6.0; the price is lower recall of mutations
carried by a single read (70% against Geneious's 93% on JH4).

Two things follow from this that matter for interpretation. The modes differ
most on single-read calls, which is where Geneious's false calls come from. And
no mode can confirm a single-read mutation at 500 reads without some false
calls, so a single-read difference between two samples is not evidence on its
own.

## Controls

A control is a sample that should carry no somatic mutations: an *Aicda*−/−
mouse, an empty-vector transfection, a germline or untransfected sample. Tick
the **control** box for each one. Everything the app knows about the sequencing
background comes from those samples.

Every sample needs controls from the same sequencing run. Nanopore error depends
on the flow cell, the chemistry and the basecaller version, and the basecaller
can change between orders without being announced, so a background measured on
one run does not describe another. Put a mutation-free sample on every flow
cell. That is the same argument for running an AID-null control alongside every
experiment that the manuscript makes.

Each control is measured against the other controls, never against itself, so a
control's own mutation count is an honest estimate of the background rather than
zero by construction. The **null from** column of the samples table names the
controls each sample was measured against.

Use at least two controls, and three if you can. With one control only, there is
no second control to calibrate against, and the measured rates are used as the
quality scores give them, which the calibration document shows is too permissive
(about 29.5 false calls per JH4 sample). With no controls at all, the
quality-calibrated model has nothing to measure and falls back silently to a
test on the base qualities alone. The run still completes and still produces
numbers, but they carry more background and the only sign of the fallback in the
app is `none` in the **null from** column. Check that column after every run.

A group counts as the control group for the statistics only if every sample in
it is ticked as a control. If no group is entirely controls, the Groups table
has means and standard deviations but no p-values.

## Reading the results

### The figures

| Figure | What it shows |
|---|---|
| Mutation density | Unique mutations per 10⁵ analysed bases, one dot per sample, bar at the group mean with the standard deviation. The headline figure of the assay |
| Mutation frequency | Mutations per 10³ bases sequenced, which counts a mutation once per read that carries it rather than once per site |
| % of reads mutated | The share of reads carrying at least one mutation (in the exports; the Results tab does not draw it) |
| Substitution spectrum | The twelve substitution classes as a percentage of each group's unique mutations |
| Position along the amplicon | The percentage of reads mutated at each reference position, with WRC/GYW positions marked |
| Mutations per read | The distribution of mutations per read, which separates many lightly mutated reads from a few heavily mutated ones |
| Hotspot targeting | The share of G:C mutations in WRC or GYW motifs, one dot per sample |

### The samples table

One row per sample.

| Column | Meaning |
|---|---|
| `sample`, `group`, `control` | As you entered them on the Set-up tab |
| `reads` | Reads actually analysed, after filtering, phasing and depth equalising |
| `unique_mutations` | Number of distinct called mutations, counting a site and base change once however many reads carry it |
| `unique_density_x1e5` | `unique_mutations` divided by (analysed bases × mean depth), multiplied by 10⁵. This is the paper's unique-variant density. It saturates with depth, so it is only comparable between samples analysed at the same depth |
| `mutation_frequency_x1e3` | Total mutations summed over all reads, divided by (reads × analysed bases), multiplied by 10³. A mutation in 40 reads counts 40 times. Unlike the density this does not saturate |
| `mutated_reads_pct` | Percentage of analysed reads carrying at least one mutation |
| `mean_mutations_per_read` | Total mutations divided by reads |
| `pct_GC` | Percentage of the unique mutations whose reference base is G or C. AID targets G:C, so a wild-type library should sit well above half |
| `pct_transitions` | Percentage of the unique mutations that are transitions, that is A to G, G to A, C to T or T to C |
| `pct_GC_in_WRC_GYW` | Of the mutations at G:C bases, the percentage in a WRC or GYW motif, the canonical AID hotspot. The denominator is the G:C mutations, not all mutations |
| `pct_AT_in_WA_TW` | Of the mutations at A:T bases, the percentage in a WA or TW motif, the Pol-η signature. The denominator is the A:T mutations |
| `analysed_bp` | Width of the analysed window after the end trims |
| `null_from` | The controls this sample's background was measured from, or `none` |
| `warnings` | Messages from the run, most often that fewer reads were available than the target depth |

A sample with no reads shows `0` reads and empty cells in the derived columns.

### The groups table

One row per group, with `n`, whether it is the control group, and for each of
`unique_density_x1e5`, `mutation_frequency_x1e3`, `mutated_reads_pct`,
`mean_mutations_per_read`, `pct_GC`, `pct_transitions` and
`pct_GC_in_WRC_GYW`:

- `..._mean` and `..._sd`, the mean and sample standard deviation across the
  samples of the group. A group of one sample has a mean and no standard
  deviation.
- `..._p_welch_vs_<control group>`, Welch's unequal-variance t-test of this
  group against the control group. It needs at least two samples in each group,
  and is blank otherwise.
- `..._p_mwu_vs_<control group>`, the two-sided Mann-Whitney U test of the same
  comparison. It makes no assumption about the shape of the distribution, but
  it is also limited by how few orderings three against three allows: with
  three samples per group the smallest two-sided p it can ever return is 0.1,
  even when the two groups do not overlap at all. A Mann-Whitney p of 0.1 from a
  3-versus-3 design means "as separated as this design can show", not "not
  significant". With three per group, report the Welch p, or the effect size
  with every sample shown as a dot, and keep the Mann-Whitney as a check.

Each sample here is one mouse or one transfection, not one read, so these tests
are over animals, which is the comparison that matters.

### The other tables

- **Every mutation**: one row per called mutation per sample, with the position,
  reference and alternative base, supporting reads, depth at that position, the
  frequency, the p-value, the plus- and minus-strand read counts, the 5-mer
  context, the homopolymer run length, and whether the position is a WRC/GYW
  hotspot. A p-value displayed as `<1e-308` is smaller than the computer can
  represent, not zero.
- **Mutations per position**: for each group and position, how many reads were
  mutated there and what percentage of the group's reads that is. This is the
  table behind the positional figure.
- **Substitution spectrum**: counts per group per substitution class, plus
  percentages in the CSV.
- **Mutations per read** (`mutation_load.csv`, in the Export tab only): for each
  sample, how many reads carried 0, 1, 2 and so on mutations.

## Exporting for a paper

| Export | What is in it |
|---|---|
| `bundle.zip` ("Everything") | All six CSV tables, all seven figures as PNG, `report.html`, `project.json`, and `tables.xlsx` if the Excel writer is installed |
| `report.html` | One self-contained file: the figures, the samples and groups tables, the first 2,000 mutations, and the settings that produced them. It opens in any browser and can be emailed as it is. Good for a lab meeting or a supplementary file |
| `tables.xlsx` | One workbook: a Run sheet with the date and the settings, then Samples, Groups, Variants, Positional and Spectrum sheets |
| Individual CSVs | `samples.csv`, `groups.csv`, `variants.csv`, `positional_frequency.csv`, `spectrum.csv`, `mutation_load.csv` |
| Figures | Each figure as PNG, SVG or PDF |
| `project.json` | Every number in the result, plus the settings. The machine-readable record to deposit alongside a paper |

For a figure you will edit in Illustrator, Inkscape or Affinity, download the
**SVG** or the **PDF**. Both keep the text as text and the lines as lines, so
fonts, labels and colours can be changed and nothing goes blurry when it is
scaled. Use the PNG for slides and for a quick look. The PNGs are also what the
ZIP and the HTML report contain, so download the SVG or PDF separately for the
figures that will go into the manuscript.

If you prefer to redraw the figures yourself in Prism, the CSVs hold the
underlying numbers: `samples.csv` for the per-sample dots, `spectrum.csv` for
the spectrum, `positional_frequency.csv` for the map and `mutation_load.csv` for
the mutations-per-read distribution.

## Sanger, Illumina and PacBio data

A sample can be one file, a whole folder, or a pair of Illumina files. The app
works out which from what you drop in.

### Sanger traces

Set the platform to `sanger` and the region to `sanger_clones` (or your usual
region preset). Leaving the mode on `recommended` is correct: with the Sanger
platform, `recommended` means the Sanger caller.

`.ab1` and `.abi` traces are read directly, including the per-base quality
values stored in the file, and `.seq` text files are read too. Each read is
trimmed to its best-quality segment, and every mismatch of quality 20 or better
inside the window is a mutation, with no statistical test, because a clone's
mutation is seen in exactly one read by construction.

When you drop a folder of traces on the sample box, the whole folder becomes **one sample**, which is what a library of clones is. The Choose a folder button does the same. Drop the individual files instead if you want each trace treated as its own sample.

```bash
nanoshm --fastq sanger_clones/ --reference ref.fasta --preset sanger_clones \
        --platform sanger --out results/library1
```

Depth equalising and phasing are both off under `sanger`, so the number of
clones per library is whatever you sequenced. Compare libraries with similar
numbers of clones, or compare the per-clone frequency rather than the density.

**Give it the other haplotype.** If your animals or cells can carry a second
allele of the locus, pass it with `--contaminant` alongside your reference.
On thousands of nanopore reads nanoshm finds such a haplotype by itself,
because the same block of differences recurs across reads; on a few dozen
clones there is nothing for it to recur across, and a clone from the wrong
allele is then read as one carrying forty or more mutations. On the paper's
JH4 clones, adding the Agouti allele removed exactly the nineteen sequences
the authors removed by hand and left the same ninety-two
(`sanger_vs_ont/README.md`).

### Illumina

Set the platform to `illumina`. Paired reads are expected as an R1 and an R2
file whose names differ only in the read-number token, written `_R1`/`_R2` or
`_1`/`_2`, with `.`, `-` or `_` before it, for example
`lib1_S1_L001_R1_001.fastq.gz` and `lib1_S1_L001_R2_001.fastq.gz`. Each pair is
merged on its overlap into a single read, which needs at least 20 bp of overlap;
pairs that do not overlap are dropped, so a long insert with a short read length
loses most of its reads.

Drop both files of a pair on the sample box and the app pairs them for you: it
recognises the read-number token, merges them into one sample and marks the row
"paired reads". A lone R1 with no R2 is treated as a single-end sample. On the
command line, write the pair as one argument:

```bash
nanoshm --fastq lib1_R1.fastq.gz+lib1_R2.fastq.gz --reference ref.fasta \
        --preset generic_amplicon --platform illumina --out results/lib1
```

If your reads are single-end, or already merged, one file per library is all the
app needs and nothing special applies.

Set an analysis window the reads actually cover. A 2 × 150 bp run cannot span
the 529 bp JH4 window, and reads that do not reach the window contribute
nothing. The quality gate is looser than for nanopore (0.02 expected errors per
base rather than 0.01) because read tails are poor while the merged overlap is
not, and multiple-testing correction counts every analysed position times every
alternative base, which it must on accurate reads.

### PacBio

Set the platform to `pacbio`. HiFi reads are long enough for the full window and
clear the quality gate easily, so nothing else changes except the aligner and
the same correction change as for Illumina. One BAM or FASTQ per sample. BAM,
SAM and CRAM need the BAM reader, which the built app and the `[app]` install
already include; secondary and supplementary records are skipped, reverse reads
are turned back to the orientation they were sequenced in, and each read is used
once.

Geneious's own exported assemblies are BAM files, so an old Geneious analysis can
be re-run through nanoshm on exactly the reads Geneious used.

## Troubleshooting

**A sample's `warnings` column says "only 312 reads available for a target depth
of 500".** The sample had fewer usable reads than the run was aiming at, so it
was analysed at what it had. Unique density saturates with depth, so that
sample's density is not comparable with the others. Either leave that sample out,
or set **Reads per sample** in the Advanced box to a number every sample can
reach (the shallowest sample's read count) and run again. Note that when you
leave the box blank the app already lowers the common depth to the shallowest
sample, so this warning normally means you typed a depth yourself.

**One shallow sample has dragged every other sample's depth down.** That is
deliberate: with the box blank the common depth is the smaller of the preset's
default and the shallowest sample's usable reads, because densities at different
depths cannot be compared. Remove the shallow sample and run again if you would
rather keep the depth.

**A sample shows 0 reads and empty metric columns.** Nothing survived the
filters. The usual causes, in order of likelihood: the reference is not the
amplicon that was sequenced, so nothing aligns; the reads are shorter than the
preset's minimum alignment length (800 bp for `jh4_intron`), for example because
they are Illumina or Sanger reads analysed with the platform left on `ont`; or
the reads are too noisy for the quality gate. Check the platform first, then
whether the reference is right, then try `generic_amplicon`, which asks for only
200 bp of alignment.

**The `null from` column says `none` for every sample.** No control was ticked,
or no control other than the sample itself. The quality-calibrated model had
nothing to measure and quietly fell back to a test on the base qualities alone,
which carries more background. Nothing in the app says so beyond that column.
Go back to Set up, tick the controls, and run again. If the run genuinely has no
mutation-free sample, use `pbinom`, which is designed for that case, and say so
when you report the numbers.

**Only one control is ticked.** The run works, but the calibration step needs at
least two controls, so the measured error rates are used uninflated, which the
calibration document shows is too permissive. Tick two, and preferably three.

**The Groups table has no p-values.** The statistics need a group in which every
sample is ticked as a control. If your controls are spread across groups, or
only some samples of the control group are ticked, there is no reference group to
test against. Tick every sample of the control group.

**A group's `..._sd` columns are empty and it has no p-values.** That group has
one sample. A standard deviation needs two, and Welch's test needs two in each
group. Mann-Whitney will still report a value for a group of one, but it cannot
go below 0.5 for one against three, so do not read anything into it.

**Every Mann-Whitney p is 0.1 and nothing looks significant.** With three
samples per group that is the smallest two-sided p the test can produce, however
wide the separation. Report the Welch p instead, and show the individual samples
as dots.

**Short reads cannot span the window.** With `illumina`, if most reads are
dropped or the densities are implausibly low, the analysis window is probably
wider than the reads reach. Set a window in the Advanced box inside the region
the reads cover, written as `start-end`, for example `200-400`.

**The reference is the wrong length for the preset.** Each preset's window is
written in the coordinates of its own reference: 101-629 for the 1,230 bp JH4
intron, 21-635 for the 655 bp N114 V region. Give a shorter reference and the
window is clipped to what exists, which can leave a very narrow analysed region
or, with a reference shorter than the start of the window, nothing sensible at
all. The headline block on the Results tab shows the window and the reference
length; check that line first. Either use the matching reference or switch to
`generic_amplicon` and set the window yourself.

**The `.xlsx` download fails or `tables.xlsx` is missing from the ZIP.** The
Excel writer (`openpyxl`) is not installed. The built app and
`pip install -e ".[app]"` include it; `pip install -e ".[web]"` does not. Either
install it (`pip install openpyxl`) or use the CSV files, which contain the same
numbers.

**"Two samples have the same name."** Sample names must be unique, and they come
from the file names with the extension removed, so two files of the same name
from different folders collide. The run refuses to start. Edit the group boxes
all you like, but the names themselves are fixed once the file is added, so
remove one of the two with the ✕, rename the file, and add it again.

**The contaminant box is empty after you dropped your own reference.** The
bundled JH4 contaminant is added for you only when you pick the built-in JH4
intron from the menu. With your own reference FASTA, drop the contaminant FASTA
on the second box yourself. Without it, contaminant reads compete for the
alignment and inflate the background.

**Clicking the sample box opens a file picker and not a folder picker.** To add
a folder, drag it onto the box; clicking only ever selects files.

**The progress bar reaches the end but the page shows a red message.** The
message is the error the analysis raised. The most common are a reference file
that is not FASTA, a read file that is none of the formats the app recognises,
and duplicate sample names.

## What the app does not do

- **Indels.** Insertions and deletions are excluded at every step. Only
  substitutions are called, counted and plotted. An indel-rich amplicon needs a
  different tool.
- **Clonality and lineage.** There is no clone calling, no lineage tree and no
  mutation-sharing analysis. The app flags a block of co-occurring positions
  that could be a second haplotype or an expanded clone, and leaves it in unless
  the same block also appears in a control or in another animal, but it does not
  tell you how many clones there are or how they are related.
- **AIRR export.** There is no AIRR Community rearrangement output, and no V, D
  or J gene assignment. The reference is whatever FASTA you supply, not a
  germline database.
- **UMI-based template counting.** Reads are not grouped by molecular barcode,
  so the counts are read counts, not template counts, and PCR duplicates are not
  collapsed.
