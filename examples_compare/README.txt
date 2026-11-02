SHMtrace test dataset (mouse JH4 intron, Nanopore R10.4.1, 500 reads each)

Six libraries, one FASTQ per mouse:

  WT2.fastq WT3.fastq WT4.fastq   wild type, AID-proficient: the samples
  KO1.fastq KO2.fastq KO3.fastq   Aicda-/-, AID-null: the controls

The Aicda-/- mice cannot carry out somatic hypermutation, so whatever is called
in them is error. SHMtrace uses them to fit its error model and to set the
calling threshold, which is why they are marked control=yes in sheet.csv.

  Reference:   JH4_reference.fasta
  Contaminant: JH4_contaminant.fasta
  Preset:      jh4_intron

The reads are a Geneious export of reads assembled to J꜀4, so Geneious has
already filtered them and turned every read to the reference: they all show as
forward, and the run notes say so. They are kept for the tests and the benchmark
against Geneious. The app's nanopore sample data is the same six libraries as
sequenced, in examples_nanopore/.

Two sheets, for two different tools:

  sheet.csv    the SHMtrace sample sheet (file,sample,group,control), read by
               shmtrace-sheet
  samples.csv  the benchmark sheet (sample,fastq,null,run,geneious_calls), read
               by shmtrace-benchmark; its geneious_calls column holds the
               Geneious Prime counts the manuscript compares against

Run it:

  shmtrace-sheet run sheet.csv \
    --reference JH4_reference.fasta \
    --contaminant JH4_contaminant.fasta \
    --region jh4_intron -o results

Then open results/report.html.

Expected (SHMtrace 0.2.2):

  sample  genotype    unique mutations
  WT2     wild type   180
  WT3     wild type   362
  WT4     wild type    68
  KO1     Aicda-/-      0
  KO2     Aicda-/-      2
  KO3     Aicda-/-      0

Separation between the groups is about 247-fold by mutation frequency and
308-fold by unique-variant density, with the wild-type spectrum strongly biased
to G:C and strong AID-hotspot targeting.

What the controls buy you: run KO3 on its own, with nothing to calibrate
against,

  shmtrace --fastq KO3.fastq --reference JH4_reference.fasta \
    --contaminant JH4_contaminant.fasta --preset jh4_intron --out results/KO3_alone

and it reports 5 variants instead of 0. Those five are sequencing error that
looks significant on base qualities alone. (Note --out here: the single-sample
command takes --out, while shmtrace-sheet takes -o.)
