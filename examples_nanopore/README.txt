SHMtrace nanopore sample data (mouse JH4 intron, Nanopore R10.4.1, reads as sequenced)

Six libraries, one gzipped FASTQ per mouse, exactly as the sequencing service delivered
them (run DHJGKB): about 5,000 reads each, in both orientations, unfiltered.

  WT2.fastq.gz WT3.fastq.gz WT4.fastq.gz   wild type, AID-proficient: the samples
  KO1.fastq.gz KO2.fastq.gz KO3.fastq.gz   Aicda-/-, AID-null: the controls

Nothing has been done to these reads. A library as sequenced also holds reads
too short or too noisy to use, reads of the JH4-like paralogue the primers also
amplify, and other off-target products; finding and setting those aside is part
of what SHMtrace does, and the Reads tab shows what became of each one.

  Reference:   SHMtrace's bundled mouse_jh4_intron.fasta
  Contaminant: SHMtrace's bundled mouse_jh4_contaminant.fasta (the deposit's
               contaminant and the JH4-like paralogue)
  Preset:      jh4_intron (positions 101-629, 500 reads analysed per library)

In the app: Set up -> Try it with nanopore sample data. From the command line:

  shmtrace-sheet run samples.csv \
    --reference ../src/shmtrace/references/mouse_jh4_intron.fasta \
    --contaminant ../src/shmtrace/references/mouse_jh4_contaminant.fasta \
    --region jh4_intron -o results

Expected (SHMtrace 0.2.2): about 3,000 reads of each library pass the length
and quality gates; contaminant and off-target reads go next, leaving roughly
1,000-2,800 usable J꜀4 reads, of which 500 are drawn for the analysis (every
call also gives its count over all of them).

  sample  genotype    unique mutations  per 10^5 bp
  WT2     wild type   221                83.7
  WT3     wild type   344               133.2
  WT4     wild type    73                27.8
  KO1     Aicda-/-      0                 0
  KO2     Aicda-/-      2                 0.8
  KO3     Aicda-/-      0                 0

Wild type averages 81.6 unique mutations per 10^5 bases against 0.25 in
Aicda-/- (Welch's test on log density, p = 0.003).

examples_compare/ holds the 500-read subsets of the same libraries that the
manuscript analysed (exported from Geneious after it had assembled them to JH4,
so already filtered and turned to the reference). The tests and the benchmark
against Geneious use those; the app's sample data is this folder.
