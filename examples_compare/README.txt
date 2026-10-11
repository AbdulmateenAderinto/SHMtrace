SHM-webtool test dataset (mouse JH4 intron, Nanopore R10.4.1, 500 reads each)

COMPARE-GROUPS test:
  CONTROL (wild-type, AID-proficient): WT2.fastq WT3.fastq WT4.fastq
  CASE    (Aicda-/-, negative control): KO1.fastq KO2.fastq KO3.fastq
  Reference:   JH4_reference.fasta
  Contaminant: JH4_contaminant.fasta
  Preset: jh4_intron

Expected: wild-type shows ~7x higher, highly significant mutation frequency
and strong AID-hotspot targeting vs the Aicda-/- knockout.
