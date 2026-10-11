#!/usr/bin/env bash
# Reproduce the bundled demo: wild type vs Aicda-/- in the JH4 intron.
# Expect a large AID-dependent difference in unique variant count.
set -euo pipefail
cd "$(dirname "$0")/.."

for s in WT3 KO3; do
  nanoshm \
    --fastq "examples/${s}_JH4_sub500.fastq" \
    --reference examples/JH4_reference.fasta \
    --contaminant examples/JH4_contaminant.fasta \
    --preset jh4_intron \
    --out "results/${s}" --plots
done

echo
echo "Results in results/WT3 and results/KO3 (summary.json, variants.csv, PNGs)."
