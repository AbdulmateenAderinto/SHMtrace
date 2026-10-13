#!/bin/zsh
# Chain A: the manuscript analyses with the control-based FDR (recommended mode)
cd ~/Nanoshm/nanoshm_v0.2
P=.venv/bin/python
set -x
$P docs/calibration/jh4_shmtrace.py
$P docs/calibration/jh4_shmtrace.py --full-length
$P docs/calibration/n114_shmtrace.py
$P docs/calibration/depth_shmtrace.py jh4
$P docs/calibration/depth_shmtrace.py n114
$P docs/calibration/sanger_vs_nanopore.py --ont-mode recommended --sanger-dir ~/Nanoshm/sanger_ab1_original --suffix _manuscript --ont-dir "/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ"
$P docs/calibration/mutation_tree_figure.py
$P docs/calibration/crm/jh4_extras.py
$P docs/calibration/crm/targeting_jh4.py
echo CHAIN_A_DONE
