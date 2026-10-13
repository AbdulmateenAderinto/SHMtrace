#!/bin/zsh
cd ~/Nanoshm/nanoshm_v0.2
until grep -q CHAIN_A_DONE docs/calibration/crm/rerun_A.log 2>/dev/null; do sleep 30; done
P=.venv/bin/python
set -x
$P docs/calibration/crm/spikein.py jh4 n114 su_lps pim1 --profiles recommended --no-lofreq --tag _rec
$P docs/calibration/crm/merge_spikein.py jh4 n114 su_lps pim1
rm -rf docs/calibration/crm/out/lineage_sim_jh4
.venv/bin/shmtrace-simulate-lineages --region jh4_intron --control "/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ/DHJGKB_5_AID_KO_1.fastq" --control "/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ/DHJGKB_6_AID_KO_2.fastq" --control "/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ/DHJGKB_7_AID_KO_3.fastq" --reference src/shmtrace/references/mouse_jh4_intron.fasta --contaminant src/shmtrace/references/mouse_jh4_contaminant.fasta --targeting-model ~/Nanoshm/models/HH_S5F.csv --out docs/calibration/crm/out/lineage_sim_jh4
$P docs/calibration/crm/single_read_spike.py 0,25,50,100 3
cd docs/calibration && ../../.venv/bin/python quant_accuracy.py --amplicons JH4 --jobs 3 --out crm/out/quant_accuracy_jh4 && ../../.venv/bin/python quant_accuracy.py --amplicons N114 --jobs 3 --out crm/out/quant_accuracy_n114
echo CHAIN_C_DONE
