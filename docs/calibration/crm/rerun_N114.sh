#!/bin/zsh
# N114 with empty-vector subclone C3's second library (C3_Empty_March_2024_2) in place of the first
cd ~/Nanoshm/nanoshm_v0.2
P=.venv/bin/python
set -x
$P docs/calibration/n114_shmtrace.py
$P docs/calibration/depth_shmtrace.py n114
$P docs/calibration/mutation_tree_figure.py
$P docs/calibration/crm/jh4_extras.py
$P docs/calibration/crm/benchmark_real.py n114
$P docs/calibration/crm/spikein.py n114 --profiles recommended,geneious,recommended_bh
(cd docs/calibration && ../../.venv/bin/python quant_accuracy.py --amplicons N114 --jobs 3 --out crm/out/quant_accuracy_n114)
echo CHAIN_N114_DONE
