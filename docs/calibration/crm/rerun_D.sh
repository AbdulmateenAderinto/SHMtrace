#!/bin/zsh
cd ~/Nanoshm/nanoshm_v0.2
.venv/bin/python docs/calibration/crm/public_illumina.py ch12 hmces ldha
.venv/bin/python docs/calibration/crm/targeting_ch12.py
echo CHAIN_D_DONE
