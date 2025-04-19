#!/bin/bash
conda activate ccdc_roche_scoring3.3

target=$1
echo $target
i=$LSB_JOBINDEX
cd docking_job_${i}
ccdc_roche_3.9/bin/python /ccdc_roche_devel/ccdc_roche_scoring/stat_potential_inference.py -t $target -g --task taskname
