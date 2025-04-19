#!/bin/bash


conda activate ccdc_roche_scoring3.3

target=$1
echo $target
i=$LSB_JOBINDEX
cd docking_job_${i}
python /ccdc_roche_devel/ccdc_roche_scoring/ccdc_roche_scoring/descriptor_dfs.py -t $target --docking_dir .
