#!/bin/bash

conda deactivate
conda activate ccdc_roche_3.9

target=$1
echo $target
i=$LSB_JOBINDEX
cd scoring_job_$i
ccdc_roche_3.9/bin/python /ccdc_roche_devel/ccdc_roche_scoring/stat_potential_inference.py -t $target -d . -l input_ligands.sdf -p ../protein.pdb
