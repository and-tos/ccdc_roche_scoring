#!/bin/bash
echo "running Scoring workflow.."
conda deactivate
conda activate ccdc_roche_3.9

echo "Scoring..."
target=$1

file_count=$(ls -d scoring_job_* | wc -l)
out=$(bsub -J rf_scoring[1-$file_count] -q preempt -o lsf_%J_scoring.out -W 01:00 "/ccdc_roche_devel/ccdc_roche_scoring/moe/rf_score/lsf_input_scoring.sh" $target)
echo "waiting for all HPC jobs to finish..."
bwait -w "ended(${out//[!0-9]/})"

cat scoring_job_*/rescored_ligands.sdf >>rescored_ligands.sdf
