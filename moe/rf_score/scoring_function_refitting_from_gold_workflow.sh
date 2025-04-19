#!/bin/bash
echo "running Scoring workflow.."

conda activate ccdc_roche_scoring3.3

echo "Calculating RF descriptors..."
target=$1

file_count=$(ls -d docking_job_* | wc -l)
out=$(bsub -J rf_scoring[1-$file_count] -q preempt -o docking_job_%I/lsf_%J_score.out -W 01:00 "/ccdc_roche_devel/ccdc_roche_scoring/moe/rf_score/lsf_input_scoring_function_refitting.sh" $target)
echo $out
echo "waiting for all HPC jobs to finish..."
