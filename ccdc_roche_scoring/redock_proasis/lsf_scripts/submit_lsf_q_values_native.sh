#!/bin/bash
# file_count=$(find . -type f -name "ligand_native.sdf" | wc -l)
# echo $file_count
# 190255
# preempt
#bsub -J q_values[1-$file_count] -n 1 -M 1G -o lsf_%J_q_values.out -W 01:00 -q preempt -app snakemake "./lsf_q_values.sh"

# short
bsub -J q_values[1-99999] -n 1 -M 1G -o lsf_%J_q_values.out -W 01:00 -q short "redock_proasis/lsf_q_values.sh"
bsub -J q_values[99999-190255] -n 1 -M 1G -o lsf_%J_q_values.out -W 01:00 -q short "redock_proasis/lsf_q_values.sh"
