#!/bin/bash
file_count=$(find . -type f -name "api_gold.conf*" | wc -l)
echo $file_count

# preempt
#bsub -J pose_check[1-$file_count] -n 1 -M 1G -o lsf_%J_pose_check.out -W 01:00 -q preempt -app snakemake "./lsf_pose_check.sh"

# short
bsub -J rescore[1-$file_count] -n 1 -M 1G -o lsf_%J_rescore.out -W 01:00 -q short "redock_proasis/lsf_rescore.sh"
