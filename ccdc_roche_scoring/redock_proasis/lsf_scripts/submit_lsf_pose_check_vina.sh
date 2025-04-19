#!/bin/bash
# file_count=$(find ~/scratch/LoS/protonate3d_2024-12-03/redocking/ -type d | wc -l)
# echo $file_count
# 190257
# short
bsub -J pose_check[1-99999] -n 1 -M 1G -o lsf_%J_pose_check.out -W 01:00 -q short "redock_proasis/lsf_pose_check_vina.sh"
bsub -J pose_check[100000-190257] -n 1 -M 1G -o lsf_%J_pose_check.out -W 01:00 -q short "redock_proasis/lsf_pose_check_vina.sh"
