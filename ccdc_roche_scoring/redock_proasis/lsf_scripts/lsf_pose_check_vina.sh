#!/bin/bash



conda activate ccdc_roche_scoring3.3

dirs=(~/scratch/LoS/protonate3d_2024-12-03/redocking/*)
dir=${dirs[$((LSB_JOBINDEX - 1))]}
echo $dir
cd $dir
pose-check -p protein.mol2 -l vina_results.sdf -o pose_check_vina.sdf
pose-check -p protein.mol2 -l ligand_native.sdf -o pose_check_native.sdf

#find /home/tosstora/scratch/LoS/protonate3d_2024-12-03/redocking/ -maxdepth 2 -type f -name 'pose_check_vina.sdf' -print0 | xargs -0 cat > /LoS/molecular_recognition_vina_2024-12-03/pose_check_vina_cat.sdf
