#!/bin/bash
#
# $1: target
# $2: 2

conda activate ccdc_roche_scoring3.3

target=$1
flexible_residues=$2
echo $target
echo $flexible_residues
i=$LSB_JOBINDEX

ccdc_roche_devel=$(dirname $(dirname $(dirname $(dirname -- "$0"))))

cd docking_job_$i
python ${ccdc_roche_devel}/ccdc_roche_scoring/ccdc_roche_scoring/docking.py --input_ligands input_ligands/input_ligands_rdkit_quacpac_moka_split_$i.sdf -t $target -fr=$flexible_residues
best_soln=$(ls ./*/best_soln*.sdf)
gold_conf=$(ls ./*/*/gold.conf)
echo $best_soln
pose-check --ligands $best_soln --gold_conf $gold_conf --output $best_soln

if [ "$target" != "default" ]; then
    python ${ccdc_roche_devel}/ccdc_roche_scoring/ccdc_roche_scoring/stat_potential_inference.py -t $target -g --task magl_h_ms_pic50
fi
