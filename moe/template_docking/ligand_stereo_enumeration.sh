#!/bin/bash

export input_ligand_file=$1

conda activate ccdc_roche_scoring3.3

ccdc_roche_devel=$(dirname $(dirname $(dirname $(dirname -- "$0"))))

#ligands
#minimize in rdkit to ensure stereochemistry perception in CSD API and OE

ccdc_roche_3.9/bin/python ${ccdc_roche_devel}/ccdc_roche/python/database_utils/split_sdf.py -i $input_ligand_file -o input_ligands/input_ligands_rdkit_quacpac_moka_split -m -N 1 -es all --nproc $LSB_DJOB_NUMPROC
