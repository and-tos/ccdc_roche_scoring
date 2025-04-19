export input_ligand_file=$1


conda activate ccdc_roche_scoring3.3

ccdc_roche_devel=$(dirname $(dirname $(dirname $(dirname -- "$0"))))

#ligands
#minimize in rdkit to ensure stereochemistry perception in CSD API and OE

ccdc_roche_3.9/bin/python ${ccdc_roche_devel}/ccdc_roche_scoring/getreasonabletautomers.py $input_ligand_file input_ligands/input_ligands_rdkit_quacpac.sdf
/apps/rocs/restricted/pRED/2020.08/cascadelake/software/moka/3.2.2/blabber_sd --pH=7.4 -t 20 -e -o input_ligands/input_ligands_rdkit_quacpac_moka_1.sdf --load-model=/apps/rocs/restricted/pRED/2020.08/cascadelake/software/moka/3.2.2/active_model.mkd input_ligands/input_ligands_rdkit_quacpac.sdf
ccdc_roche_3.9/bin/python ${ccdc_roche_devel}/ccdc_roche/python/database_utils/split_sdf.py -i input_ligands/input_ligands_rdkit_quacpac_moka_1.sdf -o input_ligands/input_ligands_rdkit_quacpac_moka_split -N 1 -m --nproc $LSB_DJOB_NUMPROC
