#!/bin/bash

# Needs to be run in a directory with the following structure
#
# template_docking_dir
# |-tmp_aligned_3d_sdf_sanitized
# | |-native_ligands.sdf (use corresponding protein ID in template ID 1edeg.A)
# |-input_ligands
# |  |-input_ligands_rdkit_quacpac_moka_split_*.sdf
# |-docking
# |-tmp_aligned_for_MOE_sanitized
# | |-1edeg_dry.pdb
# |-apo_proteins
# | |-1edeg_dry.mol2

#
# $1: Protonate ligands [0, 1]
# $2: Enumerate ligand stereocenters [0, 1]
# $3: target name [default, egfr, magl, pde-10]
# $4: verbose output [0,1], currently irrelevant
# $5: flexible residue name [ASP123]
#
# usage examples:
# docking_workflow.sh 1 default

echo "running docking workflow.."

ccdc_roche_devel=$(dirname $(dirname $(dirname $(dirname -- "$0"))))

# Enumerate stereocenters and protonate
if [ $1 -eq 1 ] && [ $2 -eq 1 ]; then
  echo "preparing input ligands..."
  out=$(bsub -J ligprep -q preempt -o input_ligands/lsf_%J_ligprep.out -W 01:00 -n 24 "${ccdc_roche_devel}/ccdc_roche_scoring/moe/template_docking/ligand_proton_and_stereo_enum.sh" input_ligands.sdf)
  echo $out
  bwait -w "ended(${out//[!0-9]/})"
fi

# Only protonate
if [ $1 -eq 1 ] && [ $2 -eq 0 ]; then
  echo "preparing input ligands..."
  out=$(bsub -J ligprep -q preempt -o input_ligands/lsf_%J_ligprep.out -W 01:00 -n 24 "${ccdc_roche_devel}/ccdc_roche_scoring/moe/template_docking/ligand_protonation.sh" input_ligands.sdf)
  echo $out
  bwait -w "ended(${out//[!0-9]/})"
fi

# Only enumerate stereocenters
if [ $1 -eq 0 ] && [ $2 -eq 1 ]; then
  echo "preparing input ligands..."
  out=$(bsub -J ligprep -q preempt -o input_ligands/lsf_%J_ligprep.out -W 01:00 -n 24 "${ccdc_roche_devel}/ccdc_roche_scoring/moe/template_docking/ligand_stereo_enumeration.sh" input_ligands.sdf)
  echo $out
  bwait -w "ended(${out//[!0-9]/})"
fi

# No enumeration and no protonation
if [ $1 -eq 0 ] && [ $2 -eq 0 ]; then
  echo "Not preparing input ligands."
  obabel input_ligands.sdf -osdf -O input_ligands/input_ligands_rdkit_quacpac_moka_split_.sdf -m
fi

if [ $3 == "default" ]; then
  mkdir tmp_aligned_3d_sdf_sanitized
  for file in ligand_*.sdf; do
    cat "${file}" >>tmp_aligned_3d_sdf_sanitized/native_ligands.sdf
    printf "\n\$\$\$\$\n\n" >>tmp_aligned_3d_sdf_sanitized/native_ligands.sdf
  done
  rm ligand_*.sdf
fi

# docking

cd docking

target=$3
flexible_residues=$5
echo $target
echo $flexible_residues

file_count=$(ls ../input_ligands | grep input_ligands_rdkit_quacpac_moka_split | wc -l)

for i in $(seq 1 $file_count); do
  mkdir docking_job_$i
done

out=$(bsub -J template_dock[1-$file_count] -q short -o docking_job_%I/lsf_%J_dock.out -W 02:00 -app snakemake "${ccdc_roche_devel}/ccdc_roche_scoring/moe/template_docking/lsf_input_gold_docking.sh" $target $flexible_residues)
echo "bsub -J template_dock[1-$file_count] -q short -o docking_job_%I/lsf_%J_dock.out -W 02:00 -app snakemake "${ccdc_roche_devel}/ccdc_roche_scoring/moe/template_docking/lsf_input_gold_docking.sh" $target $flexible_residues"
echo $out
echo "waiting for all HPC jobs to finish..."
bwait -w "ended(${out//[!0-9]/})"

echo "Collecting results..."

conda activate ccdc_roche_scoring3.3

python ${ccdc_roche_devel}/ccdc_roche_scoring/ccdc_roche_scoring/join_docking_results.py
chmod +rwx best_docking_solutions.sdf

echo "cleaning up..."
# rm -r docking_job_*
# rm -r ../input_ligands
