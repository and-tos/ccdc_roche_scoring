#!/bin/bash



conda activate ccdc_roche_scoring3.3


files=(./*/ligand_native.sdf)
conf=${files[$((LSB_JOBINDEX - 1))]}
echo $conf
outfile=$(dirname "$conf")/q_values_native.csv
q-values -isdf $conf -o $outfile
