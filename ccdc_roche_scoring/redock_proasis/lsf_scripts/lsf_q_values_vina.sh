#!/bin/bash

conda activate ccdc_roche_scoring3.3

files=(./*/vina_results.sdf)
conf=${files[$((LSB_JOBINDEX - 1))]}
echo $conf
outfile=$(dirname "$conf")/q_values_vina.csv
q-values -isdf $conf -o $outfile
