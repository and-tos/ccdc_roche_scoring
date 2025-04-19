#!/bin/bash

conda activate ccdc_roche_scoring3.3

sdf=$1
out=$2
gmean=$3

ccdc_roche_devel=$(dirname $(dirname $(dirname $(dirname -- "$0"))))

q-values -isdf $sdf -o $out $gmean
