#!/bin/bash

conda activate ccdc_roche_scoring3.3

sdf=$1
out=$2

ccdc_roche_devel=$(dirname $(dirname $(dirname $(dirname -- "$0"))))

python ${ccdc_roche_devel}/ccdc_roche_scoring/ccdc_roche_scoring/pose_processor.py -isdf $sdf -o $out
