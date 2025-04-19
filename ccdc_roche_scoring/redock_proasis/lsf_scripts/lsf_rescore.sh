#!/bin/bash



conda activate ccdc_roche_scoring3.3


folders=(./*)
folder=${folders[$((LSB_JOBINDEX - 1))]}
echo $folder
cd $folder
rescore
