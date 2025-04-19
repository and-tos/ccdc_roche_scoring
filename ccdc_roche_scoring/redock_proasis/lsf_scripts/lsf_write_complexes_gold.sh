#!/bin/bash

#BSUB -J write_gold_complexes
#BSUB -n 32
#BSUB -W 24:00
#BSUB -q long
#BSUB -R "rusage[mem=256G/host]"
#BSUB -o lsf_write_complexes_%J.out

conda activate ccdc_roche_scoring3.3

python redock_proasis/write_docked_complexes.py -i gold -db roche
python redock_proasis/write_docked_complexes.py -i gold -db pub
