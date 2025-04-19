#!/bin/bash

#BSUB -J refitting
#BSUB -n 1
#BSUB -W 12:00
#BSUB -q long
#BSUB -M 8G
#BSUB -o lsf_refitting_%J.out

python /ccdc_roche_devel/ccdc_roche_scoring/ccdc_roche_scoring/stat_potential.py --target magl --task taskname
