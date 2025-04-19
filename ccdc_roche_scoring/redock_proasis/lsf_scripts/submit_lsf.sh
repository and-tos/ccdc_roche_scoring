#!/bin/bash
# file_count=99999 #431564
# echo $file_count

# short
# mkdir lsf_out

# bsub -J redock_proasis[1-20811] -n 1 -M 16G -o lsf_out/lsf_%J_%I_redock_proasis.out -W 01:00 -q short "redock_proasis/lsf_redock_proasis.sh"

echo "Submit batch 1"
bsub -J redock_proasis[1-99999] -n 1 -M 16G -o lsf_out/lsf_%J_%I_redock_proasis.out -W 01:00 -q short "redock_proasis/lsf_redock_proasis.sh"

echo "Submit batch 2"
bsub -J redock_proasis[99999-199997] -n 1 -M 16G -o lsf_out/lsf_%J_%I_redock_proasis.out -W 01:00 -q short "redock_proasis/lsf_redock_proasis.sh"

echo "Submit batch 3"
bsub -J redock_proasis[199997-299995] -n 1 -M 16G -o lsf_out/lsf_%J_%I_redock_proasis.out -W 01:00 -q short "redock_proasis/lsf_redock_proasis.sh"

echo "Submit batch 4"
bsub -J redock_proasis[299995-399993] -n 1 -M 16G -o lsf_out/lsf_%J_%I_redock_proasis.out -W 01:00 -q short "redock_proasis/lsf_redock_proasis.sh"

echo "Submit batch 5"
bsub -J redock_proasis[399993-431564] -n 1 -M 16G -o lsf_out/lsf_%J_%I_redock_proasis.out -W 01:00 -q short "redock_proasis/lsf_redock_proasis.sh"
