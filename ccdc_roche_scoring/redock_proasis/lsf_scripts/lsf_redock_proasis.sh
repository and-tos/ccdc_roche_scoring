#!/bin/bash

conda activate ccdc_roche_scoring3.3

entry_index=$((${LSB_JOBINDEX} - 1))
echo $entry_index
redock-proasis -i ${entry_index} --proasis_csdsql /LoS/protonate3d_2024-12-03/full_p2cq_pub_2024-12-03.csdsql
