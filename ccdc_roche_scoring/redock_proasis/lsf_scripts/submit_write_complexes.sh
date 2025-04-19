#!/bin/bash
bsub <redock_proasis/lsf_write_complexes_gold.sh
bsub <redock_proasis/lsf_write_complexes_vina.sh
bsub <redock_proasis/lsf_write_complexes_reference.sh
