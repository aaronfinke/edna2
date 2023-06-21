#!/bin/bash
module purge
module add fosscuda/2020b AlphaFold
export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b
EDNA2_PATH=/data/staff/biomax/domfas/edna2_alphafold
conda activate /home/domfas/miniconda3/envs/edna2       
export PATH=/home/domfas/miniconda3/envs/edna2/bin:/home/domfas/miniconda3/condabin:$PATH
export EDNA2_SITE=MAXIV_BIOMAX
export PATH=${EDNA2_PATH}/bin:$PATH
export PYTHONPATH=${EDNA2_PATH}/src