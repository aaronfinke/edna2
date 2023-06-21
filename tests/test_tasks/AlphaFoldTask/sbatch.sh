#!/bin/bash
#SBATCH --job-name=AF_7QRZ
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive

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

cd  /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask
python /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask/AlphaFold_exec_test.py
