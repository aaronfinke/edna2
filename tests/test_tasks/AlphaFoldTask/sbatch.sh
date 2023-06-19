#!/bin/bash
#SBATCH --job-name=AF_7QRZ
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive

source /data/staff/biomax/domfas/start_sbatch.sh
cd  /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask
python /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask/AlphaFold_exec_test.py
