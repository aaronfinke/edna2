#!/bin/bash
#SBATCH --job-name="EDNA2_SlurmTestTask"
#SBATCH --partition=all
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/offline1/staff/biomax/aarfin/edna2/tests/test_tasks/SlurmTasks/SlurmTestTask_gaawybp7
#SBATCH --output=EDNA2_SlurmTestTask_%j.out
#SBATCH --error=EDNA2_SlurmTestTask_%j.err
echo $HOSTNAME 
sleep 60
