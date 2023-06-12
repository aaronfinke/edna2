#!/bin/bash
#SBATCH -p fujitsu
#SBATCH --exclusive
#SBATCH -t 02:00:00


source /home/aarfin/start_edna2.sh
cd  /gpfs/offline1/staff/biomax/aarfin/edna2/tests/test_tasks/AutoPROCTask
python /gpfs/offline1/staff/biomax/aarfin/edna2/tests/test_tasks/AutoPROCTask/AutoPROCTask_exec_test.py
