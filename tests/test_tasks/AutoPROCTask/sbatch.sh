#!/bin/bash
#SBATCH --exclusive
#SBATCH -t 02:00:00


source /gpfs/offline1/staff/biomax/aarfin/edna2_playground/start.sh
cd  /gpfs/offline1/staff/biomax/aarfin/edna2_maxiv/edna2/tasks/test/AutoPROCTask
python /gpfs/offline1/staff/biomax/aarfin/edna2_maxiv/edna2/tasks/test/AutoPROCTask/AutoPROCTask_exec_test.py
