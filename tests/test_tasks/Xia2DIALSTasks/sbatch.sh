#!/bin/bash
#SBATCH --exclusive
#SBATCH -p fujitsu
#SBATCH -t 01-00:00
#SBATCH -N3
#SBATCH --mem=0

source /home/aarfin/start_edna2.sh
cd /gpfs/offline1/staff/biomax/aarfin/edna2/tests/test_tasks/Xia2DIALSTasks
python /gpfs/offline1/staff/biomax/aarfin/edna2/tests/test_tasks/Xia2DIALSTasks/Xia2DialsTaskTest_problems.py
#python /gpfs/offline1/staff/biomax/aarfin/edna2/tests/test_tasks/Xia2DIALSTasks/Xia2DialsTaskTest.py

