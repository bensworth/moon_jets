#!/bin/bash
rm hostfile_${SLURM_JOB_ID}.txt
echo $(head -n 1 $PBS_NODEFILE) slots=2 >> hostfile_${SLURM_JOB_ID}.txt
L=$(cat $PBS_NODEFILE | wc -l)
for i in $(seq 2 $L); do
    echo $(sed -n ${i}p $PBS_NODEFILE) slots=1 >> hostfile_${SLURM_JOB_ID}.txt
done