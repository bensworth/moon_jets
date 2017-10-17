#!/bin/bash
scontrol show hostnames >> hostname_${SLURM_JOB_ID}.txt
echo $(head -n 1 hostname_${SLURM_JOB_ID}.txt) slots=3 >> hostfile_${SLURM_JOB_ID}.txt
L=$SLURM_JOB_NUM_NODES
for i in $(seq 2 $L); do
    echo $(sed -n ${i}p hostname_${SLURM_JOB_ID}.txt) slots=2 >> hostfile_${SLURM_JOB_ID}.txt
done