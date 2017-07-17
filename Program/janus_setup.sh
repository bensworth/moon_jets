#!/bin/bash

module load slurm/slurm
module load gcc/gcc-4.9.2
module load openmpi/openmpi-1.8.4_gcc-4.9.2
module load gsl/gsl-1.16_gcc-4.9.2

make -f parallel
