#!/bin/bash
#SBATCH --nodes 1
#SBATCH -t 06:00:00
#SBATCH -o log.out
#SBATCH --partition skylake-gold

source setup.sh
for S in {0..27}
do
    for R in {0..15}
    do
        ./LowAlttiudeSim -s ${S} -r ${R} >> output1_cos2.txt
    done
done
