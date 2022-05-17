#!/bin/bash
#SBATCH --nodes 1
#SBATCH -t 06:00:00
#SBATCH -o log.out
#SBATCH --partition skylake-gold

source setup.sh
for S in {0..27}
do
    for R in {0..29}
    do
        ./LowAlttiudeSim -s ${S} -r ${R} -uni >> output.txt
    done
done
