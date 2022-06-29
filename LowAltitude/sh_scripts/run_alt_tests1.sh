#!/bin/bash
#SBATCH --nodes 1
#SBATCH -t 06:00:00
#SBATCH -o log.out
#SBATCH --partition skylake-gold

source setup.sh
for S in {0..30}
do
    ./AltitudeDensitySim -s ${S} -r 20 >> output1_alt.txt
done
