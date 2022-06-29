#!/bin/bash
#SBATCH --nodes 1
#SBATCH -t 06:00:00
#SBATCH -o log.out
#SBATCH --partition skylake-gold

source setup.sh
for S in {31..49}
do
    ./AltitudeDensitySim -s ${S} -r 20 >> output2_alt.txt
done
