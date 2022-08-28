#!/bin/bash
#PBS -r n
#PBS -N danpos
#PBS -q default
#PBS -m e
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2000mb

module load python/3.5.1
module load R/3.5.1

cd /data/simulation/bin
Rscript caculate_distance.R -i /data -a 1 -b 100 -o /data/window100.txt -w 200
