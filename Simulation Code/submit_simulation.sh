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

cd /archive/tmhgxw22/Hi-C_strip/data/simulation/bin
Rscript simulation_stripe.R -a 51 -b 100 -o /archive/tmhgxw22/Hi-C_strip/data/simulation/data
