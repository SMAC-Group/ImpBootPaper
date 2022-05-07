#!/bin/sh

#SBATCH --job-name=recombin
#SBATCH --partition=debug-cpu
#SBATCH --time=0-00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2

source ~/jime/setting.sh

srun R CMD BATCH --no-save R/${MODEL}_recomb.R
