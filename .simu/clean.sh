#!/bin/sh

#SBATCH --job-name=clean
#SBATCH --partition=debug-cpu
#SBATCH --time=0-00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

rm *.Rout
rm slurm-*
