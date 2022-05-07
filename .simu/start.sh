#!/bin/sh

ID=$(sbatch --parsable --array=1-1000 launch.sh)
ID=$(sbatch --parsable --dependency=afterany:${ID} recomb.sh)
sbatch --dependency=after:${ID} clean.sh
