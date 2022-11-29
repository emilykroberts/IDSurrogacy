#!/bin/bash
#SBATCH --job-name=IDsims

#SBATCH --time=0-5:50:00
#SBATCH --mail-user=ekrobe@umich.edu
#SBATCH --mail-type=END
#SBATCH --array=1-1600
#SBATCH -o /dev/null

R CMD BATCH --no-save --no-restore IDsims.R $SLURM_ARRAY_TASK_ID
