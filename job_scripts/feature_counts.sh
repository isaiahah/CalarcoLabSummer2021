#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=03:00:00
#SBATCH --job-name=feature_counts
#SBATCH --output=output.out
#SBATCH --error=output.err

cd $SLURM_SUBMIT_DIR

module load NiaEnv/2019b
module load gcc/9.2.0
module load r/4.0.3

Rscript feature_counts.R
