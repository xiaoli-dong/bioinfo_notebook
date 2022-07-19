#!/bin/bash
# ---------------------------------------------------------------------
# An example SLURM script for running a job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --job-name=SLURM_example

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00

#SBATCH --mem=100mb

#SBATCH --partition=parallel

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
module load megahit

megahit -1 test_data/r1.il.fa.gz -2 test_data/r2.il.fa.bz2 -o out 

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
