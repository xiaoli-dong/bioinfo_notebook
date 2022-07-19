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
module load java

./scientific_code input1.dat input2.dat input3.dat 

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
