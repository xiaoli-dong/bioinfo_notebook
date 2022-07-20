#!/bin/bash
# ---------------------------------------------------------------------
# An example SLURM script for running a job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --job-name=eukrep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --partition=parallel
#SBATCH --mail-user=xdong@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
module load miniconda3/eukrep

EukRep -i final_contgs.fasta -o euk_contigs.fa --min 1000 --prokarya prok_contigs.fasta

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
