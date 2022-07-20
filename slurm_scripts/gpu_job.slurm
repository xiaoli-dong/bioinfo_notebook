#!/bin/bash
# ---------------------------------------------------------------------
# An example SLURM script for running a job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --job-name=SLURM_example

#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:2
#SBATCH --mem=40GB
#SBATCH --time=8:00:00
#SBATCH --partition=gpu-v100 

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------

module load guppy/3.6.1
guppy_basecaller -i fast5_pass -s fastq_pass --num_callers 4 --config dna_r9.4.1_450bps_hac.cfg --device cuda:all:100% --qscore_filtering --min_qscore 7



# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
