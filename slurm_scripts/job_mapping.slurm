#!/bin/bash
# ---------------------------------------------------------------------
# An example SLURM script for running a job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --job-name=mapping
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3-01:00:00
#SBATCH --partition=parallel
#SBATCH --mail-user=xdong@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
module load bwa/0.7.17
module load miniconda3/samtools
module load miniconda3/metabat2


echo "indexing contig files"

bwa index ../euk_contigs.fa
echo "mapping"
bwa mem -x ont2d -t $SLURM_CPUS_PER_TASK ../euk_contigs.fa /work/ebg_lab/eb/ciliates/nanopore/qc/ciliates_all.qc.fastq | samtools view -bS - | samtools sort - -o nanopore.sorted.bam
echo "samtools indexing sorted bam"
samtools index nanopore.sorted.bam

samtools merge merged.bam illumina.bbmap_sorted.bam nanopore.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth  depth.txt merged.bam


# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
