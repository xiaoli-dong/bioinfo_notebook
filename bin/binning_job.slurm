#!/bin/bash
# ---------------------------------------------------------------------
# An example SLURM script for running a job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --job-name=binning
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
mkdir concoct maxbin metabat dastool

cd metabat
module load miniconda3/metabat2
metabat -i ../../euk_contigs.fa -a ../../mapping/depth.txt -o metabat_out -t $SLURM_CPUS_PER_TASK -v
grep ">" *.fa | perl -ne 'if(/(^\S+?)\.fa:>(\S+)/){print "$2\t$1\n";}' > metabat_contigs2bin.tsv
cp metabat_contigs2bin.tsv ../dastool
module unload miniconda3/metabat2


################## maxbin2 binning
cd ../maxbin
cp ../../mapping/depth.txt .
perl ~/bin/depth2abund.pl depth.txt
mkdir maxbin_output
module load miniconda3/maxbin2

run_MaxBin.pl -contig ../../euk_contigs.fa -out maxbin_output/maxbin -abund_list abund_list.txt -thread 8
cd maxbin_output
grep ">" *.fasta | perl -ne 'if(/(^\S+?)\.fasta:>(\S+)/){print "$2\t$1\n";}' > maxbin_contigs2bin.tsv
cp maxbin_contigs2bin.tsv ../../dastool
module unload miniconda3/maxbin2

################## concoct binning
cd ../../concoct
module load miniconda3/concoct
cut_up_fasta.py ../../euk_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed ../../mapping/merged.sorted.bam > coverage_table.tsv
mkdir concoct_output
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/concoct -t $SLURM_CPUS_PER_TASK
cd concoct_output
merge_cutup_clustering.py concoct_clustering_gt1000.csv > clustering_merged.csv
mkdir fasta_bins
extract_fasta_bins.py ../../../euk_contigs.fa clustering_merged.csv --output_path fasta_bins
cd fasta_bins
grep ">" *.fa | perl -ne 'if(/(^\S+?)\.fa:>(\S+)/){print "$2\t$1\n";}' > concoct_contigs2bin.tsv
cp concoct_contigs2bin.tsv ../../../dastool
module unload miniconda3/concoct

################# run dastool
cd ../../../dastool
module load miniconda3/das_tool
DAS_Tool -i concoct_contigs2bin.tsv,maxbin_contigs2bin.tsv,metabat_contigs2bin.tsv -l concoct,maxbin,metabat -c ../../euk_contigs.fa -o br --search_engine diamond  -t $SLURM_CPUS_PER_TASK --write_bin_evals 1 --create_plots 1 --write_bins 1 > DAS_Tool.log.txt
module unload miniconda3/das_tool
# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
