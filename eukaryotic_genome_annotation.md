The first step is to separate eukaryotic contigs from prokaryotic contigs

```
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

```

After separation, we do the gene prediction with the following iterative approaches and the example command is included in the script. When running script, please make sure to unload unused the modules which are installed as conda packages. Otherwise, there can be confliction happens. For example, different conda installed software may depend on different perl versions. When you load two programs which are using different perl versions. The one loaded earlier will be the perl you are using when you run the program. This may cause the second program to fail. 
```
#!/bin/bash
# ---------------------------------------------------------------------
# AN EXAMPLESLURM script for running a job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --job-name=SLURM_example

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-01:00:00


#SBATCH --partition=parallel

# $1=input fasta file name, $2=base_path $3=busco dataset $4=august test training set

# example command: sbatch my_job_maskrepeat.slurm input_fasta genomeid busco lineage

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"

#input genome fasta file to be analyzed including path
fasta=$1
#uniqu genome identifier, can be binids...
gid=$2

#busco lineage to be used. for example "eukaryota_odb10"
lineage=$3

testset=$4

##example: sbatch euk_annot.slurm ../euk_contigs.fa euk eukaryota_odb10 100
#-----------------------------------------------------------------
echo "Starting to masking repeates with RepeatModeler and RepeatMasker"
# -----------------------------------------------------------------
module load repeatmasker/4.1.1
module load repeatmodeler/2.0.1

echo "mkdir ${gid}_out"
mkdir ${gid}_out
cd ${gid}_out
ln -s ../${fasta} genome.fasta 
mkdir RepeatModeler
cd RepeatModeler

#Build the database
echo "BuildDatabase -name genome_db -engine rmblast ../genome.fasta"
BuildDatabase -name genome_db -engine rmblast ../genome.fasta

# De novo model the repeats
echo "RepeatModeler -pa $SLURM_CPUS_PER_TASK -database genome_db"
RepeatModeler -pa $SLURM_CPUS_PER_TASK -database genome_db

echo "RepeatClassifier -consensi RM*/consensi.fa  -stockholm RM*/families.stk -engine rmblast"
RepeatClassifier -consensi RM*/consensi.fa  -stockholm RM*/families.stk -engine rmblast

# Mask the assembly
echo "RepeatMasker -lib RM*/consensi.fa.classified -no_is -xsmall -gff -pa $SLURM_CPUS_PER_TASK ../genome.fasta"
RepeatMasker -lib RM*/consensi.fa.classified -no_is -xsmall -gff -pa $SLURM_CPUS_PER_TASK ../genome.fasta
module rm repeatmasker/4.1.1
module rm repeatmodeler/2.0.1
#get out of repeatmodeler directory

cd ..

# ---------------------------------------------------------------------------
echo "End of masking repeats"
# --------------------------------------------------------------------------
```

Euk phylogenetic tree building
I was building phylogenetic tree mainly referred to the following two papers:
•	[A new view of the tree of life](https://www.nature.com/articles/nmicrobiol201648)
•	[Genome-reconstruction for eukaryotes from complex natural microbial communities](https://genome.cshlp.org/content/early/2018/03/01/gr.228429.117)

In the papers, they used predicted, aligned, and concatenated 16 Ribosomal protein (RPs) genes to build the trees. The downloaded hmm for euk and arc, and the scripts I wrote to fetch hmm from contigs are all located at: /gpfs/ebg_data/programs/euk_phylogenetic_tree_building

The following is the procedures I used to build the algae tree years ago:

- Downloaded 7 genomes from ncbi. Those genomes belong to Chlorophyta nut they were not included in the laura hug's reference tree. 
- Extract rps genes from sodalake bin and the ncbi genomes using 16 euk rps hmm modelsls 
- Concatenated 16 rps genes in order for each genomes (in total 8 gene set)
- Extracted euk + archaea alignments from laura hug's total alignments (euk+archaea+bac)
- Do the reference alignment using sodalake rps genes + ncbi genome rps genes as input  and the extracted euk+archaea alignment as references:
```
mafft --addfragments  sodaLake_ncbi.concat.faa --keeplength --reorder --thread 8  16rps_align.archaea.eukaryota.fasta  > all.aligned.fasta
```
- build the tree using raxml: 
```
bsub -o out.txt -e err.txt raxmlHPC-PTHREADS  -f a -s all.aligned.fasta -n whole -m PROTGAMMAAUTO -x 0123 -# 100 -p 012345  -T 20
RAxml selected GAMMA + VT likelihood model as the best-scoring AA model  
```
After extracted rps genes from contigs, you need to manually inspect the alignment qualities. If some the genes’ alignment qualities are poor, you should exclude them from the  further analysis

