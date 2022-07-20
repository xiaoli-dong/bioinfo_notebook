# Eukaryotic genome and metagenome annotation
When the input is metagenome, we need to seperate eukaryotic contigs from the the prokaryotic contigs

## Separate eukaryotic contigs from prokaryotic contigs

```
#!/bin/bash
# ---------------------------------------------------------------------
# Scripts used to seperate eukaryotic contigs from prokaryotic contigs 
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

## Eukaryotic genome and metagenome annotation

After separations, I used the following procedures to annotate eukaryotic contigs. 

* Repeat identification and masking
* Gene prediction using Maker:
 * First round - gene pediction based on evidence
 * Second round - use ab-initio predictors to improve its predictions

When running script, please make sure to unload unused the modules which are installed as conda packages. Otherwise, there can be confliction happens. For example, different conda installed software may depend on different perl versions. When you load two programs which are using different perl versions. The one loaded earlier will be the perl you are using when you run the program. This may cause the second program to fail. 

```

#!/bin/bash
# ---------------------------------------------------------------------
# Scripts used to  annotate eukaryotic contigs
# ---------------------------------------------------------------------
#SBATCH --job-name=SLURM_example

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-01:00:00


#SBATCH --partition=parallel

# $1=input fasta file name, $2=base_path $3=busco dataset $4=august test training set
# sbatch my_job_maskrepeat.slurm input_fasta genomeid busco lineage
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


# ---------------------------------------------------------------------------                                                 
echo "Starting first round of Maker"
# --------------------------------------------------------------------------
# first round we use GeneMark-ES predicted model, Busco trained Augustus model
# and the swiss-prot db to do the gene structure prediction

mkdir 1stRoundTrain
cd 1stRoundTrain

# --------------------------------------------------------------------------- 
echo "First round of maker -- train GeneMark-ES"
# --------------------------------------------------------------------------    

#train GeneMark-ES
module load GeneMark/GeneMark-ES/v4
mkdir genemark_es
cd genemark_es
gmes_petap.pl --ES -min_contig 5000 --cores $SLURM_CPUS_PER_TASK --sequence ../../genome.fasta.masked
module rm GeneMark/GeneMark-ES/v4
#move out genemark_es directory
cd ../
module rm GeneMark/GeneMark-ES/v4
# ---------------------------------------------------------------------                                                          
echo "First round of maker -- train Augustus using busco"
# ---------------------------------------------------------------------    
#train augustus with busco
module load miniconda3/busco/4.1.4

mkdir augustus
cd augustus
output="1rnd_${gid}"
#busco lineage to be used
busco -i ../../genome.fasta.masked -o ${output} -l ${lineage}  -m genome -c $SLURM_CPUS_PER_TASK --long
echo "cp -r ${output}/run_${lineage}/augustus_output/retraining_parameters/BUSCO_${output} /work/ebg_lab/software/augustus_config/config/species/" 
cp -r ${output}/run_${lineage}/augustus_output/retraining_parameters/BUSCO_${output} /work/ebg_lab/software/augustus_config/config/species/

#move out augustus directory
cd ../
# ---------------------------------------------------------------------                               
echo "First round of maker -- start"
# ---------------------------------------------------------------------                              
module load miniconda3/maker/2.31.10
mkdir maker
cd maker
cp ../../../maker_config/maker1/* .
sed -i "s/EBGGID/${gid}/" maker_opts.ctl
#copy the maker configuration files over
maker -c $SLURM_CPUS_PER_TASK
cd *.maker.output
gff3_merge -d *_master_datastore_index.log
fasta_merge -d *_master_datastore_index.log -o genome

#look for single copy orthologous genes within the bin and estimate completeness
module load miniconda3/busco/4.1.4
busco -i *.maker.proteins.fasta -l ${lineage} -o busco_eval -m prot -f

#move out of the maker directory
cd ../../

#move out of the 1stRoundTrain directory
cd ../

# ---------------------------------------------------------------------------                                                    
echo "Starting second round of Maker"
# --------------------------------------------------------------------------                                                      
# first round we use GeneMark-ES predicted model, Busco trained Augustus model                                                   # and the swiss-prot db to do the gene structure prediction                                                                      

mkdir 2ndRoundTrain
cd 2ndRoundTrain

# ---------------------------------------------------------------------------                                                    
echo "Second round of maker -- train snap"
echo "Training snap using maker model generated from the 1stRoundTrain"
# --------------------------------------------------------------------------

module load snap/snap
mkdir snap
cd snap

# export 'confident' gene models from MAKER:
#select model with aed <= 0.5 && with >= 50 aa 
maker2zff  -c 0 -e 0 -x 0.25 -l 50 ../../1stRoundTrain/maker/genome.fasta.maker.output/genome.fasta.all.gff 

# gather some stats and validate
fathom genome.ann genome.25.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1

# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1

# move out of the params directory
cd ../

# use snap created expot.ann and export.dna file as input to create augustus input
# the following script needs to be run in the director containing export.ann and export.dna
# the output can be used to train augustus later on

/work/ebg_lab/software/scripts/zff2augustus_gbk.pl > augustus.gbk

# Training SNAP, assembly the HMM
hmm-assembler.pl genome params > genome.length50_aed0.5.hmm

#come out of the snap directory
cd ..


# ---------------------------------------------------------------------------                                                    
echo "Second round of maker -- train augustus using the 1st maker model"
# --------------------------------------------------------------------------
module load miniconda3/augustus/3.3.3                                                 
mkdir augustus
cd augustus
ln -s ../snap/augustus.gbk

# Like in many machine learning approaches, we will split the now created augustus.gbk
# file into a training and a test set: This generates a file *.gbk.test with xx randomly
# chosen loci and a disjoint file *.gbk.train with the rest of the loci from genes.gb:

randomSplit.pl augustus.gbk ${testset}

# We now have to create a new species for our AUGUSTUS training. The parameter files will
#be write to AUGUSTUS_CONFIG_PATH species directory
output="2rnd_${gid}"
new_species.pl --species=${output}

# Lets now train AUGUSTUS with the training set file, evaluate the training and save the 
# output of this test to a file for later reference:

etraining --species=${output} augustus.gbk.train
augustus --species=${output} augustus.gbk.test | tee first_training.out

# Once this is done, we are ready to improve prediction parameters of the models using 
# the optimize_augustus.pl script again located in the AUGUSTUS/scripts directory: 
# this step is very slow. it can take days to run

optimize_augustus.pl --cpus=$SLURM_CPUS_PER_TASK --kfold=$SLURM_CPUS_PER_TASK  --species=${output}  augustus.gbk.train 

# Retrain and test AUGUSTUS with the optimized paramters and compare the results to the first run:

etraining --species=${output} augustus.gbk.train
augustus --species=${output} augustus.gbk.test | tee second_training.out

#come out of the augustus directory
cd ..

mkdir maker
cd maker

#copy maker files over to the directory
cp ../../../maker_config/maker2/* .
sed -i "s/EBGGID/${gid}/" maker_opts.ctl

#run maker with the 1st round of maker model trained SNAP and augustus models
maker -c $SLURM_CPUS_PER_TASK
cd *.maker.output
gff3_merge -d *_master_datastore_index.log
fasta_merge -d *_master_datastore_index.log -o genome

#look for single copy orthologous genes within the bin and estimate completeness
module load miniconda3/busco/4.1.4
busco -i *.maker.proteins.fasta -l ${lineage} -o busco_eval -m prot -f

#come out the maker directory
cd ../..

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------

```

# Build an eukaryotic phylogenetic tree

I was building phylogenetic tree mainly referred to the following two papers:

  *	[A new view of the tree of life](https://www.nature.com/articles/nmicrobiol201648)
  *	[Genome-reconstruction for eukaryotes from complex natural microbial communities](https://genome.cshlp.org/content/early/2018/03/01/gr.228429.117)

In the papers, they used predicted, aligned, and concatenated 16 Ribosomal protein (RPs) genes to build the trees. 

The following are the procedures I used to build the algae tree years ago:

- Downloaded 7 genomes from ncbi. Those genomes belong to Chlorophyta they were not included in the laura hug's reference tree. 
- Extracted RPs genes from eukaryotic contig bins and the ncbi genomes using 16 euk rps hmm models

- Concatenated RPs genes in order for each genomes (in total 8 gene set)
- Extracted euk and arc alignments from laura hug's total alignments (euk+archaea+bac)
- Aligned concatenated RPs genes from my eukaryotic bin and ncbi downlaed genomes using the the extracted euk+archaea alignment from the preivous step as a reference:

```
mafft --addfragments  mycontigbin_ncbi.concat_rps.faa --keeplength --reorder --thread 8  16rps_align.archaea.eukaryota.fasta  > all.aligned.fasta
```
- build the tree using raxml
```
raxmlHPC-PTHREADS  -f a -s all.aligned.fasta -n whole -m PROTGAMMAAUTO -x 0123 -# 100 -p 012345  -T 20 -r maximun_likelihood_tree
```
After extracting RPs genes from contigs, you need manually inspect the alignment qualities. If some genes’ alignment qualities are poor, you should exclude them from the  further analysis

