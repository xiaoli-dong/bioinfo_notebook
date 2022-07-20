
# Nanopore reads analysis
Here, I am sharing some of the step by step procedures I used in my past project procedures

### Create various plots for long read sequencing data
```
NanoPlot --plots dot -t 8 --fastq my_read.fastq -o my_outdir -p my_prefix
```
### Finding and remove adapters from nanopore reads
```
porechop -i my_read.fastq -o porechop.fastq.gz --threads 8
filtlong --min_length 1000 porechop.fastq.gz | gzip > qc.fastq.gz
```

### Nanopore read assembly
De novo assembly using flye

```
flye --nano-raw qc.fastq.gz --meta --genome-size 450m --out-dir assembly_fly -i 3 -t 12
```

### contig polish using nanopore long reads
Here, I used racon to do four rounds of the polish and then did consensus calling using medaka

```
#!/usr/bin/env bash

#BSUB -J polish_racon_medaka
#BSUB -q coursework
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -W 10:00
#BSUB -o polish_racon_medaka_%J.OUT
#BSUB -e polish_racib_medaka_%J.ERR

echo "Start: `date`"

################################ four rounds of long reads polish ##############################
bwa index ../assembly.fasta
bwa mem -t 20 -x ont2d ../assembly.fasta ../../qc/qc.fastq.gz -o assembly.mapping.sam
racon -t 20 ../../qc/qc.fastq.gz assembly.mapping.sam ../assembly.fasta > racon1.fasta

bwa index racon1.fasta
bwa mem -t 20 -x ont2d racon1.fasta ../../qc/qc.fastq.gz -o racon1.mapping.sam
racon -t 8 ../../qc/qc.fastq.gz racon1.mapping.sam racon1.fasta > racon2.fasta

bwa index racon2.fasta
bwa mem -t 20 -x ont2d racon2.fasta ../../qc/qc.fastq.gz -o racon2.mapping.sam
racon -t 20 ../../qc/qc.fastq.gz racon2.mapping.sam racon2.fasta > racon3.fasta

bwa index racon3.fasta
bwa mem -t 20 -x ont2d racon3.fasta ../../qc/qc.fastq.gz -o racon3.mapping.sam
racon -t 20 ../../qc/qc.fastq.gz racon3.mapping.sam racon3.fasta > racon4.fasta

################################ Medaka generate consensus ##############################
mkdir medaka
cd medaka
medaka_consensus -i ../../../qc/qc.fastq.gz -d ../racon4.fasta -o . -t 14 -m r941_min_high_g303

echo "End: `date`; RC=$?"

```

### short reads polish using pilon
Nanopore assembled contig polish - short reads polish
```
#!/usr/bin/env bash

#BSUB -J pilon_polish
#BSUB -q coursework
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -o pilon_polish_%J.OUT
#BSUB -e pilon_polish_%J.ERR

echo "Start: `date`"

################################ four rounds of short reads polish ##############################
# illumina pilon polish
#round1
bwa index ../racon_polish/medaka/consensus.fasta
bwa mem -t 14 ../racon_polish/medaka/consensus.fasta ../../../illumina/qc/illumina_qc.R1.fastq ../../../illumina/qc/illumina_qc.R2.fastq  | samtools view - -Sb | samtools sort - -@14 -o illumina.mapping.sorted.bam
samtools index illumina.mapping.sorted.bam

pileup.sh in=illumina.mapping.sorted.bam ref=../racon_polish/medaka/consensus.fasta out=illumina_cov.txt
java -Xmx100G -jar pilon-1.23.jar --genome ../racon_polish/medaka/consensus.fasta --fix all --changes --frags illumina.mapping.sorted.bam --threads 14 --output pilon_round1 | tee round1.pilon

#round2
bwa index pilon_round1.fasta
bwa mem -t 14 pilon_round1.fasta ../../../illumina/qc/illumina_qc.R1.fastq ../../../illumina/qc/illumina_qc.R2.fastq   | samtools view - -Sb | samtools sort - -@14 -o illumina.mapping1.sorted.bam
samtools index illumina.mapping1.sorted.bam

pileup.sh in=illumina.mapping1.sorted.bam ref=pilon_round1.fasta out=map1_cov.txt
java -Xmx100G -jar /gpfs/ebg_data/programs/nanopore/pilon/pilon-1.23.jar --genome pilon_round1.fasta --fix all --changes --frags illumina.mapping1.sorted.bam --threads 14 --output pilon_round2 | tee round2.pilon


#round3
bwa index pilon_round2.fasta
bwa mem -t 14 pilon_round2.fasta ../../../illumina/qc/illumina_qc.R1.fastq ../../../illumina/qc/illumina_qc.R2.fastq | samtools view - -Sb | samtools sort - -@14 -o illumina.mapping2.sorted.bam
samtools index illumina.mapping2.sorted.bam
pileup.sh in=illumina.mapping2.sorted.bam ref=pilon_round2.fasta out=map2_cov.txt

java -Xmx100G -jar /gpfs/ebg_data/programs/nanopore/pilon/pilon-1.23.jar --genome pilon_round2.fasta --fix all --changes --frags illumina.mapping2.sorted.bam --threads 14 --out
put pilon_round3 | tee round3.pilon


#round4
bwa index pilon_round3.fasta
bwa mem -t 14 pilon_round3.fasta ../../../illumina/qc/illumina_qc.R1.fastq ../../../illumina/qc/illumina_qc.R2.fastq | samtools view - -Sb | samtools sort - -@14 -o illumina.mapping3.sorted.bam
samtools index illumina.mapping3.sorted.bam
pileup.sh in=illumina.mapping3.sorted.bam ref=pilon_round3.fasta out=map3_cov.txt
java -Xmx100G -jar /gpfs/ebg_data/programs/nanopore/pilon/pilon-1.23.jar --genome pilon_round3.fasta --fix all --changes --frags illumina.mapping3.sorted.bam --threads 14 --output pilon_round4 | tee round4.pilon

echo "End: `date`; RC=$?"

###########################
