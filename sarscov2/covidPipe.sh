#!/bin/bash

##  Written by Xiaoli Dong
##  Last modified Dec 1, 2020
##  Description:  Calls SARS-CoV-2 variants from Illumina amplicon data.
##                This script assumes reads data is paired-end.

##  Usage:  covidPipe.sh <prefix>
##  For example, "covidPipe.sh Sample1" if the data is in Sample1_R*.fastq.gz
## sh covidPipe.sh AB-99850 ../raw | tee -a log.txt
##  Grab the sample name from the command line.
NAME="$1"
INPUT_DIR="$2"

##  Set minimum coverage for genotype calls.
##  Areas below this depth will be set to N in the consensus genome.
MINCOV=5

##  Specify the viral reference file.
##  NC_045512.fasta contains the SARS-CoV-2 genome, equivalent to bbmap/resources/Covid19_ref.fa
REF="../ref/nCoV-2019.reference.fasta"
GFF="../ref/nCoV-2019.gff"
##  PCR primers
##  artic3 primers v3 in bed format 
## https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3
PRIMERS_BED="../ref/nCoV-2019.primer.bed"
PRIMERS_TSV="../ref/nCoV-2019.tsv"
PRIMERS_FASTA="../ref/nCoV-2019.primer.fasta"
PRIMERS_PAIR="../ref/nCoV-2019.primer.pair.tsv"
## create primer fasta
#perl -ne '@a = split(/\t/, $_); print ">$a[0]\n$a[2]\n";' ${PRIMERS_TSV} > ${PRIMERS_FASTA}

cpus=8
#iVar Params

#Length of illumina reads to keep after primer trimming
illuminaKeepLen=50
# Sliding window quality threshold for keeping reads after primer trimming (illumina)
illuminaQualThreshold=20
# Mpileup depth for ivar
mpileupDepth=100000
# iVar frequency threshold for consensus variant (ivar consensus: -t)
ivarFreqThreshold=0.75
# Minimum coverage depth to call variant (ivar consensus: -m; ivar variants -m) 
ivarMinDepth=20
# iVar frequency threshold to call variant (ivar variants: -t )
ivarMinFreqThreshold=0.03
# iVar minimum mapQ to call variant (ivar variants: -q)
ivarMinVariantQuality=20

## prepare output directory 
mkdir ${NAME}
cd ${NAME}
ln -s ../${INPUT_DIR}/${NAME}_R1.fastq.gz
ln -s ../${INPUT_DIR}/${NAME}_R2.fastq.gz

for dir in aligned2consensus  aligned2ref  atrim  consensus   primer2consensus
do
mkdir $dir
done


##  Discover adapter sequence for this library based on read overlap.
##  You can examine the adapters output file afterward if desired;
##  If there were too few short-insert pairs this step will fail (and you can just use the default Illumina adapters).
bbmerge.sh in1=${NAME}_R1.fastq.gz in2=${NAME}_R2.fastq.gz outa=atrim/${NAME}_adapters.fa ow 

##  Perform adapter-trimming on the reads.
##  Also do quality trimming and filtering.
##  If desired, also do primer-trimming here by adding, e.g., 'ftl=20' to to trim the leftmost 20 bases.
##  If the prior adapter-detection step failed, use "ref=adapters"
bbduk.sh in1=${NAME}_R1.fastq.gz in2=${NAME}_R2.fastq.gz out=atrim/${NAME}_atrimmed_R1.fastq.gz out2=atrim/${NAME}_atrimmed_R2.fastq.gz minlen=${illuminaKeepLen} ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=atrim/${NAME}_adapters.fa altref=adapters maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx1g ftm=5

##  Align reads to the reference.
##  Local flag is due to primer-amplification-induced anomalies at read ends;
##  for randomly-sheared, unamplified data, "local" should be omitted.
bwa index ${REF}
bwa mem -t ${cpus} ${REF} atrim/${NAME}_atrimmed_R1.fastq.gz atrim/${NAME}_atrimmed_R2.fastq.gz | samtools view -F 4 -Sb | samtools sort -o aligned2ref/${NAME}.sorted.bam
samtools index aligned2ref/${NAME}.sorted.bam

ivar trim -b ${PRIMERS_BED} -p aligned2ref/${NAME}.ptrimmed -i aligned2ref/${NAME}.sorted.bam -m ${illuminaKeepLen} -q ${illuminaQualThreshold} -s 4
samtools sort aligned2ref/${NAME}.ptrimmed.bam -o aligned2ref/${NAME}.ptrimmed.sorted.bam
samtools index aligned2ref/${NAME}.ptrimmed.sorted.bam


## call variants
echo "Calling variants before remove reads having mismatches in the primer part"
samtools mpileup -A -d 0 --reference ${REF} -B -Q 0 aligned2ref/${NAME}.ptrimmed.sorted.bam | ivar variants -p ${NAME}.ivar -q ${ivarMinVariantQuality}  -t ${ivarMinFreqThreshold} -m ${ivarMinDepth} -r ${REF} -g ${GFF}

## call consensus
ehco "calling consensus"
samtools mpileup -A -d ${mpileupDepth} -Q 0 -F 0 aligned2ref/${NAME}.ptrimmed.sorted.bam | ivar consensus -p consensus/${NAME}.consensus -t ${ivarFreqThreshold} -m ${ivarMinDepth} -n N
bwa index consensus/${NAME}.consensus.fa


## align primer to consensus to create primer bam and primer bed
echo "align primer to consensus and identify the mismached primer"
bwa mem -k 5 -T 16 consensus/${NAME}.consensus.fa ${PRIMERS_FASTA} | samtools view -bS -F 4 -o primer2consensus/primers.bam
samtools sort -o primer2consensus/primers.sorted.bam primer2consensus/primers.bam
bedtools bamtobed -i primer2consensus/primers.sorted.bam > primer2consensus/primers.bed

## call_variants_in_primer
echo "Calling variants in primer"
samtools mpileup -A -d 0 --reference consensus/${NAME}.consensus.fa -Q 0 -F 0 primer2consensus/primers.sorted.bam | ivar variants -p primer2consensus/primer_mismatches.tsv -t ${ivarMinFreqThreshold} 

## get masked
ivar getmasked -i primer2consensus/primer_mismatches.tsv -b primer2consensus/primers.bed  -f ${PRIMERS_PAIR} -p primer2consensus/${NAME}.masked_primer_names.txt


################### call variants with known reference ################
## remove reads
ivar removereads -i aligned2ref/${NAME}.ptrimmed.sorted.bam -p aligned2ref/${NAME}.masked.bam  -t primer2consensus/primer_mismatches.masked.txt -b primer2consensus/primers.bed
samtools sort -o aligned2ref/${NAME}.masked.sorted.bam aligned2ref/${NAME}.masked.bam
samtools index aligned2ref/${NAME}.masked.sorted.bam

#call_variants_post_removal:
echo "calling variants with known reference"
samtools mpileup -A -d 0 --reference ${REF} -Q 0 -F 0 aligned2ref/${NAME}.masked.sorted.bam | ivar variants -p ${NAME}_with_known_ref.masked.variants.tsv  -t ${ivarMinFreqThreshold}


######## with unknown reference ##########

#align adapter trimmed reads to consensus
bwa mem -t ${cpus} consensus/${NAME}.consensus.fa atrim/${NAME}_atrimmed_R1.fastq.gz atrim/${NAME}_atrimmed_R2.fastq.gz | samtools view -F 4 -Sb | samtools sort -o aligned2consensus/${NAME}.sorted.bam
samtools index aligned2consensus/${NAME}.sorted.bam

####### trim_primer_after_realign
ivar trim -b primer2consensus/primers.bed -p aligned2consensus/${NAME}.ptrimmed.bam -i aligned2consensus/${NAME}.sorted.bam
samtools sort -o aligned2consensus/${NAME}.ptrimmed.sorted.bam aligned2consensus/${NAME}.ptrimmed.bam
samtools index aligned2consensus/${NAME}.ptrimmed.sorted.bam

## remove reads with mismatches in the primer region
ivar removereads -i aligned2consensus/${NAME}.ptrimmed.sorted.bam -p aligned2consensus/${NAME}.masked.bam  -t primer2consensus/${NAME}.masked_primer_names.txt -b primer2consensus/primers.bed
samtools sort -o aligned2consensus/${NAME}.masked.sorted.bam aligned2consensus/${NAME}.masked.bam
samtools index aligned2consensus/${NAME}.masked.sorted.bam

#call_variants_post_removal:
samtools mpileup -A -d 0 --reference consensus/${NAME}.consensus.fa -Q 0 -F 0 aligned2consensus/${NAME}.masked.sorted.bam | ivar variants -p ${NAME}_with_unknown_ref.masked.variants.tsv  -t ${ivarMinFreqThreshold} 

