Download sra files using sratoolkit:

mkdir sra

cd sra

 ~/software/sratoolkit/sratoolkit.2.11.0-centos_linux64/bin/prefetch --option-file acc.txt
 
for var in sra/*.sra; do fastq-dump --outdir fastq --split-3 $var; done
