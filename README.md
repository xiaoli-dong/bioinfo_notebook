Download sra files:
for var in sra/*.sra; do fastq-dump --outdir fastq --split-3 $var; done
