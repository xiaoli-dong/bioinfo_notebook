Download sra files using sratoolkit:

mkdir sra

cd sra

 ~/software/sratoolkit/sratoolkit.2.11.0-centos_linux64/bin/prefetch --option-file acc.txt
 
for var in sra/*.sra; do fastq-dump --outdir fastq --split-3 $var; done

Download metadata
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&retmode=json&id=ERR1034587


How to mount windows network drives in wsl:
https://www.public-health.uiowa.edu/it/support/kb48568/
