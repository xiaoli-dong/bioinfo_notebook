# How to use [prefectch](https://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format/prefetch) from sratoolkit to download data 

Download sra files using sratoolkit:

```
mkdir sra
cd sra
prefetch --option-file acc.txt
for var in sra/*.sra; do fastq-dump --outdir fastq --split-3 $var; done
```

# Url used to download metadata
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&retmode=json&id=ERR1034587

# Stats tools
[GUide to STatistical Analysis in Microbial Ecology](http://mb3is.megx.net/gustame)
