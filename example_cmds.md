### How to use [prefectch](https://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format/prefetch) from sratoolkit to download data 

Download sra files using sratoolkit:

```
mkdir sra
cd sra
prefetch --option-file acc.txt
for var in sra/*.sra; do fastq-dump --outdir fastq --split-3 $var; done
```

### Url used to download metadata
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&retmode=json&id=ERR1034587

### Stats tools
[GUide to STatistical Analysis in Microbial Ecology](http://mb3is.megx.net/gustame)

### Guppy basecall example
```
nohup /opt/ont/guppy/bin/guppy_basecaller -i fast5_pass -s guppy_v_345/fast5_pass --num_callers 4 --config dna_r9.4.1_450bps_hac.cfg --device cuda:all:100% --qscore_filt
ering --min_qscore 7 >& debug.pass.txt&
```
