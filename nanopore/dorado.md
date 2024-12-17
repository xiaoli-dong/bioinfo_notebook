Basecalling:
```
# reads are classified into their barcode groups during basecalling
dorado basecaller -x cuda:all --recursive --barcode-both-ends  --kit-name  SQK-NBD114-96 --min-qscore 9  --
no-trim  ../dorado-0.8.3-linux-x64/dorado-models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0  COVID/20241126_1409_MN24257_FAX85983_d5fd5795/pod5_p
ass > dorado_out/calls.bam

dorado demux --no-classify -t 8 -v --emit-fastq --emit-summary calls.bam fastq
```
