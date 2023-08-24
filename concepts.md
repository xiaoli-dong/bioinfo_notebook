## What is decoy human genome database 

The reference genome is incomplete, particularly around the centromeres, so often reads which truly belong elsewhere are wrongly mapped to a particular place in the genome because the true match is missing from the reference. These cause false positive calls, which were bothering us in the 1000 Genomes Project. The decoy is a pragmatic solution to this - it contains known true human genome sequence that is not in the reference genome, and will "suck up" reads that would otherwise map with low quality in the reference. The decoy was built by Heng Li, at the Broad, working with Richard Durbin (Sanger) and Deanna Church of the Genome Reference Consortium (who maintain the reference genome).

- [the decoy genome](https://www.cureffi.org/2013/02/01/the-decoy-genome/)

## primary, secondary, and supplementary alignment in MiniMap2
A secondary alignment occurs when a given read could align reasonably well to more than one place. One of the possible reported alignments is termed "primary" and the others will be marked as "secondary".When the two alignments are equally good. Which one gets flagged as "primary" is entirely random. In the minimap2, you can use --secondary=no to suppress the secondary alignment

supplementary alignment is the split or chimeric alignment. You can't suppress supplementary alignment  You can use samtools to filter out these alignments:
```
minimap2 -ax map-out ref.fa reads.fq | samtools view -F0x900
```

However, this is discouraged as supplementary alignment is informative.

## Measuring the quality of reads?
"Note that FastQC is designed for fixed-length short reads; it doesn't perform well with the dynamic and long read lengths seen in nanopore sequencing. Because nanopore read quality tends to be fairly consistent throughout a read, a biplot (or contour plot) of mean quality vs read length would be more representative of the actual error distribution." --from nanopore community

"The quality of the very ends is always of lower quality due to the relatively uncontrolled translocation speed. But the adapters act as kind of quality buffer. Once the motor settles into a fairly stable translocation speed you are then progressing into the read proper. " -- from nanopore community 
