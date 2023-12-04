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

## Nanopore 16s rRNA amplicon 
  * Emu algorithm <p>"Emu was given its name because it uses an expectation-maximization (EM) approach for microscopic (mu, 𝛍) organisms, yielding EM + mu or Emu for short. This algorithm is constructed upon the idea that an unclassified sequence is more likely to come from an organism that is expected to be in the sample at high abundance rather than an organism that is either detected in low abundance or not at all. Emu first produces read classifications and a direct proportion community profile that is expected to have some inaccuracies due to error-prone sequences. The initial guess profile is then used to update the read classification likelihood giving more weight to higher abundance species. These updated read classification likelihoods then update the community profile directly. This process continues until only marginal changes are made to the community profile between iterations. A clean-up step is completed to remove low likelihood species and the final community profile estimation is ultimately returned." --[Emu: a novel computational tool for more precise community profiles from full-length 16S sequences](https://microbiologycommunity.nature.com/posts/emu-a-novel-computational-tool-for-more-precise-community-profiles-from-full-length-16s-sequences)</p>
<div id="image-table">
    <table>
	    <tr>
    	    <td style="padding:10px">
        	    <img src="https://github.com/xiaoli-dong/bioinfo_notebook/assets/52679027/f693e227-3127-47fe-adf9-6b1af601716d"/>
      	    </td>
            <td style="padding:10px">
            	<img src="https://github.com/xiaoli-dong/bioinfo_notebook/assets/52679027/4b054101-65c9-4464-90d9-58b87374f9a6">
            </td>
        </tr>
    </table>
</div>
 
## Sample-Barcode Bleeding

"As nanopore sequencing is a molecule-based approach, it is easy to “over- or underload” cDNA molecules for library preparation. This overloading often happens if the majority (>50%) of all fragments are “super small” (<100 nucleotides) in the sample. Molecule-overloading usually leads to unspecific barcode ligation during the adapter ligation step after barcode pooling. Free barcodes from samples with a low amount of DNA molecules might ligate to non-barcoded DNA from samples with a highly overloaded number of molecules. A negative control can visualize this issue (poreCov highlights negative controls in the report if specified). Furthermore, using only reads with barcodes present on both ends mitigates the problem (the default setting for poreCov). If sample bleeding occurred, a false-positive negative control usually inherits a similar SNP proportion pattern as the barcoded samples used on the same run (see Figure 3). This can be visually inspected by opening the .bam file in the result dir “2.Genomes” via, e.g., Unipro UGENE (Okonechnikov et al., 2012) or IGV viewer (Robinson et al., 2011) "  -- [poreCov-An Easy to Use, Fast, and Robust Workflow for SARS-CoV-2 Genome Reconstruction via Nanopore Sequencing](https://www.frontiersin.org/articles/10.3389/fgene.2021.711437/full)

