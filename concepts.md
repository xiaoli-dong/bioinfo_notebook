## What is decoy human genome database 

The reference genome is incomplete, particularly around the centromeres, so often reads which truly belong elsewhere are wrongly mapped to a particular place in the genome because the true match is missing from the reference. These cause false positive calls, which were bothering us in the 1000 Genomes Project. The decoy is a pragmatic solution to this - it contains known true human genome sequence that is not in the reference genome, and will "suck up" reads that would otherwise map with low quality in the reference. The decoy was built by Heng Li, at the Broad, working with Richard Durbin (Sanger) and Deanna Church of the Genome Reference Consortium (who maintain the reference genome).

[here is the link](https://www.biostars.org/p/73100/#:~:text=The%20decoy%20is%20a%20pragmatic,low%20quality%20in%20the%20reference)
[The decoy genome](https://www.cureffi.org/2013/02/01/the-decoy-genome/)
