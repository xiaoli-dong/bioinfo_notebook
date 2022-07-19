#!/usr/local/bin/perl


# mask the tRNA genes from input sequences
$#ARGV == -1

    and die "Usage: $0 <hmm against nt database search result>\n";

open(INPUT, "$ARGV[0]");
open(OUTPUT, ">$ARGV[0].summary.txt");

$/ = "\n\[QUERY\ SEQUENCE\]";

while(<INPUT>){
    if(/\[QUERY LOCUS\]\s+?(\S+.*?)\n.*?E_VALUE DESCRIPTION\s+(\S.*?)\s+\[BEGIN ALIGNMENTS TIME\]/s){
	print STDERR "$1\n";
	print  OUTPUT "Model Name= $1\n\n$2\n";
	my $list = $2;
	open(SEQS, ">$ARGV[0]\_$1.seqlist.txt");
	foreach my $record (split(/\n/, $list)){
	    #1    236.28  1 G1BJ6HH01CETH7                               D TP62010R  3.6e-066 length=506 xy=0871_3149 region=1 run=R_2011_04_20
	    if($record =~ /\d+\s+\S+\s+\d+\s+(\S+)\s+.*/){
		
		print SEQS"$1\n";
	    }
	    
	}
	close(SEQS);
	
    }
}
close(INPUT);
