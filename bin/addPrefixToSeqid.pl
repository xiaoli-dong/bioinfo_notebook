#!/usr/bin/env perl
use strict;
use warnings;
$#ARGV == -1

    and die "Usage: $0 <fasta file>  <prefix name>\n";

my $fasta = $ARGV[0];
my $prefix = $ARGV[1];


open(FASTA, "<$fasta") or die "cannot open $fasta for reading: $!";

while(<FASTA>){
    if (/^>(\w+)/){
	print ">$prefix\_$1\n"
	}
    else{
	print $_;
    }
}

close(FASTA);


