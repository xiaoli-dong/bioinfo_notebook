#!/usr/local/bin/perl

use strict;
use warnings;

@ARGV == 3 or die "Usage: $0 <fastq read1> <fastq read2> <prefix of output>\n",
                  "trim read2 and read1 to to the same length, the lenght will be the shortest between R1 and R2\n";

open(R1, $ARGV[0]) or die "Cannot open $ARGV[0] for reading: $!\n";
open(R2, $ARGV[1]) or die "Cannot open $ARGV[1] for reading: $!\n";
open(R1W, ">$ARGV[2].R1.fastq") or die "Cannot open $ARGV[2].R1.fastq for writing: $!\n";
open(R2W, ">$ARGV[2].R2.fastq") or die "Cannot open $ARGV[2].R2.fastq for writing: $!\n";

while(<R1>){
    
    my $head1 = $_;
    my $seq1 = <R1>;
    chomp($seq1);
    <R1>;
    my $qual1 = <R1>;
    chomp($qual1);
    
    my $head2 = <R2>;
    my $seq2 = <R2>;
    chomp($seq2);
    <R2>;
    my $qual2 = <R2>;
    chomp($qual2);
    
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    
    if($len1 < $len2){
	print R1W $head1, $seq1, "\n+\n", $qual1, "\n";
	print R2W $head2, substr($seq2, 0, $len1),  "\n+\n", substr($qual2, 0, $len1), "\n";
	print STDERR "R2Adapter\t", substr($seq2,$len1-1), "\n";
    }
    elsif($len1 > $len2){
	print R2W $head2, $seq2, "\n+\n", $qual2, "\n";
	print R1W $head1, substr($seq1, 0, $len2),  "\n+\n", substr($qual1, 0, $len2), "\n";
	print STDERR "R1Adapter\t", substr($seq1,$len2-1), "\n";
    }
    else{
	print R1W $head1, $seq1, "\n+\n", $qual1, "\n";
	print R2W $head2, $seq2, "\n+\n", $qual2, "\n";
    }
    
}

close(R1);
close(R2);

close(R1W);
close(R2W);
