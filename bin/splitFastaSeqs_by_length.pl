#!/usr/bin/env perl
use strict;
use warnings;

$#ARGV == -1

  and die
"Usage: $0 <fasta file> <length> <output1 file> <output2 file> spliting fasta sequence file by length, the sequences less than or equal to length will go to the first file, otherwise go to second file\n";

$/ = "\n>";
my $len   = $ARGV[1];
my $le    = 0;
my $g     = 0;
my $total = 0;

open( FASTA, "$ARGV[0]" )  or die "Could not open $ARGV[0] to read, $!\n";
open( OUT1,  ">$ARGV[2]" ) or die "Could not open $ARGV[2] to write, $!\n";
open( OUT2,  ">$ARGV[3]" ) or die "Could not open $ARGV[3] to write, $!\n";

while (<FASTA>) {
    chomp;
    if ( !/>?(.*?)\n(.+)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    $total++;
    my $id  = $1;
    my $seq = $2;

    $seq =~ tr/ \r\n\t//d;

    if ( length($seq) <= $len ) {

        $le++;
        print OUT1 ">$id\n$seq\n";
    }
    else {
        $g++;
        print OUT2 ">$id\n$seq\n";
    }

}

print STDERR
"Total sequences: $total\nThe number of sequences shorter or equal to $len bp is: $le\nThe number of sequences longer than $len bp is: $g\n";
close(FASTA);
close(OUT1);
close(OUT2);

