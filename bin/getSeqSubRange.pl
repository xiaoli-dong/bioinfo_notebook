#!/usr/bin/env perl
use strict;
use warnings;
$#ARGV == -1

  and die "Usage: $0 <fasta file> <contigid><startPos 0 based><endPos>\n";

my $fasta = $ARGV[0];
my $seqid = $ARGV[1];
my $start = $ARGV[2];
my $end = $ARGV[3];
open( FASTA, "<$fasta" ) or die "cannot open $fasta for reading: $!";

my $i = 0;
$/ = "\n>";
while (<FASTA>) {
    chomp;

    if ( my ( $seq_name, $seq ) = /^>?(\S+).*?\n(.*)/s ) {
        $seq =~ tr/ \r\n\t//d;
        if ( $seq_name eq $seqid ) {
            my $len = $end - $start + 1;
            my $range = substr( $seq, $start, $len );
            print STDERR "length=", length($range), ", $len\n";
            print ">$seq_name\n$range\n";
            last;
        }
    }
}
$/ = "\n";
close(FASTA);

