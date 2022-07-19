#!/usr/bin/env perl

use warnings;
use strict;

$#ARGV == -1
  and die
"Get sequece ending pattern in a fasta file\nUsage: $0 <fasta file>  <cutoff_len>\n when the cufoff_len is negative, it extract the sequencing ending pattern; otherwise, it's checking the sequence starting pattern\n";

open( FASTA, "<$ARGV[0]" )
  or die "Could not open input fasta file $ARGV[0] to read, $!\n";

$/ = "\n>";
my %pattern      = ();
my %pattern2seqs = ();

while (<FASTA>) {
    chomp;
    if ( !/>?(\S+).*?\n(.+)/s ) {

        #die "Could not read FastA record #$.: $_\n";
        next;
    }
    my $name = $1;
    my $seq  = $2;

    $seq =~ s/\s+//g;

    #print STDERR "name=$name\n$seq\n";

    if ( $ARGV[1] > 0 ) {

        my $start = substr $seq, 0, $ARGV[1];

        #print STDERR "$ARGV[1], start=$start\n";
        $pattern{$start}++;
        if ( exists $pattern2seqs{$start} ) {
            $pattern2seqs{$start} .= ":$name";
        }
        else {
            $pattern2seqs{$start} = $name;
        }

    }
    else {
        my $end = substr $seq, $ARGV[1];
        $pattern{$end}++;

        if ( exists $pattern2seqs{$end} ) {
            $pattern2seqs{$end} .= ":$name";
        }
        else {
            $pattern2seqs{$end} = $name;
        }
    }
}

$/ = "\n";

close(FASTA);

foreach ( sort { $pattern{$b} <=> $pattern{$a} } keys %pattern ) {

    print "$_\t$pattern{$_}\n";

    my @tmp = sort split( /:/, $pattern2seqs{$_} );

    #$tmp =~ s/:/\n/g;
    #print  join ("\n", @tmp), "\n";

}
