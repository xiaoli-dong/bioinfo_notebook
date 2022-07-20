#!/usr/bin/env perl

use warnings;
use strict;

$#ARGV == -1
  and die
"Get sequece ending pattern in a fasta file\nUsage: $0 <fastq file>  <cutoff_len>\n when the cufoff_len is negative, it extract the sequencing ending pattern; otherwise, it's checking the sequence starting pattern\n";

open( FASTQ, "<$ARGV[0]" )
  or die "Could not open input fasta file $ARGV[0] to read, $!\n";

my %pattern      = ();
my %pattern2seqs = ();

while (<FASTQ>) {

    my $line1 = $_;
    my $line2 = <FASTQ>;
    chomp($line2);
    $_ = <FASTQ>;
    $_ = <FASTQ>;
    my ($name) = $line1 =~ /^(\S+)/;
    my $seq = $line2;

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

close(FASTQ);

foreach ( sort { $pattern{$b} <=> $pattern{$a} } keys %pattern ) {

    print "$_\t$pattern{$_}\n";

    my @tmp = sort split( /:/, $pattern2seqs{$_} );

    #$tmp =~ s/:/\n/g;
    #print  join ("\n", @tmp), "\n";

}
