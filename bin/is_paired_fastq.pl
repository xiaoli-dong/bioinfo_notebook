#!/usr/bin/env perl

use strict;
use warnings;
$#ARGV == -1

  and die "Usage: $0 <fastq1>  <fastq2> \n";

open( FASTQ1, "$ARGV[0]" ) or die "Could not open $ARGV[0] to read, $!\n";
open( FASTQ2, "$ARGV[1]" ) or die "Could not open $ARGV[1] to read, $!\n";

my $count = 0;
my $total = 0;

while (<FASTQ1>) {

    my $seqhead1 = $_;
    $_ = <FASTQ1>;
    $_ = <FASTQ1>;
    $_ = <FASTQ1>;

    $seqhead1 =~ /^@(\S+?)\/\d+$/;
    my $id1       = $1;
    my $seqhead2  = <FASTQ2>;
    my $seq2      = <FASTQ2>;
    my $qualhead2 = <FASTQ2>;
    my $qual2     = <FASTQ2>;

    $seqhead2 =~ /^@(\S+?)\/\d+/;
    my $id2 = $1;
    $count++;
    if ( $id1 ne $id2 ) {

        print "$count pair: $id1 and $id2 are not matching\n";
        exit();
    }

}
print "they are perfact pairs\n";
close(FASTQ1);
close(FASTQ2);

