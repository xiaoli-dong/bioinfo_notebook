#!/usr/bin/perl -w

use strict;

my $infasta = $ARGV[0];

open FH1, "$infasta" or die "$!";
$/ = "\n>";

my %seq2gc = ();
while (<FH1>) {
    chomp;
    if ( !/>?(\w+?)\s+.*?\n(.+)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    my $seqid = $1;
    my $seq   = $2;
    $seq =~ tr/ \r\n\t//d;
    my $UC        = uc($seq);
    my $length    = ( $UC =~ tr/ATGCN// );
    my $numsGCs   = ( $UC =~ tr/GC/gc/ );
    my $percentGC = ( $numsGCs / $length ) * 100;
    my $rounded   = sprintf( "%.2f", $percentGC );

    #print "$seqid\t$rounded\n";
    $seq2gc{$seqid} = $rounded;
}
$/ = "\n";
close(FH1);

my $len       = 0;
my %len2count = ();

foreach ( sort { $seq2gc{$a} <=> $seq2gc{$b} } keys %seq2gc ) {

    if ( $seq2gc{$_} < $len ) {
        $len2count{$len}++;
    }
    else {
        $len += 5;
        $len2count{$len}++;
    }
}

#produce the data for chart
print "var data = [\n";

my $i        = 1;
my $keyCount = keys %len2count;

foreach ( sort { $a <=> $b } keys %len2count ) {

    if ( $i < $keyCount ) {
        print "[$_, $len2count{$_}],\n";
    }
    else {
        print "[$_,$len2count{$_}]\n";
    }
    $i++;

}
print "];\n";
