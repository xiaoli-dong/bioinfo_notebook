#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my ( $fasta, $minlen, $maxlen );

&GetOptions(
    "i=s" => \$fasta,
    "min=i" => \$minlen,
    "max=i" => \$maxlen
);

($fasta)
  || die "Filter sequence by length
    usage: $0 OPTIONS
    where options are:
        -i  <fasta sequence file>
        -min <sequence min length> 
        -max <sequence max length, default is -1>\n";

$minlen ||= 0;
$maxlen ||= -1;


my $i = 0;

open( FASTA, "<$fasta" ) or die "Could not open $fasta to read, $!\n";
$/ = "\n>";
while (<FASTA>) {
    chomp;

    if ( my ( $seq_name, $seq ) = /^>?(\S+.*?)\n(.*)/s ) {

        $seq =~ tr/ \r\n\t//d;
        my $seqlen = length($seq);
        if ( $seqlen >= $minlen ) {

            if ( $maxlen != -1 ) {
                if ( $seqlen <= $maxlen ) {
                    print ">$seq_name\n$seq\n";
                }
            }
            else {
                print ">$seq_name\n$seq\n";

            }
        }
    }
}

$/ = "\n";

close(FASTA);

