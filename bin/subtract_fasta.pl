#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ( $seqFile, $queryFile );

&GetOptions(
    "s=s" => \$seqFile,
    "q=s" => \$queryFile
);

($seqFile)
  || die
  "usage: get sequences which was not included in the queryFile
  $0 OPTIONS
    where options are:
    -s  <input total sequence file>
    -q < sequences will be excluded from the input total sequence file>\n";

my %names = ();
open( QUERY, "<$queryFile" ) or die $!;

while (<QUERY>) {
    if (/^>?(\S+)/) {
        $names{$1} = 1;
    }
}
close(QUERY);

open( SEQ, "<$seqFile" ) or die "cannot open $seqFile for reading: $!";

my $i = 0;

$/ = "\n>";
while (<SEQ>) {
    chomp;

    if ( my ( $seq_name, $other, $seq ) = /^>?(\S+)(.*?)\n(.*)/s ) {

        if ( not exists $names{$seq_name} ) {
            print ">$seq_name$other\n$seq\n";
        }

    }

}

$/ = "\n";

close(SEQ);
