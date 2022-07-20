#!/usr/bin/perl

$infile      = $ARGV[0];
$matchString = $ARGV[1];

chomp($infile);

open( INFILE, "<$infile" );

$/     = "\n>";
$count = 0;

while ( $entry = <INFILE> ) {

    chomp($entry);

    if ( ( $head, $seq ) = $entry =~ /^>?(.+?)\n(.*)/s ) {

        if ( $head =~ /$matchString/i ) {
            $count++;
            print ">$head\n$seq\n";
        }

    }

}

close(INFILE);

print STDERR "$count sequences has been retrieve\n";

