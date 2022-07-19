#!/usr/local/bin/perl

$#ARGV == -1

  and die "Usage: $0  hmm search nt database output\n";

open( HMMOUT, "$ARGV[0]" ) or die "Could not open $ARGV[0] to read, $!";

$/ = "\n\___";

while (<HMMOUT>) {

    if (/T\s+\=\s+(\w+).*?TS =\s+?(\d+)\s+TE =\s+(\d+).*?\n(Q .*)$/so) {

        print ">$1\_$2\_$3\n";

        my @seqs = split( /\n/, $4 );
        foreach my $seq (@seqs) {
            if ( $seq =~ /^T\s+\d+/ ) {
                $seq =~ s/[T\s+\-\d+]//g;
                print "$seq";
            }
        }
        print "\n";
    }

}
close(HMMOUT);
