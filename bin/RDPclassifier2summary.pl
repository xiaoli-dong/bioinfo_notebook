#!/usr/local/bin/perl
@ARGV == 3
  or die
"Usage: $0 <RDPClassifier output> <simplified RDP classifier> <taxon summary>\n",
  "generate taxonomic summary from RDP classifier output\n";

use strict;
use warnings;
open( RDP, "$ARGV[0]" )  or die "Could not open $ARGV[0] to read, $!\n";
open( SIM, ">$ARGV[1]" ) or die "Can't open $ARGV[1] for writting: $1\n";
open( OUT, ">$ARGV[2]" ) or die "Can't open $ARGV[2] for writting: $1\n";

my %taxons = ();
my $count  = 0;
while (<RDP>) {
    if (/^(\S+?\s+?)(\S+?)(unclassified|uncultured)/) {
        my $tax   = $2;
        my $seqid = $1;
        print SIM "$seqid$tax\n";
        $tax =~ s/\(\d+\)//g;
        $taxons{$tax} += 1;
        $count++;
    }
}
close(SIM);
close(RDP);
print OUT "#Total classified sequecnes: $count\n";
print OUT "#The first colum is Taxon_name\n";
print OUT "#The second colum is taxon_count\n\n";

foreach ( sort { $taxons{$b} <=> $taxons{$a} } keys %taxons ) {

    print OUT "$_\t$taxons{$_}\n";

}
close(OUT);
