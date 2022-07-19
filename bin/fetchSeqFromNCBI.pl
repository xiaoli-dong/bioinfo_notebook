#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

use Bio::DB::GenBank;

my ($gilistFile);
&GetOptions(
    "gi=s" => \$gilistFile

);

($gilistFile)
  || die "usage: $0 OPTIONS where options are:  -gi <gi number list file>\n";
my $gb = Bio::DB::GenBank->new();
open( GI, "$gilistFile" ) or die "could not read $gilistFile to read\n";
my $list  = "";
my $count = 0;
while (<GI>) {
    if (/^(\d+)/) {
        my $gi = $1;
        if ($list) {
            $list .= ",  \'$gi\'";
        }
        else {
            $list = "\'$gi\'";
        }
        $count++;

        if ( $count == 100 ) {

            my $seqio = $gb->get_Stream_by_acc( [$list] );
            while ( my $clone = $seqio->next_seq ) {

         #print "cloneid is ", $clone->display_id, " ", $clone->primary_id, " ",

                #$clone->accession_number, "\n";
                print ">gi|", $clone->primary_id(), "\| ", $clone->desc(),
                  "\n", $clone->seq(), "\n";

            }
            $count = 0;
            $list  = "";

        }
    }
}
print STDERR "outside, list = $list\n";
my $seqio = $gb->get_Stream_by_acc( [$list] );
while ( my $clone = $seqio->next_seq ) {

    print ">gi|", $clone->primary_id(), "\| ", $clone->desc(), "\n",
      $clone->seq(), "\n";
}
