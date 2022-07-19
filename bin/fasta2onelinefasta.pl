#!/usr/bin/perl

# AUTHOR: Joseph Fass
# LAST REVISED: October 2009
#
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2009 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;

my $usage =
    "\nusage: cat <input fasta file> | $0 > <output fasta file>\n\n"
  . "For each sequence, puts nt/aa all on one line following header line.\n"
  . "Reads from and writes to STDOUT, so use redirection ... \n\n\n";

my $block = "";
my $line;

while (<>) {
    $line = $_;
    if (m/^>/) {    # header line
        print $block;    # which is either empty, or capped by a newline
        $block =
          $line . "\n"; # adds another newline, because it will be chomped below
    }
    else {
        chomp $line;
        if ( substr( $line, -1 ) =~ m/[0-9]/ ) {
            chomp $block;       # remove newline before adding more sequence
            $block .=
              $line . " \n";    # must be fasta qualities, and needs a space
        }
        else {
            chomp $block;           # remove newline before adding more sequence
            $block .= $line . "\n"; # not fasta qualities, so no space needed
        }
    }
}
print $block;    # print last sequence block
