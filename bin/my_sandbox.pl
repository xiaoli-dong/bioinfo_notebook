#!/usr/bin/env perl

use strict;
use warnings;

my $seq1  = "AAAAAACGATCTTTTTTTTTCCC";
my $qual1 = "GGGGGGFFFFFFFFFFFFFFFFF";
print "$seq1\n$qual1\n";
$seq1 =~ s/^A{3,}|T{3,}|C{3,}|G{3,}//;
$qual1 = substr( $qual1, $+[0] );

$seq1 =~ s/(A{3,}|T{3,}|C{3,}|G{3,})$//;
$qual1 = substr( $qual1, 0, $-[0] );

print "$seq1\n$qual1\n";
