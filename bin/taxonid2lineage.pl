#!/usr/local/bin/env perl

BEGIN {
    # MAGPIE is required for taxonomic rank determination
    $ENV{MAGPIEHOME} = "/export/home/xdong/deve/metagenome/magpie";
    push @INC, "$ENV{MAGPIEHOME}/lib";
}

use Taxonomy;

$#ARGV == -1
  and die "Usage: $0  <tab delimited file, first column is txonid>\n";

my $tax = new Taxonomy();

open FILE, $ARGV[0];

while (<FILE>) {
    if (/^#/) {
        print $_;
        next;
    }
    chomp;

    my ( $taxida, $other ) = $_ =~ /^(\d+?)_.*?\t(\d+.*)$/;

    #print STDERR "$taxida\n";
    my @taxids = $tax->lineage($taxida);

    #print $_;
    foreach my $id (@taxids) {
        my $rank    = $tax->rank($id);
        my $sysname = $tax->taxid2sysname($id);
        $sysname =~ s/^\s+//;
        $sysname =~ s/\s+$//;
        next if $sysname eq "root" || $sysname eq "cellular organisms";
        print "$sysname;";

        #print  "$_=>$rank=>$sysname;\n";
    }
    print "\t$other\n";
}
