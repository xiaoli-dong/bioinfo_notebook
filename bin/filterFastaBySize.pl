#!/usr/bin/env perl

use strict;
use warnings;

@ARGV == 3 or die "Usage: $0 <src_dir><target_dir><suffix>\n";

my $src_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $suffix = $ARGV[2];
my @files = glob "$src_dir/*\.$suffix";

for my $file (@files){
    open(SEQ, "$file") or die "Could not open $file file to read, $!\n";
    my $tp = 0;
    my $seq_dir = $file;
    $seq_dir =~ s/\.$suffix$//;
    while(<SEQ>){
        next if /^>/;
        s/\s+//g;
        $tp += length($_);
    }
    if($tp >= 2000000){
        system("mkdir $target_dir/$seq_dir; cp $src_dir/$file $target_dir/$seq_dir");
    }
    else{
        print STDERR "$src_dir/$file is less than 2M bp\n";
    }
    close(SEQ);
}