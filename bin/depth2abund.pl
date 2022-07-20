#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 1 or die "Usage: $0 <depth file>\n";
my $depthf = $ARGV[0];


open(DEPTH, $depthf) or die "Could not open $depthf, !$\n";
open (LIST, ">abund_list.txt");
my @header = ();
my %abund = ();

while(<DEPTH>){
    chomp;
    if(/^contigName.*totalAvgDepth\s+(\S.*)/){
	@header = split(/\t/, $1);
    }
    else{
	my @l = split(/\t/, $_);
	my $cid = shift @l;
	shift @l;
	shift @l;
	
	for my $i (0 .. $#l) {
	    if($i % 2 == 0){
		$abund{$cid}->{$header[$i]} = $l[$i];
	    }
	}
    }
}

close(DEPTH);

my %fhs = ();

for my $i (0 .. $#header) {
    
    if($i % 2 == 0){
	
	my $file_out = "./$header[$i].abund.txt";
	open $fhs{$header[$i]},'>', $file_out or die $!;
	print LIST $file_out, "\n";
    }
}

close(LIST);

for my $cid (keys %abund){
    
    for my $h (keys %{$abund{$cid}}){
	
	print {$fhs{$h}} "$cid\t$abund{$cid}->{$h}\n";
    }
}

foreach (keys %fhs){

    close($fhs{$_});
}
