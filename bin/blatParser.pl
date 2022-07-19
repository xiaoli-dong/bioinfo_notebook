#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $qlen = 0.90;

my ($blat_output);
&GetOptions(
	    "b=s" =>\$blat_output,
	    "ql=s" =>\$qlen
	    );

($blat_output ) ||
    die "usage: $0 OPTIONS

where options are:  -b  <psl format blat output> -ql <query length threadhold>\n";

open(PSL, "<$blat_output") or die "Could not open $blat_output file to read, $!\n";

my %hits = ();
my %len = ();
print "QueryName,QueryLength,Matches,Strand,TargetName,TargetStart,TargetEnd\n";
while (<PSL>){
    next if !/^\d+/;
    my @arr=split("\t", $_);
    if((not exists $hits{$arr[9]}) || (exists $hits{$arr[9]} && $len{$arr[9]} < $arr[0]) ){
	if($arr[0]/$arr[10] > $qlen){
	    $len{$arr[9]} = $arr[0];
	    $hits{$arr[9]} = "$arr[10],$arr[0],$arr[8],$arr[13],$arr[14],$arr[15]\n";
    }
    }
    
}
close(PSL);

foreach (sort keys %hits){

    print "$_,$hits{$_}\n";
    
}
