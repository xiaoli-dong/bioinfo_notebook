#!/usr/local/bin/perl

use strict;
use warnings;
use threads;
use Benchmark;
  
@ARGV == 2 or die "Usage: $0 <sff file> <sample Name>\n";


my $sff = $ARGV[0];
my $sampleName = $ARGV[1];

my $t0 = Benchmark->new;
my $cmd = "";
my $sffinfo = "/export/data/linux/sffinfo";
#my $sffinfo = "/opt/454/bin/sffinfo";
if(! -e "$sampleName.fna" || ! -e "$sampleName.qual"){
    $cmd = "$sffinfo -s $sff > $sampleName.fna; $sffinfo -q $sff > $sampleName.qual";
    print STDERR "$cmd\n";
    system $cmd; 
}


my $mothur = "/data/metatools/Mothur_1_18/Mothur.source/mothur";
#get rid of short sequecnes, low quality reads and reads containing long homopolyers
if(! -e "$sampleName.trim.fasta" || ! -e "$sampleName.trim.qual"){
    $cmd = "$mothur \"\#trim.seqs(fasta=$sampleName.fna,qfile=$sampleName.qual,qaverage=25, maxambig=0, minlength=100, maxhomop=6)\"";
    print STDERR "$cmd\n";
    system $cmd;
}



#~/tools/usearch4.2.66_i86linux32 --sort ../GXTPUQ202.GV_TDS_TP6_30ft_2008.trim.fasta --output GXTPUQ202.GV_TDS_TP6_30ft_2008.trim.sorted.fasta
#get rid of artificial duplicates
#my $usearch = "/data/metatools/usearch6.0/usearch6.0.152_i86linux32";
my $usearch = "/data/metatools/usearch4.2.66_i86linux32";
if(! -e "$sampleName.seed.fasta"){
    $cmd = "$usearch -sort $sampleName.trim.fasta --output $sampleName.trim.sorted.fasta; ";
    $cmd .= "$usearch --cluster  $sampleName.trim.sorted.fasta --id 0.90 --seedsout  $sampleName.seed.fasta --uc $sampleName.seed.uc --idprefix 5; ";
    print STDERR "$cmd\n";
    system $cmd;
}


#perl -ne 'if(/^>(\w+)/){print "$1\n";}' Pond6_2010.seed.fasta > includes.txt
if(! -e "includes.txt"){
    $cmd = "perl -ne \'if(\/\^>(\\w+)\/){print \"\$1\\n\";}\' $sampleName.seed.fasta > includes.txt; ";
    print STDERR "$cmd\n\n";
    system $cmd;
}

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print STDERR "$sampleName took:",timestr($td)," to finish quality control\n";
