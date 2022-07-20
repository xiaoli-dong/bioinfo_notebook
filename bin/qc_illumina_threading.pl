#!/usr/local/bin/perl

use strict;
use warnings;
use threads;
use Benchmark;
use Getopt::Long;
  
my $threads = 10;
my $min_len = 50;
my $headcrop_len = 12;
my $qcutoff = 20;
my $window_size = 5;
my $phred = 64;

my ($file1, $file2, $sampleName);
&GetOptions(

	    "s1=s" =>\$file1,
	    "s2=s" =>\$file2,
	    "sn=s" =>\$sampleName,
	    "t=i" =>\$threads,
	    "minlen=i"    => \$min_len,
	    "hdclen=i"    => \$headcrop_len,
	    "q=i"    => \$qcutoff,
	    "w=i"  =>\$window_size,
	    "p=i" =>\$phred
	    #"db=s"    => \$blastdb,
	    
	    );

($file1 && $file2 && $sampleName) ||
    die "usage: $0 OPTIONS

where options are:  -s1  seqfile1 -s2 seqfile2 -sn sampleName options (-t <number of threads> -minlen <min length to keep> -hdclen <head cut off length> -q <quality cutoff> -w <window size> -p <phred offset 64|33>\n";

print STDERR "illumina_qc_threading -s1 $file1 -s2 $file2 -sn $sampleName -t $threads -minlen $min_len -hdclen $headcrop_len -q $qcutoff -w $window_size -p $phred\n";

my $bin_dir = "/export/home/xdong/bin";

#/data/metatools/bowtie2-2.0.0-beta3/bowtie2-build artifacts.fasta artifacts_bowtie
#/data/metatools/bowtie2-2.0.0-beta3/bowtie2  --phred64  -q  -t -p 10  -x ../../../../database/artifacts_bowtie -U s_6_1_5_sequence.txt -S test.sam

# filter out artifacts
my $artifacts_db_bowtie =  "/export/hmp_tmp/data/mdata/database/artifacts_bowtie ";
my $bowtie = "/data/metatools/bowtie2-2.0.0-beta3/bowtie2 "; 

my $t0 = Benchmark->new;


my $cmd1 = "";
my $cmd2 = "";
my $cmd = "";
my $p = "";
    
if(! -e "$file1.sam" || ! -e "$file2.sam"){
    
    $cmd1 = "$bowtie  --phred$phred  -q  -t -p $threads  -x  $artifacts_db_bowtie -U $file1 -S $file1.sam; ";
    $cmd2 = "$bowtie  --phred$phred  -q  -t -p $threads  -x  $artifacts_db_bowtie -U $file2 -S $file2.sam; ";
}



if(! -e "$file1.sam.filtered" || ! -e "$file2.sam.filtered"){
    
    $cmd1 .= "perl  -ne \'next if(\/^\\S+\\s+?4\\s+\/); print \$_; \'  $file1.sam > $file1.sam.filtered; ";
    $cmd2 .= "perl  -ne \'next if(\/^\\S+\\s+?4\\s+\/); print \$_;\'  $file2.sam > $file2.sam.filtered; ";
}

if(length($cmd1) && length($cmd2)){
    my $thr1 = threads->new(\&thrsub, $cmd1);
    my $thr2 = threads->new(\&thrsub, $cmd2);
    $thr1->join();
    $thr2->join();
}
$cmd1 = "";
$cmd2 = "";


#fiter out the artifacts accroding to the mapping results, pair sequecnes and do head hard cropping
if(! -e "$file1.artiF.txt" || ! -e "$file2.artiF.txt"){
    my $p = "perl $bin_dir/getUnalignedFastqPair.pl ";
    $cmd = "$p $file1 $file2 $file1.sam.filtered $file2.sam.filtered $headcrop_len";
    print STDERR "$cmd\n";
    system $cmd;
    
}

$cmd = "";
####################################### cut off adaptor ################################################


my $t_filter = Benchmark->new;
my $td = timediff($t_filter, $t0);
print STDERR "Illumina QC filter takes :",timestr($td)," to run\n";

#s_6_1_5_sequence.txt.artiF.atrimmed.txt
if(! -e "$file1.artiF.atrimmed.txt" || ! -e "$file2.artiF.atrimmed.txt"){
    $p = "/data/metatools/cutadapt-1.0/cutadapt";
    $cmd1 = "$p -f fastq -a AGATCGGAAGAGC -o 10  --quality-base=$phred  -o $file1.artiF.atrimmed.txt  $file1.artiF.txt; ";
    $cmd2 = "$p -f fastq -a AGATCGGAAGAGC -o 10  --quality-base=$phred  -o $file2.artiF.atrimmed.txt  $file2.artiF.txt; ";
}

if(length($cmd1) && length($cmd2)){
    my $thr1 = threads->new(\&thrsub, $cmd1);
    my $thr2 = threads->new(\&thrsub, $cmd2);
    $thr1->join();
    $thr2->join();
}

$cmd1 = "";
$cmd2 = "";

my $t_trimming = Benchmark->new;
$td = timediff($t_trimming, $t_filter);
print STDERR "Illumina QC adaptor trimming takes :",timestr($td)," to run\n";

$cmd = "";
#filter out the shortmer and shuffle the sequences
if(! -e "$sampleName.artiF.noa.shuffled.fastq"){
    $p = "perl $bin_dir/filterUnpairedReads.pl";
    $cmd = "$p $file1.artiF.atrimmed.txt  $file2.artiF.atrimmed.txt  $min_len > $sampleName.artiF.noa.shuffled.fastq; ";
    print STDERR "$cmd\n";
    system $cmd;
}
if(! -e "$sampleName.artiF.noa.shuffled.qtrim.fastq"){

    $cmd = "perl $bin_dir/fastq_qualitytrim_window.pl  -q $qcutoff -m $min_len -w $window_size -o $phred -p 1 $sampleName.artiF.noa.shuffled.fastq > $sampleName.artiF.noa.shuffled.qtrim.fastq";
    print STDERR "$cmd\n";
    system $cmd;
}

if(! -e "$sampleName.artiF.noa.shuffled.qtrim.dust.fastq"){
    if($phred == 64){
	$cmd = "/data/assemblers/sga0.9.17/bin/sga  preprocess -v -o $sampleName.artiF.noa.shuffled.qtrim.dust.fastq -m $min_len --pe-mode=2 --phred64 --dust --dust-threshold=4  $sampleName.artiF.noa.shuffled.qtrim.fastq >& $sampleName.artiF.noa.shuffled.qtrim.lowcomplexity.fastq";
    }
    else{
	
	$cmd = "/data/assemblers/sga0.9.17/bin/sga  preprocess -v -o $sampleName.artiF.noa.shuffled.qtrim.dust.fastq -m $min_len --pe-mode=2  --dust --dust-threshold=4  $sampleName.artiF.noa.shuffled.qtrim.fastq >& $sampleName.artiF.noa.shuffled.qtrim.lowcomplexity.fastq";
    }
    print STDERR "$cmd\n";
    system $cmd;
}
#if(! -e "$sampleName.artiF.noa.shuffled.qtrim.dust.rf.fastq"){
 #   $cmd = "perl $bin_dir/catRF_singlefile.pl $sampleName.artiF.noa.shuffled.qtrim.dust.fastq > $sampleName.artiF.noa.shuffled.qtrim.dust.rf.fastq; ";
  #  print STDERR "$cmd\n";
   # system $cmd;
    
#}

#if(! -e "collapsed.txt"){
 #   $cmd = " /data/metatools/brentp/fastq_removedRecord  filter --adjust 33 --unique  $sampleName.artiF.noa.shuffled.qtrim.dust.rf.fastq > collapsed.txt; ";
  #  print STDERR "$cmd\n";
   # system $cmd;
    
#}
#if(! -e "$sampleName.artiF.noa.shuffled.qtrim.dust.unique.fastq"){
 #   $cmd = "perl $bin_dir/getUniqueSeq.pl  $sampleName.artiF.noa.shuffled.qtrim.dust.fastq collapsed.txt  > $sampleName.artiF.noa.shuffled.qtrim.dust.unique.fastq";
  #  print STDERR "$cmd\n";
   # system $cmd;
#}


#$cmd = "perl $bin_dir/unshuffleSequences_fastx.pl  4 $sampleName.artiF.noa.shuffled.qtrim.dust.unique.fastq $sampleName.qc";
#system $cmd;

my $t_dust_unique = Benchmark->new;
$td = timediff($t_dust_unique, $t_trimming);
print STDERR "Illumina QC dust and unique takes :",timestr($td)," to run\n";


my $t1 = Benchmark->new;
$td = timediff($t1, $t0);
print STDERR "Illumina QC took:",timestr($td)," in total to finish\n";

sub thrsub {
    my ($cmd) = @_;
    print STDERR "Threading: $cmd\n";
    system $cmd;
    
    print STDERR "************ Threading: $cmd ended\n";
}
