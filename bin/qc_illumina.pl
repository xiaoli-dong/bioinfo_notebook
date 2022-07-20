#!/usr/local/bin/perl

use strict;
use warnings;
use threads;
use Benchmark;
use Getopt::Long;

my $minlen  = 30;
my $tcut    = 0;
my $hcut    = 0;
my $qcutoff = 10;
my $lib     = "nextera";

my ( $file1, $file2, $sampleName );
&GetOptions(

    "s1=s"     => \$file1,
    "s2=s"     => \$file2,
    "sn=s"     => \$sampleName,
    "minlen=i" => \$minlen,
    "hcut=i"   => \$hcut,
    "tcut=i"   => \$tcut,
    "q=i"      => \$qcutoff,
    "lib=s"    => \$lib
);

( $file1 && $file2 && $sampleName )
  || die
"usage: $0 -s1  seqfile1 -s2 seqfile2 -sn sampleName\noptions are:\n-minlen <min length to keep: default 30bp>\n-hcut <remove the number of bases from front the reads: default 0>\n-tcut <remove the number of bases from end of the reads: default 0>\n -q <quality trim cutoff: default 10> -lib <nextera|truseq: default nextera>\n";

print STDERR
"$0 -s1 $file1 -s2 $file2 -sn $sampleName -minlen $minlen -hcut $hcut -tcut $tcut -q $qcutoff -lib $lib\n";

my $t0  = Benchmark->new;
my $cmd = "";

my $adapter_end = "CTGTCTCTTATACA";
$adapter_end = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC";
my $adapter_end_read1 =
"-a CTGTCTCTTATACACATCT -a CCGAGCCCACGAGACNNNNNNNN -a NNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
my $adapter_end_read2 =
"-a CTGTCTCTTATACACATCT -a GACGCTGCCGACGANNNNNNNN -a NNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";

my $adapter_front = "AGATGTGTATAAGAGACAG";

# bbduk.sh in=DLM1BRNov714.trimN.1.fastq in2=DLM1BRNov714.trimN.2.fastq out=bclean.1.fastq out2=bclean.2.fastq ref=/export/data/programs/bbmap/resources/nextera.fa.gz  interleaved=auto  overwrite=true k=16 hdist=1 mink=5 ktrim=r minoverlap=18 tbo=t

if ( $lib eq "truseq" ) {
    $adapter_end   = "gatcggaagagc";
    $adapter_front = "gctcttccgatct";
}
#######Get Fastqc report#################
#$cmd = "/export/data/programs/FastQC/fastqc -o fastqc -t 2 -q -f fastq $file1 $file2";
#print STDERR "$cmd\n";
#(system $cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";

#######Trim NN at the end#################

if ( !-e "$sampleName.trimN.1.fastq" || !-e "$sampleName.trimN.2.fastq" ) {
    $cmd =
"bbduk.sh -Xmx50g in=$file1 in2=$file2 out=$sampleName.trimN.1.fastq  out2=$sampleName.trimN.2.fastq qtrim=rl trimq=1 minlength=$minlen interleaved=true overwrite=true;";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

#trim the adapter
if ( !-e "$sampleName.trimA.1.fastq" || !-e "$sampleName.trimA.2.fastq" ) {

    #trim the adapter at the front
    $cmd =
"bbduk.sh -Xmx50g in=$sampleName.trimN.1.fastq in2=$sampleName.trimN.2.fastq out=bcleanl.1.fastq out2=bcleanl.2.fastq literal=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG  interleaved=auto  overwrite=true k=16 hdist=1 edist=1 hdist2=1 edist2=1 mink=6 ktrim=l tbo=t useshortkmers=t rcomp=f minlen=30;";
    print STDERR "$cmd\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

    $cmd =
"bbduk.sh -Xmx50g in=bcleanl.1.fastq in2=bcleanl.2.fastq out=bcleanr.1.fastq out2=bcleanr.2.fastq ref=/export/data/programs/bbmap/resources/nextera.fa.gz,/export/data/programs/bbmap/resources/truseq.fa.gz interleaved=auto rcomp=t overwrite=true k=16 hdist=1 edist=1 hdist2=1 edist2=1 mink=4 ktrim=r tbo=t useshortkmers=t minlen=30;";
    print STDERR "$cmd\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

    $cmd =
"cutadapt -f fastq -a CTGTCTCTTATA -u -4 -m 30 -o tmp.1.fastq -p tmp.2.fastq bcleanr.1.fastq bcleanr.2.fastq";
    print STDERR "$cmd\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

    $cmd =
"cutadapt -f fastq -a CTGTCTCTTATA -u 4 -m 30 -o $sampleName.trimA.2.fastq  -p $sampleName.trimA.1.fastq tmp.2.fastq tmp.1.fastq;";
    print STDERR "$cmd\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

#######quality trim and filter out the reads has "N"
if (   !-e "$sampleName.trimA.qtrim.1.fastq"
    || !-e "$sampleName.trimA.qtrim.2.fastq" )
{
    $cmd =
"bbduk.sh -Xmx50g in=$sampleName.trimA.1.fastq in2=$sampleName.trimA.2.fastq out=$sampleName.trimA.qtrim.1.fastq  out2=$sampleName.trimA.qtrim.2.fastq qtrim=rl trimq=$qcutoff minlength=$minlen interleaved=true overwrite=true";
    print STDERR "$cmd\n\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

########filter out the contaminant#####################

if (   !-e "$sampleName.trimA.qtrim.fc.1.fastq"
    || !-e "$sampleName.trimA.qtrim.fc.2.fastq" )
{
    $cmd =
"bbduk.sh -Xmx50g in=$sampleName.trimA.qtrim.1.fastq in2=$sampleName.trimA.qtrim.2.fastq out=$sampleName.trimA.qtrim.fc.1.fastq  out2=$sampleName.trimA.qtrim.fc.2.fastq  overwrite=true outm=$sampleName.contaminant.match.fastq ref=/export/data/programs/ebg/contaminant_list.fasta k=31 hdist=1  maxns=2 interleaved=true";
    print STDERR "$cmd\n\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

##########trim the poly-A poly-T tail###################
if ( !-e "$sampleName.qc_1.fastq" || !-e "$sampleName.qc_2.fastq" ) {
    $cmd =
"$^X /export/data/programs/prinseq-lite-0.20.4/prinseq-lite.pl  -fastq $sampleName.trimA.qtrim.fc.1.fastq -fastq2 $sampleName.trimA.qtrim.fc.2.fastq -trim_tail_left 3  -trim_tail_right 3 -out_format 3 -out_good $sampleName.qc -min_len 30;";

    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}
###############get sequence stats #############
if (   !-e "$sampleName.qc_1.fastq.stats.txt"
    || !-e "$sampleName.qc_2.fastq.stats.txt" )
{
    $cmd =
"$^X /export/home/xdong/bin/seqStats.pl -f fastq -s $sampleName.qc_1.fastq > $sampleName.qc_1.fastq.stats.txt &";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

    $cmd =
"$^X /export/home/xdong/bin/seqStats.pl -f fastq -s $sampleName.qc_2.fastq > $sampleName.qc_2.fastq.stats.txt &";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}
#############get fastqc report ###################

$cmd =
"/export/data/programs/FastQC/fastqc -o fastqc -q -t 2 $sampleName.qc_1.fastq $sampleName.qc_2.fastq;";
print STDERR "$cmd\n";
( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

my $t1 = Benchmark->new;
my $td = timediff( $t1, $t0 );
print STDERR "Illumina QC took:", timestr($td), " in total to finish\n";

