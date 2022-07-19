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

my $adapter_end   = "CTGTCTCTTATACA";
my $adapter_front = "AGATGTGTATAAGAGACAG";

if ( $lib eq "truseq" ) {
    $adapter_end   = "gatcggaagagc";
    $adapter_front = "gctcttccgatct";
}

#######Trim NN at the end#################

if ( !-e "$sampleName.trimN.1.fastq" || !-e "$sampleName.trimN.2.fastq" ) {
    $cmd =
"bbduk.sh -Xmx1g in=$file1 in2=$file2 out=$sampleName.trimN.1.fastq  out2=$sampleName.trimN.2.fastq qtrim=rl trimq=1 minlength=$minlen interleaved=true overwrite=true;";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

#trim the adapter
if ( !-e "$sampleName.trimA.1.fastq" || !-e "$sampleName.trimA.2.fastq" ) {
    $cmd =
"cutadapt -f fastq -a $adapter_end  -g $adapter_front -n 2 -u -$tcut -u $hcut -m $minlen -o tmp.1.fastq -p tmp.2.fastq $sampleName.trimN.1.fastq $sampleName.trimN.2.fastq;";
    $cmd .=
"cutadapt -f fastq -a $adapter_end -g $adapter_front -n 2 -u -$tcut -u $hcut -m $minlen -o $sampleName.trimA.2.fastq -p $sampleName.trimA.1.fastq tmp.2.fastq tmp.1.fastq";
    print STDERR "$cmd\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
    unlink "tmp.1.fastq", "tmp.2.fastq"
      or warn "Could not unlink tmp.1.fastq tmp.2.fastq : $!\n";
}

############### do the merge #######################

if (   !-e "=$sampleName.unmerged.1.fastq"
    || !-e "$sampleName.unmerged.2.fastq"
    || !-e "$sampleName.merged.fastq" )
{
    $cmd =
"bbmerge.sh in=$sampleName.trimA.1.fastq in2=$sampleName.trimA.2.fastq out=$sampleName.merged.fastq  outu=$sampleName.unmerged.1.fastq outu2=$sampleName.unmerged.2.fastq outinsert=$sampleName.outinsert.txt qtrim=rl  trimq=5 interleaved=true overwrite=true";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

#######quality trim and filter out the reads with NN
if (   !-e "$sampleName.trimA.qtrim.1.fastq"
    || !-e "$sampleName.trimA.qtrim.2.fastq" )
{
    $cmd =
"bbduk.sh -Xmx1g in=$sampleName.unmerged.1.fastq in2=$sampleName.unmerged.2.fastq out=$sampleName.trimA.qtrim.1.fastq  out2=$sampleName.trimA.qtrim.2.fastq qtrim=rl trimq=$qcutoff minlength=$minlen interleaved=true overwrite=true;";
    $cmd .=
"bbduk.sh -Xmx1g in=$sampleName.merged.fastq out=$sampleName.merged.qtrim.fastq qtrim=rl trimq=$qcutoff minlength=$minlen overwrite=true;";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

########filter out the contaminant#####################

if (   !-e "$sampleName.trimA.qtrim.fc.1.fastq"
    || !-e "$sampleName.trimA.qtrim.fc.2.fastq"
    || !-e "$sampleName.merged.qtrim.fc.fastq" )
{
    $cmd =
"bbduk.sh -Xmx1g in=$sampleName.trimA.qtrim.1.fastq in2=$sampleName.trimA.qtrim.2.fastq out=$sampleName.trimA.qtrim.fc.1.fastq  out2=$sampleName.trimA.qtrim.fc.2.fastq  overwrite=true outm=$sampleName.contaminant.match.fastq ref=/export/data/programs/ebg/contaminant_list.fasta k=31 hdist=1  maxns=0 interleaved=true;";

    $cmd .=
"bbduk.sh -Xmx1g in=$sampleName.merged.qtrim.fastq out=$sampleName.merged.qtrim.fc.fastq   overwrite=true outm=$sampleName.merged.contaminant.match.fastq ref=/export/data/programs/ebg/contaminant_list.fasta k=31 hdist=1  maxns=0";

    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

##########trim the poly-A poly-T tail###################
if (   !-e "$sampleName.qc.1.fastq"
    || !-e "$sampleName.qc.2.fastq"
    || !-e "$sampleName.merged.qc.fastq" )
{
    $cmd =
"cutadapt -f fastq -a aaaaaaaaaaaa -a tttttttttttt -g aaaaaaaaaaaa -g tttttttttttt -n 4  -m $minlen -o $sampleName.trimA.qtrim.fc.aatt.1.fastq  -p $sampleName.trimA.qtrim.fc.aatt.2.fastq  $sampleName.trimA.qtrim.fc.1.fastq $sampleName.trimA.qtrim.fc.2.fastq;";

    $cmd .=
"cutadapt -f fastq -a aaaaaaaaaaaa -a tttttttttttt  -g aaaaaaaaaaaa -g tttttttttttt  -n 4   -m $minlen -o $sampleName.qc.2.fastq -p $sampleName.qc.1.fastq  $sampleName.trimA.qtrim.fc.aatt.2.fastq  $sampleName.trimA.qtrim.fc.aatt.1.fastq;";
    $cmd .=
"cutadapt -f fastq -a aaaaaaaaaaaa -a tttttttttttt -g aaaaaaaaaaaa -g tttttttttttt -n 4  -m $minlen -o $sampleName.merged.qc.fastq  $sampleName.merged.qtrim.fc.fastq ;";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}
###############get sequence stats #############
if (   !-e "$sampleName.qc.1.fastq.stats.txt"
    || !-e "$sampleName.qc.2.fastq.stats.txt"
    || !-e "$sampleName.merged.qc.stats.txt" )
{
    $cmd =
"$^X /export/home/xdong/bin/seqStats.pl -f fastq -s $sampleName.qc.1.fastq > $sampleName.qc.1.fastq.stats.txt;";
    $cmd .=
"$^X /export/home/xdong/bin/seqStats.pl -f fastq -s $sampleName.qc.2.fastq > $sampleName.qc.2.fastq.stats.txt;";
    $cmd .=
"$^X /export/home/xdong/bin/seqStats.pl -f fastq -s $sampleName.merged.qc.fastq  >  $sampleName.merged.qc.stats.txt";
    print STDERR "$cmd\n";
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}
#############get fastqc report ###################

$cmd =
"/export/data/programs/FastQC/fastqc $sampleName.qc.1.fastq $sampleName.qc.2.fastq $sampleName.merged.qc.fastq;";
print STDERR "$cmd\n";
( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

my $t1 = Benchmark->new;
my $td = timediff( $t1, $t0 );
print STDERR "Illumina QC took:", timestr($td), " in total to finish\n";

