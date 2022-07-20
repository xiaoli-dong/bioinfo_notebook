#!/usr/local/bin/perl

use strict;
use warnings;
use Benchmark;

@ARGV == 2 or die "Usage: $0 <fasta seq file> <contig|singleton|raw>\n";

my $contig = $ARGV[0];

my $last_period_index = rindex( $contig, "." );
my $prefix            = substr( $contig, 0, $last_period_index );

my $bin      = "/export/home/xdong/bin";
my $seqtype  = $ARGV[1];
my $t0       = Benchmark->new;
my $cmd      = "";
my $tRNA     = "/data/metatools/tRNAscan-SE-1.23/tRNAscan-SE";
my $tRNAMask = "perl $bin/mask_tRNA.pl";
my $rna =
  "/data/metatools/meta_rna/rRNA_prediction/scripts/rRNA_hmm_run_wst_v0.pl";
my $rnagff   = "perl $bin/extractRNAgff.pl";
my $rnafasta = "perl $bin/extractRNAFasta.pl";

################### tRNA prediction and mask ##########################################

if ( !-e "$contig.tRNAscan" ) {
    $cmd =
"$tRNA -G -f $contig.structure -o $contig.tRNAscan -m tRNAscan.stat.txt $contig; ";
    $cmd .= "$tRNAMask  $contig $contig.tRNAscan >  $prefix.tRNAMask.fna; ";

    print STDERR "***************************** $cmd\n\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

}
my $t_trna = Benchmark->new;
my $td     = timediff( $t_trna, $t0 );
print STDERR "tRNAScan-SE takes :", timestr($td), " to run\n";

############################# rRNA prediction and mask ##########################
$cmd = "";
my $input = "input.1";
mkdir( $input, 0777 ) unless -e "$input";
chdir $input;
symlink( "../$prefix.tRNAMask.fna", "$prefix.tRNAMask.fna" )
  unless -e "$prefix.tRNAMask.fna";
chdir "../";

if ( !-e "output.1/$prefix.tRNAMask.fna.mask" ) {
    $cmd = "$rna $input output.1 ";
    use Cwd;
    my $dir = getcwd;

    print STDERR "$dir\/$cmd\n\n";
    ( system $cmd) >> 8 and die "Cound not execute cmd=$cmd, $!\n";
}

chdir "output.1";

my $t_rna = Benchmark->new;
$td = timediff( $t_rna, $t_trna );
print STDERR "rna prediction takes :", timestr($td), " to run\n";

############################# extract rna gff and fasta file ################################

if ( !-e "5SrRNA.fasta" || !-e "5SrRNA.fasta" || !-e "5SrRNA.fasta" ) {

    $cmd = "$rnagff $prefix.tRNAMask.fna.coord; ";
    $cmd .= "$rnafasta $prefix.tRNAMask.fna.seq $seqtype";
    print STDERR "$cmd\n\n";
    ( system $cmd) >> 8 and die "Cound not execute cmd=$cmd, $!\n";
}

print STDERR "metagenome gene predicted finished\n";

my $t_end = Benchmark->new;
$td = timediff( $t_end, $t0 );
print STDERR "metagenome gene prediction takes :", timestr($td), " to run\n";

