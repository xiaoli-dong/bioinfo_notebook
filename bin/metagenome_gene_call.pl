#!/usr/local/bin/perl

use strict;
use warnings;
use threads;
use Benchmark;

@ARGV == 2
  or die "Usage: $0 <fasta seq file><sample name or sequence prefix>\n";

my $contig = $ARGV[0];

my $last_period_index = rindex( $contig, "." );
my $prefix            = substr( $contig, 0, $last_period_index );
my $seqprefix         = $ARGV[1];
my $t0                = Benchmark->new;
my $bin               = "/export/home/xdong/bin";

my $cmd      = "";
my $tRNA     = "/export/data/programs/tRNAscan-SE-1.4alpha/bin/tRNAscan-SE";
my $tRNAMask = "perl $bin/mask_tRNA.pl";
my $rna =
"/export/data/programs/meta_rna/rRNA_prediction/scripts/rRNA_hmm_run_wst_v0.pl";
my $gene_predict = "/export/data/programs/MetaGeneMark_linux64/gmhmmp";
my $m = "/export/data/programs/MetaGeneMark_linux64/MetaGeneMark_v1.mod";
my $process_gene_predict = "perl $bin/metagenemark_parser.pl";
my $rnagff               = "perl $bin/extractRNAgff.pl";
my $rnafasta             = "perl $bin/extractRNAFasta.pl";
my $getGff               = "perl $bin/getGff.pl";
my $gff2fasta            = "perl $bin/gff2fasta.pl";

################### tRNA prediction and mask ##########################################

if ( !-e "$contig.tRNAscan" ) {
    $cmd =
"$tRNA -G -f $contig.structure -o $contig.tRNAscan -m tRNAscan.stat.txt $contig; ";
    print STDERR "***************************** $cmd\n\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}
my $t_trna = Benchmark->new;
my $td     = timediff( $t_trna, $t0 );
print STDERR "tRNAScan-SE takes :", timestr($td), " to run\n";

$cmd = "";
if ( !-e "$prefix.tRNAMask.fna" ) {
    $cmd = "$tRNAMask  $contig $contig.tRNAscan >  $prefix.tRNAMask.fna; ";
    print STDERR "***************************** $cmd\n\n";

  #bitshift the return value by 8 to get the return value of the program called.
    ( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";
}

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

##################protein coding gene prediciton ################################

if ( !-e "$contig.predicted.genes.ffn" || !-e "$contig.predicted.genes.ffa" ) {

    $cmd =
"$gene_predict -a -d -f G -m $m -o $prefix.tRNAMask.rRNAMask.metagene.gff $prefix.tRNAMask.fna.mask; ";
    $cmd .=
"$process_gene_predict $prefix.tRNAMask.rRNAMask.metagene.gff $prefix.predicted.genes $seqprefix; ";

    print STDERR "$cmd\n\n";
    ( system $cmd) >> 8 and die "Cound not execute cmd=$cmd, $!\n";
}

my $t_metagene = Benchmark->new;
$td = timediff( $t_metagene, $t_rna );
print STDERR "metagenemark takes :", timestr($td), " to run\n";

##################Get GFF Gile ################################
$cmd =
"$getGff -s ../$contig -t ../$contig.tRNAscan -r $prefix.tRNAMask.fna.coord -c $prefix.tRNAMask.rRNAMask.metagene.gff > $seqprefix.gff;";
print STDERR "$cmd\n\n";
( system $cmd) >> 8 and die "Cound not execute cmd=$cmd, $!\n";

############################# fasta file ################################

$cmd =
  "$gff2fasta -s ../$contig -g $seqprefix.gff -f tRNA > $prefix.tRNA.fasta;";
$cmd .=
"$gff2fasta -s ../$contig -g $seqprefix.gff -f 5SrRNA > $prefix.5SrRNA.fasta;";
$cmd .=
"$gff2fasta -s ../$contig -g $seqprefix.gff -f 16SrRNA -l 180 > $prefix.16SrRNA.fasta;";
$cmd .=
"$gff2fasta -s ../$contig -g $seqprefix.gff -f 23SrRNA -l 180 > $prefix.23SrRNA.fasta;";
$cmd .=
"$gff2fasta -s ../$contig -g $seqprefix.gff -f CDS -l 180 > $prefix.orf.fasta;\n";

print STDERR "$cmd\n\n";
( system $cmd) >> 8 and die "Cound not execute cmd=$cmd, $!\n";

#if(! -e "5SrRNA.fasta" || ! -e "5SrRNA.fasta" || ! -e "5SrRNA.fasta"){

#   $cmd = "$rnagff $prefix.tRNAMask.fna.coord; ";
#  $cmd .= "$rnafasta $prefix.tRNAMask.fna.seq $seqtype";
# print STDERR "$cmd\n\n";
#(system $cmd) >> 8 and die "Cound not execute cmd=$cmd, $!\n";
#}

print STDERR "metagenome gene predicted finished\n";

my $t_end = Benchmark->new;
$td = timediff( $t_end, $t0 );
print STDERR "metagenome gene prediction takes :", timestr($td), " to run\n";

