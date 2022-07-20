#!/usr/local/bin/perl

use strict;
use warnings;
use threads;
use Benchmark;
use Getopt::Long;

my $threads = 10;

my ( $ref, $rid, $s1, $s2, $sampleName );
&GetOptions(

    "r=s"   => \$ref,
    "rid=s" => \$rid,
    "s1=s"  => \$s1,
    "s2=s"  => \$s2,
    "sn=s"  => \$sampleName,
    "t=i"   => \$threads
);

( $ref && $rid && $s1 && $s2 && $sampleName ) || die "usage: $0 OPTIONS

where options are:\n-r <fasta format reference contigs>\n-rid <reference index name>\n-s1  seqfile1\n-s2 seqfile2\n-sn sampleName\noptions (-t <number of threads>)\n";

print STDERR
"illumina_coverage -r $ref -rid $rid -s1 $s1 -s2 $s2 -sn $sampleName -t $threads\n";
my $t0 = Benchmark->new;

genomeFileForGenomeCoverageBed( $ref, $sampleName );
getCoverageHistogram( $ref, $s1, $s2, $sampleName, $threads );

my ( $covs, $lens ) = getContigCoverage("$sampleName.coverage.hist.txt");
my $gcs = getGC($ref);

open( PRO, ">$ref.signatures.txt" )
  or die "could not open $ref.signatures.txt to write, $!\n";

print PRO "id\tlength\tgc\tcoverage\n";
foreach my $contigid ( sort { $lens->{$b} <=> $lens->{$a} } keys %$lens ) {
    print PRO "$contigid\t", $lens->{$contigid}, "\t", $gcs->{$contigid}, "\t",
      $covs->{$contigid}, "\n";
}
close(PRO);

my $t1 = Benchmark->new;

my $td = timediff( $t1, $t0 );

print STDERR "illumina_coverage takes ", timestr($td), "to run\n";

sub getCoverageHistogram {
    my ( $ref, $s1, $s2, $sampleName, $threads ) = @_;
    my $bowtie_home  = "/data/metatools/bowtie2-2.1.0";
    my $samtool_home = "/usr/local/bin";
    my $bedtool_home = "/data/metatools/bedtools-2.17.0/bin";
    my $cmd          = "";
    if ( !-e "$rid.1.bt2" ) {
        $cmd .= "$bowtie_home/bowtie2-build $ref $rid;";
    }
    if ( !-e "$sampleName.sam" ) {

        $cmd .=
"$bowtie_home/bowtie2 -q -t -p $threads -x $rid -1 $s1 -2 $s2 -S $sampleName.sam;";
    }
    if ( !-e "$sampleName.sorted.bam" ) {
        $cmd .=
"$samtool_home/samtools view -bS $sampleName.sam |  $samtool_home/samtools sort - $sampleName.sorted;";
    }
    if ( !-e "$sampleName.coverage.hist.txt" ) {

        $cmd .=
"$bedtool_home/genomeCoverageBed -ibam $sampleName.sorted.bam  -g $sampleName.genome.txt > $sampleName.coverage.hist.txt;";
    }

#my $cmd = "$bowtie_home/bowtie2-build $ref $rid;$bowtie_home/bowtie2 -q -t -p $threads -x $rid -1 $s1 -2 $s2 -S $rid.sam; $samtool_home/samtools view -bS $rid.sam |  $samtool_home/samtools sort - $rid.sorted; $bedtool_home/genomeCoverageBed -d -ibam $rid.sorted.bam  -g $rid.genome.txt > $sampleName.coverage.hist.txt";

    print STDERR join( "\n", split( ";", $cmd ) );
    system $cmd;
}

sub getContigCoverage {

    my ($covHistFile) = @_;

    my %covs    = ();
    my %lengths = ();

    open( FILE, $covHistFile )
      or die "Could not open $covHistFile to read, $!\n";

#chromosome (or entire genome)
#depth of coverage from features in input file
#number of bases on chromosome (or genome) with depth equal to column 2.
#size of chromosome (or entire genome) in base pairs
#fraction of bases on chromosome (or entire genome) with depth equal to column 2.

    #contig1	3	1	550011	1.81815e-06
    #contig1	4	1	550011	1.81815e-06

    my $totalBases  = 0;
    my $pre_id      = "";
    my $pre_length  = 0;
    my $curr_id     = "";
    my $curr_length = 0;

    my $cov = 0;
    while (<FILE>) {

        chomp;
        my @line = split( "\t", $_ );
        $curr_id     = $line[0];
        $curr_length = $line[3];

        if ( !length($pre_id) ) {

            #print STDERR "pre_id=$pre_id, curr_id=$curr_id\n";
            $pre_id     = $line[0];
            $totalBases = $line[1] * $line[2];
            $pre_length = $curr_length;

        }
        elsif ( length($pre_id) && $pre_id ne $curr_id ) {

     #print STDERR "pre_id=$pre_id, curr_id=$curr_id, pre_length=$pre_length\n";
            $covs{$pre_id}    = sprintf( "%.0f", $totalBases / $pre_length );
            $lengths{$pre_id} = $pre_length;

            $pre_id     = $curr_id;
            $pre_length = $curr_length;
            $totalBases = $line[1] * $line[2];
        }

        else {
            $totalBases += $line[1] * $line[2];
        }

    }

    $covs{$pre_id}    = sprintf( "%.0f", $totalBases / $pre_length );
    $lengths{$pre_id} = $pre_length;

    return ( \%covs, \%lengths );
}

sub genomeFileForGenomeCoverageBed {
    my ( $fasta, $sampleName ) = @_;
    open( FASTA,  "$fasta" ) or die "Could not open $fasta to read, $!\n";
    open( GENOME, ">$sampleName.genome.txt" )
      or die "Could not open $sampleName.genome.txt to write, $!\n";

    $/ = "\n>";
    while (<FASTA>) {
        chomp;
        if ( my ( $seqName, $seq ) = /^>?(\S+).*?\n(.*)/s ) {
            $seq =~ s/\s//g;
            print GENOME "$seqName\t", length($seq), "\n";
        }
    }

    $/ = "\n";

    close(GENOME);

    close(FASTA);

}

sub getGC {

    my ($fasta) = @_;
    my %gcs = ();

    open( IN, "$fasta" ) or die "Could not open $fasta to read, $!\n";

    $/ = "\n>";
    while (<IN>) {
        chomp;
        if ( my ( $seqName, $DNA ) = /^>?(\S+).*?\n(.*)/s ) {
            $DNA =~ s/\s//g;
            my $A = ( $DNA =~ tr/A// );
            my $C = ( $DNA =~ tr/C// );
            my $G = ( $DNA =~ tr/G// );
            my $T = ( $DNA =~ tr/T// );

            my $gccontent =
              sprintf( "%.2f", 100 * ( $G + $C ) / ( $A + $G + $C + $T ) );

            #print STDERR "$seqName\t$gccontent\n";
            $gcs{$seqName} = $gccontent;

        }
    }

    $/ = "\n";
    close(IN);

    return \%gcs;
}

