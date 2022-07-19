#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ( $fastq, $fastq2, $sampleName, $threads, $bowtieIndex );
$threads = 10;

&GetOptions(

    "s1=s" => \$fastq,
    "s2=s" => \$fastq2,
    "sn=s" => \$sampleName,
    "t=i"  => \$threads
);

#bowtie2-build -f  LTP106.dna.fasta 16SDB
#bowtie2 -p 20 -x ~/deve/database/living_tree/16SDB -1 GEM1BrNov714.qc.1.fastq -2 GEM1BrNov714.qc.2.fastq -S mapp.sam

my $bowtieIndex = "/export/home/xdong/deve/database/living_tree/16SDB";
my $cmd =
"bowtie2 -p $threads -x $bowtieIndex -1 $fastq -2 $fastq2 -S $sampleName.LTP.sam;";
( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

#filter out the unmapped reads
$cmd =
"samtools view -bS $sampleName.LTP.sam | samtools view -h -F 4 - > $sampelName.LTP.mapped.sam";
( system $cmd) >> 8 and die "Could not execute cmd=$cmd, $!\n";

sub sam2insertSizes {

    open( BOWTIE_OUTPUT, "< $bowtie_output" )
      or die "Can't open $bowtie_output\n";
    my $prev_seq_name;
    my $prev_pos;
    my $prev_base_name;
    my $prev_lane;
    my $prev_dir;
    my $first_bool = 0;
    my $count      = 0;
    my @diffs;

    while (<BOWTIE_OUTPUT>) {

        #skip sam file header
        next if (/^@/);
        my @line = split( /\t/, $_ );
        die "not enough fields in line $_\n" if ( @line < 11 );
        my $seq_name  = @line[0];
        my $pos       = @line[3];
        my @tmp       = split( /\//, $seq_name );
        my $base_name = @tmp[0];
        my $dir       = @line[1];

        #print STDERR "seqname=$base_name pos= $pos, prev_pos=$prev_pos\n";
        if ( $dir eq $prev_dir ) {

#print STDERR "seqname=$base_name prev_base_name=$prev_base_name,pos= $pos, prev_pos=$prev_pos dir=$dir, prev_dir=$prev_dir\n";
#die "$_ read pairs were mapped in the same orientation, check on this\n";
            next;
        }

#get the length of the sequence b/c we are going to add it to the difference since
#bowtie reports the starts of both pairs in reference space, but velvet wants the
#length of the sequenced fragment which includes the reads
        my $sequence   = @line[9];
        my $seq_length = length($sequence);
        if ( $base_name eq $prev_base_name ) {
            $count++;
            my $diff;
            if ( $prev_pos > $pos ) {
                $diff = $prev_pos - $pos + $seq_length;
            }
            elsif ( $pos > $prev_pos ) {
                $diff = $pos - $prev_pos + $seq_length;
            }
            else {
                print STDERR
                  "seqname=$seq_name pos= $pos, prev_pos=$prev_pos\n";

#die "something looks wrong with the positions of the mapped reads, $seq_name\n";
                next;
            }
            print STDERR "$base_name\t$diff\n";
            push( @diffs, $diff );
        }
        $prev_base_name = $base_name;
        $prev_dir       = $dir;

        #$prev_lane = $lane;
        $prev_pos = $pos;

    }

    for my $i (@diffs) {
        print "$i\n";
    }

    #print "@diffs\n";
    #print "$count\n";

}
