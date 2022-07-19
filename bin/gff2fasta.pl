#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my ( $seq, $gffFile, $feature, $lencutoff );

&GetOptions(
    "s=s" => \$seq,
    "g=s" => \$gffFile,
    "f=s" => \$feature,
    "l=s" => \$lencutoff
);
if ( !$seq ) {
    die
"Usage: $0 -s <fasta sequences of contigs+singletons> -g <gff file> -f <type> -l <length cutoff for genes>\n";
}

if ( not defined $lencutoff ) {

    $lencutoff = 0;
}
my $verbose = 0;
my %gff     = ();

open( GFF, "$gffFile" ) or die "Could not open $gffFile file to read, $!\n";

while (<GFF>) {
    my ( $seqid, undef, $type, $start, $end, undef, $strand, undef, $attrs ) =
      split( "\t", $_ );
    push @{ $gff{$seqid} }, [ $start, $end, $type, $strand, $attrs ];
}

## Do the fasta
my $seqio = Bio::SeqIO->new(
    -file   => $seq,
    -format => 'fasta'
) or die "double fail\n";

while ( my $sobj = $seqio->next_seq ) {
    my $seqid = $sobj->id;

    unless ( defined( $gff{$seqid} ) ) {

        #warn "no features for $seqid\n";
        next;
    }

    my $seq = $sobj->seq;

    for ( @{ $gff{$seqid} } ) {

        my ( $start, $end, $type, $strand, $attrs ) = @$_;
        warn join( "\t", $start, $end, $type, $attrs ), "\n" if $verbose > 0;

        if ( $feature eq $type ) {
            my %attrs  = split( /=|;|\s+/, $attrs );
            my $seqlen = abs( $end - $start ) + 1;
            if ( $seqlen >= $lencutoff ) {
                print ">", $attrs{"ID"}, " length=",
                  ( abs( $end - $start ) + 1 ), " strand=$strand\n";

                my $range = "";

                if ( $start < $end ) {
                    $range = substr( $seq, $start - 1, $end - $start + 1 );
                }
                else {
                    $range = substr( $seq, $end - 1, $start - $end + 1 );
                }

                if ( $feature =~ m/SrRNA/i && $strand eq "-" ) {
                    print reverse_complement($range), "\n";
                }
                else {
                    print "$range\n";
                }
            }
        }
    }

    #exit;
}

sub reverse_complement {
    my ($dna) = @_;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

