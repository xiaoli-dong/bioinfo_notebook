#!/usr/local/bin/perl
use warnings;
use strict;

$#ARGV == 2

  or die
"Usage: $0 <gff format output from metagenemark> <outputname> <seqid prefix>\n";

my $prefix = $ARGV[2];

open( GFF, "$ARGV[0]" )      or die "Could not open $ARGV[0] to read, $!";
open( DNA, ">$ARGV[1].ffn" ) or die "Could not open ARGV[1].ffn to read, $!";
open( PROTEIN, ">$ARGV[1].ffa" )
  or die "Could not open ARGV[1].ffa to read, $!";
open( DNAFORMATED, ">$ARGV[1].org.ffn" )
  or die "Could not open ARGV[1].org.ffn to write, $!";
open( PROTEINFORMATED, ">$ARGV[1].org.ffa" )
  or die "Could not open ARGV[1].org.ffa to write, $!";
my $seq           = "";
my $id            = 0;
my $type          = 0;
my $dna_index     = 0;
my $protein_index = 0;
my $contig        = "";

my @starts  = ();
my @ends    = ();
my @strands = ();
my @ids     = ();

my $start = 0;

while (<GFF>) {

    if ( !/^(#)/ ) {
        $start = 1;
    }
    next until ($start);

    if (/^\s+/) {
        $protein_index = 0;
        $dna_index     = 0;
        $contig        = "";
        @starts        = ();
        @ends          = ();
        @strands       = ();
        @ids           = ();
        $type          = 0;
    }
    elsif (/^(\w+)\s+.*?(\d+)\s+(\d+)\s+\S+\s+(\S+).*?gene_id\s+(\d+)$/) {

        $contig = $1;
        push( @starts, $2 );
        push( @ends,   $3 );

        #if($4 eq "-"){
        push( @strands, "-" );

        #}
        #elsif($4 eq "+"){
        push( @strands, "+" );

        #}

        push( @ids, $5 );
    }
    elsif (/^##Protein\s+(\d+)/) {

        $id   = $1;
        $type = 1;
    }
    elsif (/^##DNA\s+(\d+)/) {
        $id   = $1;
        $type = 2;
    }
    elsif (/^##end/) {

        if ( $type == 1 ) {
            if ( length($seq) >= 60 ) {

                print PROTEINFORMATED
">$prefix\_$contig\_$starts[$protein_index]\_$ends[$protein_index] /start=$starts[$dna_index] /end=$ends[$dna_index] strand=$strands[$protein_index] id=$id len=",
                  length($seq), "\n$seq\n";
                print PROTEIN ">gene_id\_$id\n$seq\n";
                $protein_index++;
            }
            else {
                $protein_index++;
            }
        }
        elsif ( $type == 2 ) {
            if ( length($seq) >= 180 ) {
                print DNAFORMATED
">$prefix\_$contig\_$starts[$dna_index]\_$ends[$dna_index] /start=$starts[$dna_index] /end=$ends[$dna_index] strand=$strands[$dna_index] id=$id len=",
                  length($seq), "\n$seq\n";
                print DNA ">gene_id\_$id\n$seq\n";
                $dna_index++;
            }
            else {
                $dna_index++;
            }
        }

        $seq = "";
    }
    elsif (/^##(.*)/) {

        $seq .= $1;
    }

}

close(GFF);
close(DNA);
close(PROTEIN);
