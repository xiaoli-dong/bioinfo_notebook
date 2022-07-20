#!/usr/local/bin/perl

use strict;
use warnings;

my $qfile = $ARGV[0];
my $sfile = $ARGV[1];

open( QFILE, $qfile );
open( SFILE, $sfile );

while (<QFILE>) {
    open( QFW, ">qfile.fasta" );
    open( SFW, ">sfile.fasta" );
    my $qhead = $_;
    my $qseq  = <QFILE>;
    <QFILE>;
    <QFILE>;
    print QFW ">$qhead", $qseq;

    my $shead = <SFILE>;
    my $sseq  = <SFILE>;
    <SFILE>;
    <SFILE>;
    print SFW ">$shead", $sseq;
    close(SFW);
    close(QFW);
    my $cmd =
"blastn -query qfile.fasta -subject sfile.fasta -outfmt \"6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore\"  >> blast2.results.txt";
    system $cmd;

}
