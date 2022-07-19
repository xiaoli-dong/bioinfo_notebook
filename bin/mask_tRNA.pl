#!/usr/local/bin/perl

# mask the tRNA genes from input sequences
$#ARGV == -1

  and die "Usage: $0 fasta_seq tRNAscan_output\n";

print STDERR "running maks_tRNA script\n";
open( FEATURES, "$ARGV[1]" ) or die $!;

my %features = ();

while (<FEATURES>) {

    if (/^(Sequence|Name|\---)/) {
        print STDERR "$_";
        next;
    }
    else {
        my @feature = split( /\s+/, $_ );
        $features{ $feature[0] }->{ $feature[1] } = \@feature;
    }
}

close(FEATURES);

open( SEQ, "$ARGV[0]" ) or die $!;

$/ = "\n>";
my $seq2quals = ();

while (<SEQ>) {

    chomp;
    if ( !/^>?(\w+).*?\n(.+)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    my $seqid = $1;
    my $seq   = $2;

    $seq =~ tr/ \r\n\t//d;

    my @chars = split( //, $seq );
    if ( exists $features{$seqid} ) {

        foreach my $id ( keys %{ $features{$seqid} } ) {
            my $start = $features{$seqid}->{$id}[2];
            my $end   = $features{$seqid}->{$id}[3];

            if ( $start < $end ) {
                for ( my $index = $start - 1 ; $index < $end ; $index++ ) {
                    $chars[$index] = "N";
                }
            }
            else {

                for ( my $index = $end - 1 ; $index < $start ; $index++ ) {
                    $chars[$index] = "N";
                }

            }
        }

    }
    print ">$seqid\n", join( '', @chars ), "\n";
}

close(SEQ);
