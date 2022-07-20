#!/usr/bin/perl

@ARGV == 2 or die "Usage: $0 <fasta file><length cutoff>\n",

  "Calculate sequences tetrnucletide frequency\n";

my $fasta      = $ARGV[0];
my $len_cutoff = $ARGV[1];

chomp($fasta);
open( FASTA, "$fasta" );
open( FREQ,  ">$fasta.freq" );

my $words = generate_k_mer_nucleotides(4);
print FREQ "id";
foreach (@$words) {
    print FREQ "\t$_";
}
print FREQ "\n";

$/ = "\n>";
while (<FASTA>) {
    chomp;
    if ( !/>?(\S+?)\n(.*)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    my $name = $1;
    my $seq  = $2;

    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    if ( length($seq) < $len_cutoff ) {

        next;
    }

    my $freqs = freq_Kmers_gen_symm( $seq, 4, 1 );

    print FREQ "$name";

    foreach my $key ( sort keys %$freqs ) {
        print FREQ "\t", $freqs->{$key};
    }
    print FREQ "\n";

}

close(FREQ);

sub reverse_complement {
    my ($dna) = @_;
    my $revcom = reverse $$dna;
    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    return \$revcom;
}

sub generate_k_mer_nucleotides {
    my ($k)   = @_;
    my $kmers = ();
    my @bases = ( 'A', 'C', 'G', 'T' );
    my @words = @bases;
    for ( my $i = 1 ; $i < $k ; $i++ ) {
        my @newwords;
        foreach my $w (@words) {
            foreach my $b (@bases) {
                push( @newwords, $w . $b );
            }
        }
        undef @words;
        @words = @newwords;
    }
    return ( \@words );
}

sub freq_Kmers_gen {
    my ( $seq, $k ) = @_;
    my $km_a = generate_Kmers($k);
    my $kmers;
    map { $kmers->{$_} = 0 } @$km_a;
    for ( my $i = 0 ; $i <= length($$seq) - $k ; $i++ ) {
        $kmers->{ substr( $$seq, $i, $k ) }++
          if exists( $kmers->{ substr( $$seq, $i, $k ) } );
    }
    my $total = 0;
    $total += $kmers->{$_} for keys %$kmers;
    $kmers->{$_} /= $total for keys %$kmers;
    return ($kmers);
}

sub freq_Kmers_gen_symm {
    my ( $seq, $k, $symm ) = @_;
    my $km_a = generate_Kmers($k);
    my $kmers;
    map { $kmers->{$_} = 0 } @$km_a;
    for ( my $i = 0 ; $i <= length($$seq) - $k ; $i++ ) {
        $kmers->{ substr( $$seq, $i, $k ) }++
          if exists( $kmers->{ substr( $$seq, $i, $k ) } );
    }
    if ( $symm == 1 ) {
        my $rev = reverse_complement($seq);

        for ( my $i = 0 ; $i <= length($$rev) - $k ; $i++ ) {
            $kmers->{ substr( $$rev, $i, $k ) }++
              if exists( $kmers->{ substr( $$rev, $i, $k ) } );
        }
    }
    my $total = 0;
    $total += $kmers->{$_} for keys %$kmers;
    $kmers->{$_} /= $total for keys %$kmers;
    return ($kmers);
}

