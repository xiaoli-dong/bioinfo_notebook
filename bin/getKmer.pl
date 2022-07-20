#!/usr/bin/perl

my $result = generate_unique_Kmers(4);

my $size = @$result;

foreach my $word ( sort @$result ) {

    print "\t$word";
}

print "\nsize=$size\n";

sub revComp {
    my ($seq) = @_;
    $reverse_seq = reverse($seq);
    $reverse_seq =~ tr/ATGC/TACG/;
    return $reverse_seq;
}

#136 unique tetranucletide
sub generate_unique_Kmers {
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

    my @unique_words = ();
    foreach my $word ( sort @words ) {

        my $find = 0;

        foreach my $unique (@unique_words) {

            if ( $word eq $unique ) {
                $find = 1;
                last;
            }
        }

        foreach my $unique (@unique_words) {

            if ( revComp($word) eq $unique ) {
                $find = 1;
                last;
            }
        }

        if ( $find == 0 ) {
            push( @unique_words, $word );
        }
    }

    return ( \@unique_words );
}
