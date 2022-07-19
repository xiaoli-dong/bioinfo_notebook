#!/usr/bin/perl

@ARGV == 3
  or die
"Usage: $0 <fasta file><case sensitive: 1 sensitive, 0 no. when it's 1, low case means the masked bases><length cutoff>\n",

"Calculate the z-score of the tetrnucletide frequency: assume the input is big case\n";

my $fasta      = $ARGV[0];
my $sen        = $ARGV[1];
my $len_cutoff = 0;

chomp($fasta);
open( FASTA,  "$fasta" );
open( FREQ,   ">$fasta.freq" );
open( ZSCORE, ">$fasta.zscore" );

#count sequencs it has in file
my @header = ();

#	print STDERR $theseq, "\n";
my @n1a = ( 'A', 'C', 'G', 'T' );
my @n2a = ( 'A', 'C', 'G', 'T' );
my @n3a = ( 'A', 'C', 'G', 'T' );
my @n4a = ( 'A', 'C', 'G', 'T' );
my $i   = 0;
foreach my $n1 (@n1a) {
    foreach my $n2 (@n2a) {
        foreach my $n3 (@n3a) {
            foreach my $n4 (@n4a) {
                $header[$i] = "$n1$n2$n3$n4";
                $i++;
            }
        }
    }
}

print FREQ "id";
print ZSCORE "id";

foreach (@header) {

    print FREQ "\t$_";
    print ZSCORE "\t$_";
}
print FREQ "\n";
print ZSCORE "\n";

$/ = "\n>";
while (<FASTA>) {
    chomp;
    if ( !/>?(\S+?)\s.*?\n(.*)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    my $name = $1;
    my $seq  = $2;

    # print STDERR "$name\n";
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    if ( length($seq) < $len_cutoff ) {

        next;
    }

    if ( $sen == 0 ) {
        $seq = uc($seq);
    }

    my $zscore = compute_zscores_for_tetranucleotides( $name, $seq );
    my $freqs  = compute_tetranucleotides_freq( $name, $seq );

    print FREQ "$name";

    foreach my $freq (@$freqs) {
        print FREQ "\t$freq";
    }
    print FREQ "\n";

    print ZSCORE "$name";
    foreach my $score (@$zscore) {
        print ZSCORE "\t$score";
    }
    print ZSCORE "\n";

}

close(FREQ);
close(ZSCORE);

#######################################################################################

sub RevComp {
    my ($seq) = @_;
    $reverse_seq = reverse($seq);
    $reverse_seq =~ tr/ATGC/TACG/;
    return $reverse_seq;
}

#######################################################################################
sub compute_zscores_for_tetranucleotides {
    my ( $name, $theseq ) = @_;

    #$theseq = uc($theseq);

    #as described by Teeling et al. BMC Bioinformatics 2004, 5:163.
    $theseq .= RevComp($theseq);

    my %oligos2 = ();
    my %oligos3 = ();
    my %oligos4 = ();
    my @zscore  = ();
    my @header  = ();
    my $i       = 0;
    my $len     = length($theseq) - 2 + 1;

    while ( $i < $len ) {
        my $seq = substr( $theseq, $i, 2 );
        $oligos2{$seq}++;
        $i++;
    }

    $i   = 0;
    $len = length($theseq) - 3 + 1;

    while ( $i < $len ) {
        $seq = substr( $theseq, $i, 3 );
        $oligos3{$seq}++;
        $i++;
    }

    $i = 0;

    $len = length($theseq) - 4 + 1;
    while ( $i < $len ) {
        $seq = substr( $theseq, $i, 4 );
        $oligos4{"$seq"}++;
        $i++;
    }

    foreach ( sort keys %oligos4 ) {

        #print STDERR "$_\t", length($_), "\t$oligos4{$_}\n";
    }

    $i = 0;
    foreach my $n1 (@n1a) {
        foreach my $n2 (@n2a) {
            foreach my $n3 (@n3a) {
                foreach my $n4 (@n4a) {
                    if ( not exists $oligos4{"$n1$n2$n3$n4"} ) {
                        $zscore[$i] = 0;
                    }
                    else {
                        $exp{"$n1$n2$n3$n4"} =
                          ( $oligos3{"$n1$n2$n3"} * $oligos3{"$n2$n3$n4"} ) /
                          $oligos2{"$n2$n3"};

                     #print STDERR "$n1$n2$n3$n4\t", $exp{"$n1$n2$n3$n4"}, "\n";
                        $var{"$n1$n2$n3$n4"} = $exp{"$n1$n2$n3$n4"} * (
                            (
                                ( $oligos2{"$n2$n3"} - $oligos3{"$n1$n2$n3"} )
                                * (
                                    $oligos2{"$n2$n3"} - $oligos3{"$n2$n3$n4"}
                                )
                            ) / ( $oligos2{"$n2$n3"}**2 )
                        );

                        #print STDERR $var{"$n1$n2$n3$n4"}, "\n";
                        if ( sqrt( $var{"$n1$n2$n3$n4"} ) == 0 ) {
                            print STDERR "$name\t$n1$n2$n3$n4,",
                              $oligos4{"$n1$n2$n3$n4"}, ",",
                              $var{"$n1$n2$n3$n4"}, ",", $exp{"$n1$n2$n3$n4"},
                              ",", $oligos2{"$n2$n3"}, ",",
                              $oligos3{"$n1$n2$n3"}, ",",
                              $oligos3{"$n2$n3$n4"}, "\n";
                        }
                        $zscore[$i] =
                          ( $oligos4{"$n1$n2$n3$n4"} - $exp{"$n1$n2$n3$n4"} ) /
                          sqrt( $var{"$n1$n2$n3$n4"} );

                        #print STDERR $zscore{"$n1$n2$n3$n4"}, "\n";
                    }
                    $i++;
                }

            }

        }

    }

    return \@zscore;
}

#######################################################################################
sub compute_tetranucleotides_freq {
    my ( $name, $theseq ) = @_;

    #$theseq = uc($theseq);

    #as described by Teeling et al. BMC Bioinformatics 2004, 5:163.
    #$theseq .=RevComp($theseq);

    my $seqLen = length($theseq);

    my %oligos4 = ();
    my @freqs   = ();

    my $i = 0;

    my $len = length($theseq) - 4 + 1;
    while ( $i < $len ) {
        $seq = substr( $theseq, $i, 4 );
        $oligos4{"$seq"}++;
        $i++;
    }

    $revseq = RevComp($theseq);
    $i      = 0;
    while ( $i < $len ) {
        $seq = substr( $revseq, $i, 4 );
        $oligos4{"$seq"}++;
        $i++;
    }

    my $total = 0;
    foreach ( keys %oligos4 ) {

        $total += $oligos4{$_};
    }

    $i = 0;
    foreach my $n1 (@n1a) {
        foreach my $n2 (@n2a) {
            foreach my $n3 (@n3a) {
                foreach my $n4 (@n4a) {
                    if ( not exists $oligos4{"$n1$n2$n3$n4"} ) {
                        $freqs[$i] = 0;

                    }
                    else {

                        $freqs[$i] = sprintf( "%.4f",
                            256 * ( $oligos4{"$n1$n2$n3$n4"} / $total ) );

                    }
                    $i++;
                }

            }

        }

    }
    return \@freqs;
}
