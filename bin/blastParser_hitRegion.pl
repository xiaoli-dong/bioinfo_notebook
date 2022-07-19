#!/usr/local/bin/perl

$#ARGV == 1

  or die "Usage: $0 <tabular blast output> <fasta sequence file>\n";

open( TAB, $ARGV[0] ) or die "Could not open $ARGV[0] to read, $!\n";
open( SEQ, $ARGV[1] ) or die "Could not open $ARGV[1] to read, $!\n";

my %hits = ();

while (<TAB>) {

    my @tmp = split( /\s+/, $_ );

    if ( $tmp[0] =~ /tem/ && $tmp[2] >= 50 ) {

        if ( $tmp[11] < $tmp[12] ) {

            if ( exists $hits{ $tmp[1] } ) {

                if ( $hits{ $tmp[1] }->{"start"} > $tmp[11] ) {

                    #print "$tmp[1], $tmp[11], $tmp[12]\n";
                    $hits{ $tmp[1] }->{"start"} = $tmp[11];
                }
                if ( $hits{ $tmp[1] }->{"end"} < $tmp[12] ) {
                    $hits{ $tmp[1] }->{"end"} = $tmp[12];
                }
            }
            else {
                $hits{ $tmp[1] }->{"start"} = $tmp[11];
                $hits{ $tmp[1] }->{"end"}   = $tmp[12];
            }
        }
        else {

            if ( exists $hits{ $tmp[1] } ) {

                if ( $hits{ $tmp[1] }->{"start"} > $tmp[12] ) {
                    $hits{ $tmp[1] }->{"start"} = $tmp[12];
                }
                if ( $hits{ $tmp[1] }->{"end"} < $tmp[11] ) {
                    $hits{ $tmp[1] }->{"end"} = $tmp[11];
                }
            }
            else {
                $hits{ $tmp[1] }->{"start"} = $tmp[12];
                $hits{ $tmp[1] }->{"end"}   = $tmp[11];
            }
        }
    }
}

$/ = "\n>";

while (<SEQ>) {

    chomp;

    if (/>?(\S+?)\n(.*)/s) {
        my $seqid = $1;
        my $seq   = $2;

        $seq =~ s/\n//g;

        #print ">$seqid\n$seq\n";
        #print STDERR "$seqid\n";
        if ( exists $hits{$seqid} ) {

            #print ">$seqid\n$seq\n";
            print STDERR "seqid=$seqid, start=", $hits{$seqid}->{"start"},
              " end=", $hits{$seqid}->{"end"}, ", length=",
              $hits{$seqid}->{"end"} - $hits{$seqid}->{"start"} + 1, "\n";

            my $subseq = substr(
                $seq,
                $hits{$seqid}->{"start"} - 1,
                $hits{$seqid}->{"end"} - $hits{$seqid}->{"start"} + 1
            );
            print ">$seqid\n$subseq\n";
        }

    }
}

close(SEQ);
close(TAB);
