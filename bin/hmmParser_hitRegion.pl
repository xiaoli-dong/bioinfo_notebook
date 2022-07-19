#!/usr/local/bin/perl

$#ARGV == 1

  or die "Usage: $0 <hmm vs nt output> <fasta sequence file>\n";

open( HMM, $ARGV[0] ) or die "Could not open $ARGV[0] to read, $!\n";
open( SEQ, $ARGV[1] ) or die "Could not open $ARGV[1] to read, $!\n";

$/ = "\nRANK ";
my %hits = ();

my $count = 0;
while (<HMM>) {
    $count++;

    if (/E_Value =\s+(\S+).*?T = (\S+).*?TS =\s+?(\d+?)\s+TE =\s+(\d+)/os) {
        my $evalue = $1;
        my $name   = $2;
        my $tstart = $3;
        my $tend   = $4;
        $scientific_notation = $evalue;
        $decimal_notation    = sprintf( "%.10g", $scientific_notation );
        if ( $decimal_notation < sprintf( "%.10g", "1e-009" ) ) {

            #print "$name\t$tstart\t$tend\t$evalue\n";

            if ( $tstart > $tend ) {
                $hits{$name}->{"start"} = $tend;
                $hits{$name}->{"end"}   = $tstart;
            }
            else {
                $hits{$name}->{"start"} = $tstart;
                $hits{$name}->{"end"}   = $tend;
            }
        }
    }

}

close(HMM);

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
