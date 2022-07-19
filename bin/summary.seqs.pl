#!/usr/local/bin/perl

$#ARGV == -1
  and die "Usage: $0 <fasta sequence file> <start length> <bin space>\n";

open( SEQFILE, "$ARGV[0]" ) or die "$!";

my $seq_count  = 0;
my $base_count = 0;
my %seq2len    = ();

$/ = "\n>";

while (<SEQFILE>) {

    chomp;
    if ( !/>?(\w+).*?\n(.+)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    my $seqid = $1;
    my $seq   = $2;

    $seq =~ tr/ \r\n\t//d;
    $seq_length = length($seq);

    $seq_count++;
    $base_count += $seq_length;

    $seq2len{$seqid} = $seq_length;

}
close(SEQFILE);

print "input sequence file name: $ARGV[0]\n";
print "Total sequence count: $seq_count\n";
print "Total bases count: $base_count\n";

my $len       = $ARGV[1];
my %len2count = ();

foreach ( sort { $seq2len{$a} <=> $seq2len{$b} } keys %seq2len ) {

    if ( $seq2len{$_} <= $len ) {
        $len2count{$len}++;
    }
    else {

        while ( $len < $seq2len{$_} ) {
            $len += $ARGV[2];
        }

        $len2count{$len}++;
    }

}

#produce the data for chart
print "var data = [\n";

my $i        = 1;
my $keyCount = keys %len2count;

foreach ( sort { $a <=> $b } keys %len2count ) {

    if ( $i < $keyCount ) {
        print "[$_,$len2count{$_}],\n";
    }
    else {
        print "[$_,$len2count{$_}]\n";
    }
    $i++;

}
print "];\n";

$i = 1;

foreach ( sort { $a <=> $b } keys %len2count ) {

    print "\'$_\',";

}
print "\n";

foreach ( sort { $a <=> $b } keys %len2count ) {

    print "$len2count{$_},";

}
print "\n";
