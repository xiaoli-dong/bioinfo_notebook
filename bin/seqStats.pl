#!/usr/local/bin/perl

use strict;	
use warnings;
use List::Util qw(sum);
#use Statistics::Descriptive;
use Getopt::Long;

my ($seq, $format);

$format = "fasta";

&GetOptions(
    "f=s" =>\$format,
    "s=s"       => \$seq
    );

if(not defined $seq ){
    die "Usage: $0  <-f fasta|fastq> <-s fasta format sequence file> \n";
}
print STDERR "$0 -f $format -s $seq\n";
open(SEQ, "$seq") or die "Could not open $seq file to read, $!\n";
my @lengths = ();
my @gcs = ();

if($format eq "fasta"){
    
    $/ = "\n>";
    while (<SEQ>){
	chomp;
	if(! />?(\S+).*?\n(.+)/s){
	    die "Could not read fasta record #$.: $_\n";
	}
	my $id = $1;
	my $sequence = $2;
	$sequence =~ s/\s+//g;
	my $len = length ($sequence);
	push(@lengths, $len);
	push(@gcs, calcgc($sequence));
    }
    $/ = "\n";
}
elsif($format eq "fastq"){
    while (<SEQ>){
	chomp;
	
	my $id = $1;
	my $sequence = <SEQ>;
	<SEQ>;
	<SEQ>;
	$sequence =~ s/\s+//g;
	my $len = length ($sequence);
	push(@lengths, $len);
	push(@gcs, calcgc($sequence));
    }
}
close(SEQ);

my @sorted_lengths = sort { $a <=> $b }@lengths;
my $total_bases = sum(@lengths);
my $total_seqCount = @lengths;
my $min_len = $sorted_lengths[0];
my $max_len = $sorted_lengths[$#sorted_lengths];
my $mean_value = mean(\@lengths);
my $median_value = median(\@lengths);
my $stdev_value = stdev(\@lengths);
my $n50 = get_N50(\@lengths);

print "Total bases\tTotal sequence count\tMin length\tMax length\tMean length\tMedian length\tstdev\tN50\n";
print "$total_bases\t$total_seqCount\t$min_len\t$max_len\t$mean_value\t$median_value\t$stdev_value\t$n50\n";

sub mean {
    @_ == 1 or die ('Sub usage: $average = average(\@array);');
    my ($array_ref) = @_;
    my $sum;
    my $count = scalar @$array_ref;
    foreach (@$array_ref) { $sum += $_; }
    return sprintf("%.2f",$sum / $count);
}

sub median {
    @_ == 1 or die ('Sub usage: $median = median(\@array);');
    my ($array_ref) = @_;
    my $count = scalar @$array_ref;
# Sort a COPY of the array, leaving the original untouched
    my @array = sort { $a <=> $b } @$array_ref;
    if ($count % 2) {
	return $array[int($count/2)];
    } else {
	return sprintf("%.2f", ($array[$count/2] + $array[$count/2 - 1]) / 2);
    }
} 

sub stdev{
    my ($array_ref) = @_;
    my $sqsum = 0;
    my $n = @$array_ref;
    my $mean_value = sum(@$array_ref)/$n;
    
    for (@$array_ref) {
	$sqsum += ( $_ ** 2 );
    } 
    $sqsum /= $n;
    $sqsum -= ( $mean_value ** 2 );
    my $stdev = sqrt($sqsum);
    return sprintf("%.2f", $stdev);
}


sub get_N50{
    my ($arrRef) = @_;
    my @sort = sort {$b <=> $a} @$arrRef;
    my $totalLength = sum(@sort);
    my $n50 = 0;
    my $n50_value = 0;
    foreach my $val(@sort){
	$n50+=$val;
	if($n50 >= $totalLength/2){
	    #print "N50 length is $n50 and N50 value is: $val\n";
	    $n50_value = $val;
	    last;
	}
    }
    return $n50_value;
}

###################################################################
# _histogram_bins - calculates the bins usings Scott's algorithm
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Number of bins to create. If $nbins is undef
#                    the number of bins is calculated using Scott's
#                    algorithm
#
###################################################################
sub _histogram_bins {
	my ( $data, $nbins ) = @_;

	if( !defined $data ) { return; }

	my $calcBins = ( defined $nbins )? 0 : 1;
	my $cnt = 0;
	my $mean= 0;
	my $max = my $min = $data->[0];
	foreach (@$data) {
		$mean += $_;
		$min = ( $_ < $min )? $_ : $min;
		$max = ( $_ > $max )? $_ : $max;
		$cnt++;
	}
	$mean /= $cnt if( $cnt > 1 );

	my $sumsq = 0;
	$nbins = 1 if( $calcBins );
	my $s = 0;
	if( $cnt > 1 ) {
		foreach (@$data) {
			$sumsq += ( $_ - $mean )**2;
		}
		$s = sqrt( $sumsq / ($cnt - 1));
		$nbins = 3.49 * $s / $cnt**0.33 if( $s > 0 && $calcBins );
	}

	my $binwidth = ( $max - $min ) / $nbins;

	my $lower = $min;
	my $upper = $lower;

	my $bins;
	my @cutPoints;
	my $cntr = 0;
	while ( $upper <= $max && $cntr < $nbins) {
		$upper = $lower + $binwidth;
		push( @cutPoints, [int($lower), int($upper)] );
		$lower = $upper;
		$cntr++;
	}

	return \@cutPoints;
}
###################################################################
# _histogram_frequency - bins the data
#
#     Lower Boundry <= data value < Upper Boundry
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Vector containing the cutpoints to bin the data
#
###################################################################
sub _histogram_frequency {
	my ( $data, $cutPoints ) = @_;

	if( !defined $data || !defined $cutPoints ) { return; }

	my @freqs;
	foreach (@$cutPoints) {
		push( @freqs, 0 );
	}

	foreach (@$data) 
	{
		for( my $i = 0; $i < scalar( @$cutPoints ); $i++ ) 
		{
		if( ($_ >= $cutPoints->[$i]->[0] && $_ < $cutPoints->[$i]->[1])
			||
			($i == (scalar (@$cutPoints) - 1) && $_ >= $cutPoints->[$i]->[1]) ) 
			{	

				$freqs[$i]++;
			}
		}
	}
	return \@freqs;
}



sub calcgc {
    
    my ($seq) = @_;
    my $count = 0;
    my $len   = length($seq);
    for (my $i = 0;$i<$len; $i++) {
	my $base = substr $seq, $i, 1;
	$count++ if $base =~ /[G|C]/i;
    }
    return sprintf("%.2f",($count / $len) * 100);
}
