#!/usr/local/bin/perl

use strict;	
use warnings;
use List::Util qw(sum);
#use Statistics::Descriptive;
use Getopt::Long;

my ($seq, $type, $header);

&GetOptions("s=s"       => \$seq,
	    "t=s" =>\$type,
	    "h=s" =>\$header
	    
    );

if(not defined $seq || not defined $type || not defined $header){
    die "Usage: $0 <-s fasta format sequence file> <-t 454|sanger|contig> <-h header informaiton for table and page\n";
}
open(SEQ, "$seq") or die "Could not open $seq file to read, $!\n";
my @lengths = ();
my @gcs = ();

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
close(SEQ);

my @bins = ();

if($type eq "contig"){
    push (@bins, [200,900]);
    push(@bins, [900,2000]);
    push(@bins, [2000,3000]);
    push(@bins, [3000,4000]);
    push(@bins, [4000,5000]);
    push(@bins, [5000,6000]);
    push(@bins, [6000,7000]);
    push(@bins, [7000,8000]);
    push(@bins, [8000,9000]);
    push(@bins, [9000, 10000]);
    push(@bins, [10000, 20000]);
    push(@bins, [20000, 30000]);
    push(@bins, [30000, 40000]);
    push(@bins, [40000, 50000]);
    push(@bins, [50000, 60000]);
    push(@bins, [60000, 70000]);
    push(@bins, [70000, 80000]);
    push(@bins, [80000, 90000]);
    push(@bins, [90000, 100000]);
    push(@bins, [100000, 1000000]);
}
elsif($type eq "454"){
     push (@bins, [0,50]);
     push (@bins, [50,100]);
     push (@bins, [100,150]);
     push (@bins, [150,200]);
     push (@bins, [200,250]);
     push (@bins, [250,300]);
     push (@bins, [300,350]);
     push (@bins, [350,400]);
     push (@bins, [400,450]);
     push (@bins, [450,500]);
     push (@bins, [500,550]);
     push (@bins, [550,600]);
     push (@bins, [600,1000]);
}
elsif($type eq "sanger"){
     push (@bins, [0,50]);
     push (@bins, [50,100]);
     push (@bins, [100,150]);
     push (@bins, [150,200]);
     push (@bins, [200,250]);
     push (@bins, [250,300]);
     push (@bins, [300,350]);
     push (@bins, [350,400]);
     push (@bins, [400,450]);
     push (@bins, [450,500]);
     push (@bins, [500,550]);
     push (@bins, [550,600]);
     push (@bins, [600, 650]);
     push (@bins, [650,700]);
     push (@bins, [700, 750]);
     push (@bins, [750,800]);
     push (@bins, [800, 850]);
     push (@bins, [850,900]);
     push (@bins, [900, 950]);
     push (@bins, [950,1000]);
     push (@bins, [1000,1050]);
     push (@bins, [1050,1100]);
     push (@bins, [1100,1150]);
     push (@bins, [1150,1200]);
     push (@bins, [1200,1250]);
     push (@bins, [1250,1300]);
     push (@bins, [1300,1350]);
     push (@bins, [1350,1400]);
     push (@bins, [1400,1450]);
     push (@bins, [1450,1500]);
     push (@bins, [1500,1550]);
     push (@bins, [1550,2000]);
     push (@bins, [2000,2500]);
     push (@bins, [2500,3000]);
     push (@bins, [3000,10000]);
}
my $histogram_bin = 20;
my $cp = _histogram_bins( \@lengths, $histogram_bin);
my $binArrRef = _histogram_frequency( \@lengths, \@bins);
my $count = 0;
my $len_cat = "";
my $len_data ="";
my $len_table = "length_range(bp)\tcount\n";

for my $bin (@bins)
{
    # push(@labelArr, _numformat( $bin->[0] + ($bin->[1] - $bin->[0])/2 ) );
    
    
    #print "$bin->[0]-$bin->[1]\t $binArrRef->[$count]\n" unless $binArrRef->[$count] == 0;
    $len_cat .= "\'$bin->[0]-$bin->[1]\', " unless $binArrRef->[$count] == 0;
    $len_data .= "$binArrRef->[$count]," unless $binArrRef->[$count] == 0;
    $len_table .= "$bin->[0]-$bin->[1]\t $binArrRef->[$count]\n" unless $binArrRef->[$count] == 0;
    #print "$bin->[0]-$bin->[1]\t $binArrRef->[$count]\n" unless $binArrRef->[$count] == 0;
    $count++;
}
$len_cat =~ s/,$//;
$len_data =~ s/,$//;


#print join("\n", @gcs);
my @sorted_gcs = sort { $a <=> $b } @gcs;
my $gc_min = $sorted_gcs[0];
my $gc_max = $sorted_gcs[$#sorted_gcs];
$histogram_bin = $gc_max -$gc_min + 1;
$cp = _histogram_bins( \@gcs,$histogram_bin );
$binArrRef = _histogram_frequency( \@gcs, $cp );
$count = 0;
my $gc_cat = "";
my $gc_data = "";
my $gc_table = "g+c%\tcount\n";
for my $bin (@$cp)
    {
	# push(@labelArr, _numformat( $bin->[0] + ($bin->[1] - $bin->[0])/2 ) );
	#print "$bin->[1]\t$binArrRef->[$count]\n";
	$gc_cat .= "\'$bin->[1]\', ";
	$gc_data .= "$binArrRef->[$count],";
	$gc_table .= "$bin->[1]\t$binArrRef->[$count]\n";
	$count++;

    }
$gc_cat =~ s/,$//;
$gc_data =~ s/,$//;

my @sorted_lengths = sort { $a <=> $b }@lengths;
my $total_bases = sum(@lengths);
my $total_seqCount = @lengths;
my $min_len = $sorted_lengths[0];
my $max_len = $sorted_lengths[$#sorted_lengths];
my $mean_value = mean(\@lengths);
my $median_value = median(\@lengths);
my $stdev_value = stdev(\@lengths);
my $n50 = get_N50(\@lengths);

    
print  <<END;

<html> 
  <head> 
    <title>$header</title>
    
    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>
    <script type="text/javascript" src="/HMP/Highcharts/js/highcharts.js"></script>
    <script type="text/javascript" src="/HMP/Highcharts/js/modules/exporting.js"></script>
    <script type="text/javascript" src="/HMP/js/column.config.js"></script>
    
    <script>
    var charts = [];
\$(document).ready(function() {
    var len_cat = [$len_cat];
    var len_data = [{
      name: 'length',
      data: [$len_data]
	  }];
    
    var gc_cat = [$gc_cat];
    var gc_data = [{
      name: 'G+C%',
      data: [$gc_data]
	  
		   }];
  
    var title_stats = "Taxonomic Classification Stats";
    var title_kingdom = "Superkingdom Distribution";
    var title_bphylum = "Bacteria Phylum Distribution";
    var title_aphylum = "Archaea Phylum Distribution";
    var dis = 10;
 //now, creating a new chart is easy!
     charts.push(new Highcharts.Chart(
		     getChartConfig("length", "Length Distribution", len_data, dis, len_cat)
		 ));
    
 
 charts.push(new Highcharts.Chart(
 getChartConfig("gc", "G+C% Content Distribution", gc_data, dis, gc_cat)
 ));
 
 });

</script>


 </head>

    <body topmargin="50">
    
    <table border=0 width=80% align="center">
    <tr><th colspan=2 height=50>$header</th></tr>
    <tr><td>Bases</td><td>$total_bases</td><td></tr>
    <tr><td>Reads</td><td>$total_seqCount</td><td></tr>
    <tr><td>Max Length</td><td>$max_len</td><td></tr>
    <tr><td>Min Length</td><td>$min_len</td><td></tr>
    <tr><td>Median Length</td><td>$median_value</td><td></tr>
    <tr><td>Mean Length</td><td>$mean_value</td><td></tr>
    <tr><td>N50</td><td>$n50</td><td></tr>
    <tr><td>stdev</td><td>$stdev_value</td><td></tr>  
    </table>
    <div id="length"></div>
    <p>Length summary table</p>
    <pre>
    $len_table
    </pre>
    <div id="gc"></div>
    <p>G+C% summary table</p>
    <pre>
    $gc_table
    </pre>
    </body>
    </html>
END

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
