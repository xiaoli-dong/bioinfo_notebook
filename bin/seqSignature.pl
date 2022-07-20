#!/usr/local/bin/perl

use strict;	
use warnings;
use List::Util qw(sum);
#use Statistics::Descriptive;
use Getopt::Long;

my ($seq);

&GetOptions("s=s"       => \$seq
    );

if(not defined $seq ){
    die "Usage: $0 <-s fasta format sequence file> \n";
}

my $tetras = getUniqueTetraNucleotides();
my @unique_tetras = keys %$tetras;

print "id\t", join("\t", sort @unique_tetras), "\n"; 
open(SEQ, "$seq") or die "Could not open $seq file to read, $!\n";
my %seqInfo = ();
my %lengths = ();
my @gcs = ();

$/ = "\n>";

while (<SEQ>){
    chomp;
    if(! />?(\S+)(.*?)\n(.+)/s){
	die "Could not read fasta record #$.: $_\n";
    }
    my $id = $1;
    my $desc = $2;
    my $sequence = $3;
    $sequence =~ s/\s+//g;
    my $len = length ($sequence);
    $seqInfo{$id}->{gc} = calcgc($sequence);
    $seqInfo{$id}->{length} = length($sequence);
    $seqInfo{$id}->{desc} = $desc;
    my $freq = compute_tetranucleotides_freq($sequence,\@unique_tetras);
    print "$id";
    foreach my $u (sort keys %$freq){
	print "\t", sprintf("%.4f", $freq->{$u});
    }
    print "\n";
	    
}
close(SEQ);

#######################################################################################

sub RevComp{
    my ($seq) = @_;
    my $reverse_seq =reverse($seq);
    $reverse_seq =~ tr/ATGC/TACG/;
    return $reverse_seq;
}


#######################################################################################
sub compute_zscores_for_tetranucleotides{
  
}


#######################################################################################
sub compute_tetranucleotides_freq{
  
    my($theseq, $tetras) = @_;
    
    #$theseq = uc($theseq);
    #as described by Teeling et al. BMC Bioinformatics 2004, 5:163.
    #$theseq .=RevComp($theseq);
    
    my $seqLen = length($theseq);
    
    my %oligos4 = ();
    my %freqs = ();
    my $i=0;
    
    my $len=length($theseq)-4+1;
    
    while ($i<$len){
	$seq=substr($theseq,$i,4);
	#print STDERR "seq=$seq\t", length($seq), "\n";
	$oligos4{$seq}++;
	$i++;
    }
    
    my $revseq = RevComp($theseq);
    
    $i=0;
    while ($i<$len){
	$seq=substr($revseq,$i,4);
	$oligos4{$seq}++;
	$i++;
    }
    
    my $total = 0;
    
    foreach (keys %oligos4){
	
	$total += $oligos4{$_};
    }
    
    foreach my $item(@$tetras){
	
	if((not exists $oligos4{$item}) && (not exists $oligos4{RevComp($item)}) ){
	    
	    $freqs{$item} = 0;
	}
	elsif(not exists $oligos4{$item}){
	     #print STDERR length($item), "\t$item\t$oligos4{$item}\n";
	    $freqs{$item} = 136*$oligos4{RevComp($item)}/$total;
	}
	elsif(not exists $oligos4{RevComp($item)}){
	    
	    $freqs{$item} = 136 * $oligos4{$item}/$total;
	}
	else{
	   
	    $freqs{$item} = 136 * ($oligos4{$item} + $oligos4{RevComp($item)})/$total;
	}
    }
    return \%freqs;
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

sub getUniqueTetraNucleotides{
    
    my %tetras = ();
    my @n1a = ('A', 'C', 'G', 'T');
    my @n2a = ('A', 'C', 'G', 'T');
    my @n3a = ('A', 'C', 'G', 'T');
    my @n4a = ('A', 'C', 'G', 'T');
    my $i = 0;
    foreach my $n1 (@n1a){
	foreach my $n2 (@n2a){
	    foreach my $n3 (@n3a){
		foreach my $n4 (@n4a){
		    my $tetra =  "$n1$n2$n3$n4";
		    if(not exists $tetras{RevComp($tetra)}){
			$tetras{$tetra} = 1;
		    }
		}
	    }
	}
    }
    
    return \%tetras;
}
