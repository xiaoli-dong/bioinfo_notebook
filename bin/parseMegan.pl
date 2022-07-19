#!/usr/bin/perl

@ARGV == 4
  or die
"Usage: $0 <read2taxonpath output from megan> <fasta file used of megan annotation> <fasta contigs> <rank name: domain, phylum, class, order, family>\n",
  "\n";

open( MEGAN,  $ARGV[0] ) or die "Could not open $ARGV[0] to read\n";
open( FASTA,  $ARGV[1] ) or die "Could not open $ARGV[1] to read\n";
open( CONTIG, $ARGV[2] ) or die "Could not open $ARGV[2] to read\n";
my %contig2len = ();
$/ = "\n>";
while (<CONTIG>) {
    chomp;
    if ( !/>?(\S+).*?\n(.*)/s ) {
        die "Could not read FastA record #$.: $_\n";
    }
    my $name = $1;
    my $seq  = $2;
    $seq =~ s/\s+//g;
    $contig2len{$name} = length($seq);
}
close(CONTIG);
$/ = "\n";

my %contig2taxons = ();
my %contigs       = ();
my $var           = $ARGV[3];
while (<MEGAN>) {
    chomp;
    next if /@Parameters=/;

    my ( $name, $taxon ) = $_ =~ /^(\w+?)\t(\w+.*)/;
    my @array = split( /\_/, $name );

    #print STDERR "$array[1]\_$array[2]", "\n";
    #$arry[0] = contigname, array[1,2] position
    #superkingdom::Bacteria;
    my ($mytaxon) = $taxon =~ /$var\::(\w.*?);/;
    if ( $taxon =~ /$var/ ) {
        $contig2taxons{ $array[0] }->{"$array[1]\_$array[2]"} = $mytaxon;
        $contigs{ $array[0] }->{$mytaxon}++;
    }
}
close(MEGAN);
my $total_contig_count = 0;
my %contig2genepos     = ();
while (<FASTA>) {

    if (/^>(\w+?)\_(\d+?\_\d+)/) {
        $contig2genepos{$1}->{$2} = 1;
    }
}
close(FASTA);
my $total_contig_count = keys %contig2genepos;

#foreach my $contig (sort keys %contig2taxons){
#   foreach my $pos (keys %{$contig2taxons{$contig}}){
#	print "$contig\t$pos\t$contig2taxons{$contig}->{$pos}\n";
#   }
#}
my %taxon2ids        = ();
my %taxon2totalbases = ();

open( TAXON2IDS, ">taxon2ids.$var.txt" );
open( ID2TAXON,  ">id2taxon.$var.txt" );
open( ID2NAME,   ">id2name.$var.txt" );

print ID2TAXON "contigname\tgenecount\ttaxonName\ttaxonHits\n";
print ID2NAME "id\t$var\n";

foreach my $contig ( sort keys %contigs ) {

    foreach
      my $tax ( sort { $contigs{$contig}->{$b} <=> $contigs{$contig}->{$a} }
        keys %{ $contigs{$contig} } )
    {

        #how many annotated
        #	my $size = keys  %{$contig2taxons{$contig}} ;
        #	how many gene in contigs
        my $size  = keys %{ $contig2genepos{$contig} };
        my $count = $contigs{$contig}->{$tax};
        if ( $var =~ /superkingdom/ && $count / $size >= 0.5 ) {
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};
        }
        elsif ( $var =~ /phylum/ && $count / $size >= 0.5 ) {

            #if(split(/;/, $tax) >3){
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};

            #}
        }
        elsif ( $var =~ /class/ && $count / $size >= 0.5 ) {

            #if(split(/;/, $tax) >4){
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};

            #}
        }
        elsif ( $var =~ /order/ && $count / $size >= 0.4 ) {

            #if(split(/;/, $tax) >5){
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};

            #}

        }
        elsif ( $var =~ /family/ && $count / $size >= 0.34 ) {

            #if($tax =~ /aceae;/ && split(/;/, $tax) >6 ){
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};
        }

        elsif ( $var =~ /genus/ && $count / $size >= 0.34 ) {

            #if(split(/;/, $tax) >7){
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};

            #}
        }
        elsif ( $var =~ /species/ && $count / $size >= 0.34 ) {

            #if(split(/;/, $tax) >8){
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
            $taxon2totalbases{$tax} += $contig2len{$contig};

            #}
        }

    }
}

my %taxon2idcount = ();
my $total_count   = 0;

foreach ( sort keys %taxon2ids ) {
    my $size = split( /;/, $taxon2ids{$_} );
    $Total_count += $size;
    $taxon2idcount{$_} = $size;
}

print TAXON2IDS "taxon\tcontigcount\tpercentage%\ttotalbases\tids\t\n";

foreach my $taxon ( sort { $taxon2idcount{$b} <=> $taxon2idcount{$a} }
    keys %taxon2idcount )
{
    my $size    = $taxon2idcount{$taxon};
    my $percent = sprintf( "%0.2f", 100 * ( $size / $Total_count ) );
    print TAXON2IDS "$taxon\t$size\t$percent\%\t", $taxon2totalbases{$taxon},
      "\t", $taxon2ids{$taxon}, "\n";

    if ( $percent > 1 ) {
        my @ids   = split( /;/, $taxon2ids{$taxon} );
        my @terms = split( /;/, $taxon );
        foreach my $id (@ids) {
            print ID2NAME "$id\t$terms[-1]\n";
        }
    }

}

