#!/usr/bin/perl

@ARGV == 2
  or die
"Usage: $0 <read2taxonpath output from megan> <rank name: domain, phylum, class, order, family>\n",
  "\n";

open( MEGAN, $ARGV[0] ) or die "Could not open $ARGV[0] to read\n";

my %contig2taxons = ();
my %contigs       = ();
while (<MEGAN>) {
    chomp;
    my ( $name, $taxon ) = $_ =~ /^(\w+?),\s+?(\w+.*)/;
    my @array = split( /\_/, $name );

    #print STDERR "$array[1]\_$array[2]", "\n";
    #$arry[0] = contigname, array[1,2] position

    $contig2taxons{ $array[0] }->{"$array[1]\_$array[2]"} = $taxon;
    $contigs{ $array[0] }->{$taxon}++;
}
close(MEGAN);

#foreach my $contig (sort keys %contig2taxons){
#   foreach my $pos (keys %{$contig2taxons{$contig}}){
#	print "$contig\t$pos\t$contig2taxons{$contig}->{$pos}\n";
#   }
#}
my %taxon2ids = ();

my $var = $ARGV[1];
open( TAXON2IDS, ">taxon2ids.$var.txt" );
open( ID2TAXON,  ">id2taxon.$var.txt" );
open( ID2NAME,   ">id2name.$var.txt" );

print ID2TAXON "contigname\tgenecount\ttaxonName\ttaxonHits\n";
print ID2NAME "id\tname\n";

foreach my $contig ( sort keys %contigs ) {

    foreach
      my $tax ( sort { $contigs{$contig}->{$b} <=> $contigs{$contig}->{$a} }
        keys %{ $contigs{$contig} } )
    {
        my $size  = keys %{ $contig2taxons{$contig} };
        my $count = $contigs{$contig}->{$tax};
        if ( $var =~ /domain/ && $count / $size >= 0.5 ) {
            print ID2TAXON "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
            $taxon2ids{$tax} .= "$contig;";
        }
        elsif ( $var =~ /phylum/ && $count / $size >= 0.5 ) {
            if ( split( /;/, $tax ) > 3 ) {
                print ID2TAXON
                  "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
                $taxon2ids{$tax} .= "$contig;";
            }
        }
        elsif ( $var =~ /class/ && $count / $size >= 0.5 ) {
            if ( split( /;/, $tax ) > 4 ) {
                print ID2TAXON
                  "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
                $taxon2ids{$tax} .= "$contig;";
            }
        }
        elsif ( $var =~ /order/ && $count / $size >= 0.4 ) {
            if ( split( /;/, $tax ) > 5 ) {
                print ID2TAXON
                  "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
                $taxon2ids{$tax} .= "$contig;";
            }

        }
        elsif ( $var =~ /family/ && $count / $size >= 0.34 ) {
            if ( $tax =~ /aceae;/ && split( /;/, $tax ) > 6 ) {
                print ID2TAXON
                  "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
                $taxon2ids{$tax} .= "$contig;";
            }
        }
        elsif ( $var =~ /genus/ && $count / $size >= 0.34 ) {
            if ( split( /;/, $tax ) > 7 ) {
                print ID2TAXON
                  "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
                $taxon2ids{$tax} .= "$contig;";
            }
        }
        elsif ( $var =~ /species/ && $count / $size >= 0.34 ) {
            if ( split( /;/, $tax ) > 8 ) {
                print ID2TAXON
                  "$contig\t$size\t$tax\t$contigs{$contig}->{$tax}\n";
                $taxon2ids{$tax} .= "$contig;";
            }
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

print TAXON2IDS "taxon\tcount\tpercentage%\tids\t\n";

foreach my $taxon ( sort { $taxon2idcount{$b} <=> $taxon2idcount{$a} }
    keys %taxon2idcount )
{
    my $size    = $taxon2idcount{$taxon};
    my $percent = sprintf( "%0.2f", 100 * ( $size / $Total_count ) );
    print TAXON2IDS "$taxon\t$size\t$percent\%\t", $taxon2ids{$taxon}, "\n";

    if ( $percent > 1 ) {
        my @ids   = split( /;/, $taxon2ids{$taxon} );
        my @terms = split( /;/, $taxon );
        foreach my $id (@ids) {
            print ID2NAME "$id\t$terms[-1]\n";
        }
    }

}

