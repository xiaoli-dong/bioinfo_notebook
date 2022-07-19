#!/usr/local/bin/perl

use strict;
use warnings;
use JSON;
use File::Slurp;
use Data::Dumper;
my $text = read_file( $ARGV[0] );

my $decoded = decode_json($text);

#print(from_json($decoded));

#print  Dumper($decoded);
my $perl_scalar = from_json( $text, { utf8 => 1 } );

#print  Dumper($perl_scalar);

foreach my $object ( keys %$decoded ) {

    #print $object, "\n";
    if ( $object eq "FragmentSize" ) {

        print join( "\n", @{ $decoded->{$object}->{"sizes"} } ), "\n";
    }

}

