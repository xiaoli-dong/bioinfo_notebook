#!/usr/bin/perl

@ARGV == 3
  or die "Usage: $0 <fasta file> <number of files> <output directory>\n",
  "split a fasta file into multiple files\n";

$infile     = $ARGV[0];
$num_files  = $ARGV[1];
$output_dir = $ARGV[2];

chomp($infile);

open( INFILE, "<$infile" );

#count sequencs it has in file

while ( $entry = <INFILE> ) {

    chomp($entry);

    if ( ( $name, $seq ) = $entry =~ /^>.*?/ ) {
        $count++;
    }
}

print "$count sequences we have\n";

close(INFILE);

$remainder    = $count % $num_files;
$seq_per_file = ( $count - $remainder ) / $num_files;

open( INFILE, "<$infile" );

$/ = "\n>";

$count = 0;

$file_count = 1;

$output_file = "$output_dir/$file_count.fasta";
open( OUTPUT, ">$output_file" ) || die "cannot open the file $output_file\n";

#split input file into multiple files

while ( $entry = <INFILE> ) {

    chomp($entry);

    if ( ( $name, $seq ) = $entry =~ /^>?(\S+).*?\n(.*)/s ) {
        print OUTPUT ">$name\n$seq\n";
        $count++;
    }

    if ( $file_count != $num_files && $count == $seq_per_file ) {
        close(OUTPUT);
        $count = 0;

        if ( $file_count < $num_files ) {
            $file_count++;
            $output_file = "$output_dir/$file_count.fasta";
            open( OUTPUT, ">$output_file" )
              || die "cannot open the file $output_file\n";
        }
    }

}

close(INFILE);
close(OUTPUT);

