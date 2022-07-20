#!/usr/bin/perl

$blast_result        = $ARGV[0];
$identity_threshold  = $ARGV[1];
$align_len_threshold = $ARGV[2];

%seq;

opendir( DIRHANDLE, "$blast_result" ) || die "Cannot opendir $blast_result: $!";

@dir_data = grep { !/^\.+$/ } readdir(DIRHANDLE);

foreach $name (@dir_data) {

    print STDERR "found file: $name\n";

    open( BLASTRESULT, "<$blast_result/$name" )
      || die "cannot open the file $name\n";

    while ( $entry = <BLASTRESULT> ) {

        if ( ( $query_id, $subject_id, $pct, $m_len ) =
            $entry =~ /^(\S+)\s+(\S+)\s+(\d+\.\d+)\s+(\d+).*?$/s )
        {

            #print "$query_id, $subject_id, $pct, $m_len\n";

            if (   $pct >= $identity_threshold
                && $m_len >= $align_len_threshold
                && not defined $seq{$query_id} )
            {

                $seq{$query_id} = 1;
            }

        }
    }

    close(BLASTRESULT);
}

close(DIRHANDLE);

$num_seq = scalar( keys %seq );

print STDERR
"$num_seq sequences has hit with ribosomal protein gene with identity greater than $identity_threshold and alignment length is longer than $align_threshold\n";

