#!/usr/bin/perl

$infile = $ARGV[0];
$outdir =  $ARGV[1];
$subject = $ARGV[2];
$num_procs = $ARGV[3];
$program = $ARGV[4];

chomp ($infile);

open (INFILE, "<$infile");

$/ = "\n>";

$count = 0;

$file_count = 0;



$ssaha_input = "";

@cmds = ();
#split input file into multiple files

while($entry = <INFILE>){
    
    chomp($entry);
    
    if($count == 0){
	$file_count++;
	$ssaha_input = "$outdir/input/$file_count.input.fasta";
	open (OUTPUT, ">$ssaha_input") || die "cannot open the file $ssaha_input\n";
    }
   
    if(($name, $seq) = $entry =~ /^>?(\S+).*?\n(.*)/s){
	print OUTPUT ">$name\n$seq\n";
	$count++;
    }

    
    if($count == 14060){
	close(OUTPUT);
	
	my $query,$ssaha_output;
	$query = "$outdir/input/$file_count.input.fasta";
	$ssaha_output = "$outdir/output/$file_count.ssaha.out";
	push(@cmds, "/export/data/projects/brassica/$program $query $subject -parserFriendly -mp 100 > $ssaha_output");
	$count = 0;

    }
   
}


close(OUTPUT);

my $query,$ssaha_output;
$query = "$outdir/input/$file_count.input.fasta";
$ssaha_output = "$outdir/output/$file_count.ssaha.out";
push(@cmds, "/export/data/projects/brassica/$program $query $subject -parserFriendly -mp 100 > $ssaha_output");


close(INFILE);


&run_cmds($num_procs, @cmds);


sub run_cmds {

        my ($max_cmds, @cmds) = @_;

        my ($num_children, $pid);

        return unless @cmds; # get out of the function if there is nothing in @cmds

        for($num_children = 1; $num_children < $max_cmds && @cmds; $num_children++){
        # initialize the number of child processes at 1, and increment it by one
        #while it is less than $max_cmds

                my $cmd = shift (@cmds);
                if($pid = fork) {
                        # do nothing if parent
                } elsif(defined $pid) { # $pid is zero here if defined
                        system $cmd;
                        exit;
                } else {
                        #weird fork error
                        die "Can't fork: $!\n";
                }
        }

        while(@cmds) {
                undef $pid;
                FORK: {
                        my $cmd = shift (@cmds);
                        if($pid = fork) {
                                # parent here
                                $num_children++;
                                wait;
                                $num_children--;
                                next;

                        } elsif(defined $pid) { # $pid is zero here if defined
                                system $cmd;
                                exit;

                        } else {
                                #weird fork error
                                die "Can't fork: $!\n";
                        }
                }
        }
        wait while $num_children--;
}
