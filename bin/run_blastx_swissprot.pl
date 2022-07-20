#!/usr/bin/perl

$input_dir = $ARGV[0];
$output_dir =  $ARGV[1];
$db = $ARGV[2];
$num_procs = $ARGV[3];

@cmds = ();

opendir(INPUTDIR, "$input_dir") || die "Cannot opendir $input_dir: $!";

@dir_data=grep {!/^\.+$/} readdir(INPUTDIR);


foreach $name (@dir_data) {
    
    print STDERR "file $name found\n";
    
    my $query,$megablast_output;
    
    $query = "$input_dir/$name";
    $blast_output = "$output_dir/$name.out";
    
    $blast_cmd = "blastall -p blastx -e 1e-3 -i $query  -d $db -v 100 -F F -b 50 -a $num_procs  > $blast_output";

    push(@cmds, $blast_cmd);

    #print STDERR "$query\n $megablast_output\n";
}

close(INPUTDIR);


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
