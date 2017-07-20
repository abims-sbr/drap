# $Id$

=pod
=head1 NAME

LocalScheduler - Manage local job submission.

=head1 DESCRIPTION

 This class allows to submit commands with the computer standard jobs scheduler.

=head1 VERSION

 1.0.0

=head1 AUTHOR

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package LocalScheduler ;
use strict ;
use warnings ;
use AbstractScheduler ;
our @ISA = ("AbstractScheduler");


=head2 function new

 Usage        : $cmd = new LocalScheduler( [$working_dir] )
 Function     : Creates, and returns a new scheduler.
 Returns      : [LocalScheduler] The scheduler.
 Args         : [str] The path to the working directory.
                [cpu] The number of cpu for local parallel submission.
                      If not set run as jobs as number of cpu cores.

=cut
sub new {
	my ($class, $working_dir, $cpu) = @_ ;

	if( !defined($cpu) ){ 
		$cpu = 0;
	}
	my $self = $class->SUPER::new( $working_dir );
 	$self->{cpu} = int($cpu) ;
	bless( $self, $class );
	return $self ;
}


=head2 procedure parallel_submit

 Usage        : $scheduler->parallel_submit( $cmdSet_name, $cmdSet_max_cpu, $array, @commands )
 Function     : The implementation of this abstract method is used to submit
                all the commands lines in parallel execution mode.
 Args         : [str] The name for the set of command.
                [int] The number of CPU of the greediest command of the set.
                [str] The command has or is an array command.
                [array] The list of executed commands lines.

=cut
sub parallel_submit {
	my ($self, $cmdSet_name, $cmdSet_max_cpu, $array, @commands) = @_ ;
	my $uniq_id = time."_".int(rand(10000)) ;
	my $basename = 'tmpCmd';
	if( defined($cmdSet_name) ){
		$basename = $cmdSet_name;
	}
	my $cmd_file_global = $self->{'working_dir'}."/".$basename."_".$uniq_id.".sh" ;
	my $cmd_file_wrapper = $self->{'working_dir'}."/".$basename."_".$uniq_id."_wrapper.sh" ;
	my $cmd_file_log = $self->{'working_dir'}."/".$basename."_".$uniq_id.".parallel.log" ;

	# Create command file content
	my $cmd_file_global_content = "" ;
	foreach my $current_command (@commands) {
		$cmd_file_global_content .= $current_command->{'cmd'}."\n" ;
	}

	# Write command file
	my $available_cpu = 1;
	if ( $self->{'cpu'} > 0 ){
		$available_cpu = $self->{'cpu'};
	} elsif ( $self->{'cpu'} == -1 ) {
		chomp($available_cpu = `grep -c '^processor[[:space:]]*:' /proc/cpuinfo`);
	}
	my $parallel_task;
	if ( $array eq 'has' ){
		$parallel_task = 1 ;
	} elsif ($array eq 'is') {
		$parallel_task = $available_cpu;
	} else {
		if ( $available_cpu < $cmdSet_max_cpu ){
			$parallel_task = int($available_cpu / ($cmdSet_max_cpu / scalar(@commands)));
			if ( $parallel_task < 1 ) { $parallel_task = 1; }
		}
		else {
			$parallel_task = $available_cpu;
		}
	}
	open( my $FH_cmd, ">", $cmd_file_global ) or die "Cannot create ".$cmd_file_global ;
	print $FH_cmd $cmd_file_global_content ;
	close( $FH_cmd );
	open( my $FH_cmd_wrapper, ">", $cmd_file_wrapper ) or die "Cannot create ".$cmd_file_wrapper ;
	print $FH_cmd_wrapper "#!/bin/bash\n\n" ;
	print $FH_cmd_wrapper "cat $cmd_file_global | parallel --no-notice --joblog $cmd_file_log --jobs $parallel_task\n" ;
	print $FH_cmd_wrapper "if [ \"`sed 1d $cmd_file_log | awk '\$7>0'`\" ]; then\n\techo 'Following command(s) failed:'\n\tsed 1d $cmd_file_log | awk '\$7>0' | cut -f9\n\texit 1;\nfi\n" ;
	close( $FH_cmd_wrapper );

	# Submit
	`bash $cmd_file_wrapper > $cmd_file_global.o 2> $cmd_file_global.e` ;

	# Check status
	$? and die ;

	# Delete tmp files
	unlink($cmd_file_wrapper);
	unlink($cmd_file_global, "$cmd_file_global.o", "$cmd_file_global.e");
	unlink($cmd_file_log);
}


=head2 procedure serial_submit

 Usage        : $scheduler->serial_submit( $cmdSet_name, @commands )
 Function     : The implementation of this abstract method is used to submit
                all the commands lines in serial execution mode.
 Args         : [str] The name for the set of command.
                [array] The list of executed commands lines.

=cut
sub serial_submit {
	my ($self, $cmdSet_name, @commands) = @_ ;
	my $cmd_file = $self->{'working_dir'}."/tmpCmd_".time."_".int(rand(10000)).".sh" ;
	if( defined($cmdSet_name) ){
		$cmd_file = $self->{'working_dir'}."/".$cmdSet_name."_".time."_".int(rand(10000)).".sh" ;
	}
	
	# Create command file content
	my $cmd_file_content = "#!/bin/bash\n\n" ;
	my $command_idx = 1 ;
	foreach my $current_command (@commands) {
		$cmd_file_content .= $current_command->{'cmd'}." ;\nif [ \$? -ne 0 ]; then\n\techo \"The command ".$command_idx." in ".$cmd_file." has failed\" >&2\n\texit 1;\nfi\n\n" ;
		$command_idx++ ;
	}
		
	# Write command file
	open( my $FH_cmd, ">", $cmd_file ) or die "Cannot create ".$cmd_file ;
	print $FH_cmd $cmd_file_content ;
	close( $FH_cmd );
	
	# Submit
	`bash $cmd_file > $cmd_file.o 2> $cmd_file.e` ;

	# Check status
	$? and die ;

	# Delete tmp files
	unlink($cmd_file, "$cmd_file.o", "$cmd_file.e");
}


1;