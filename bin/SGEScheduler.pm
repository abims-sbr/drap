# $Id$

=pod
=head1 NAME

SGEScheduler - Manage job submission on SGE.

=head1 DESCRIPTION

 This class allows to submit commands with SGE on HPC.

=head1 VERSION

 1.0.1

=head1 AUTHOR

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package SGEScheduler ;
use strict ;
use warnings ;
use List::Util ;
use AbstractScheduler ;
our @ISA = ("AbstractScheduler");


=head2 function new

 Usage        : $cmd = new SGEScheduler( [$working_dir] )
 Function     : Creates, and returns a new scheduler.
 Returns      : [SGEScheduler] The scheduler.
 Args         : [str] The path to the working directory.

=cut
sub new {
	my ($class, $working_dir) = @_ ;

 	my $self = $class->SUPER::new( $working_dir );
	bless( $self, $class );
	return $self ;
}


=head2 function _get_multithread_memory

 Usage        : SGEScheduler::_get_multithread_memory( $mem, $cpu )
 Function     : Return the memory string for h_vmem and mem submission 
                parameters. The global memory consumption is divided by the 
                number of CPUs.
 Returns      : [str] The value to use in h_vmem and mem submission parameters.
 Args         : [str] The global memory provided (example: "8g").
                [int] The number of CPUs used.

=cut
sub _get_multithread_memory {
	my ($mem, $nb_threads) = @_ ;
	if( $mem =~ /^(\d+)(\w+)$/ ){
		my $size = $1 ;
		my $unity = lc($2);
		if( ($size % $nb_threads) > 1 && $unity ne "b" ){
			my @unity_notation = ();
			if( length($unity) == 2 ){
				@unity_notation = ( "tb", "gb", "mb", "kb", "b" );
			} else {
				@unity_notation = ( "t", "g", "m", "k", "b" );
			}
			my $new_unity = undef ;
			for( my $idx = 0 ; $idx < scalar(@unity_notation) ; $idx++ ) {
				if( $unity_notation[$idx] eq $unity ) {
					$new_unity = $unity_notation[$idx+1] ;
				}
			}
			$mem = int(($size/$nb_threads)*1024).$new_unity ;
		} else {
			$mem = ($size/$nb_threads).$unity ;
		}
	} else {
		die "Error in memory mangement" ;
	}
	
	return $mem ;
}


=head2 function _get_ressources_opt

 Usage        : SGEScheduler::_get_ressources_opt( $mem, $virtual_mem, $cpu )
 Function     : Returns the submission parameters for to book the memory, the 
                virtual memory and the number of CPUs on the same computer used
                by the commands.
 Returns      : [str] The value to use in h_vmem and mem submission parameters.
 Args         : [str] The global memory provided (example: "8g").
                [int] The number of CPUs used.

=cut
sub _get_ressources_opt {
	my ( $mem, $virtual_mem, $cpu ) = @_ ;
	my $ressources_opt = "" ;
	
	if( $cpu ne "" ){
		$ressources_opt .= " -pe parallel_smp ".$cpu ;
		if( $virtual_mem ne "" ){
			$virtual_mem = SGEScheduler::_get_multithread_memory( $virtual_mem, int($cpu) );
		}
		if( $mem ne "" ){
			$mem = SGEScheduler::_get_multithread_memory( $mem, int($cpu) );
		}
	}
	$ressources_opt .= $virtual_mem ne "" ? " -l h_vmem=".$virtual_mem : "" ;
	$ressources_opt .= $mem ne "" ? " -l mem=".$mem : "" ;
	
	return $ressources_opt ;
}


=head2 procedure parallel_submit

 Usage        : $scheduler->parallel_submit( $cmd_set_name, @commands )
 Function     : The implementation of this abstract method is used to submit
                all the commands lines in parallel execution mode.
 Args         : [str] The name for the set of command.
                [array] The list of executed commands lines.

=cut
sub parallel_submit {
	my ($self, $cmdSet_name, @commands) = @_ ;
	my $uniq_id = time."_".int(rand(10000)) ;
	my $cmd_file_global = $self->{'working_dir'}."/tmpCmd_".$uniq_id.".sh" ;
	my $cmd_file_wrapper = $self->{'working_dir'}."/tmpCmd_".$uniq_id."_wrapper.sh" ;
	my $name_opt = "" ;
	if( defined($cmdSet_name) ){
		$name_opt = "-N ".$cmdSet_name ;
		$cmd_file_global = $self->{'working_dir'}."/".$cmdSet_name."_".$uniq_id.".sh" ;
		$cmd_file_wrapper = $self->{'working_dir'}."/".$cmdSet_name."_".$uniq_id."_wrapper.sh" ;
	}

	# Create command file content
	my $cmd_file_global_content = "" ;
	my $virtual_mem = undef ;
	my $mem = undef ;
	my $cpu = undef ;
	foreach my $current_command (@commands) {
		$cmd_file_global_content .= $current_command->{'cmd'}."\n" ;
		$virtual_mem = $current_command->{'virtual_mem'} ;
		$mem = $current_command->{'mem'} ;
		$cpu = $current_command->{'cpu'} ;
	}

	# Write command file
	open( my $FH_cmd, ">", $cmd_file_global ) or die "Cannot create ".$cmd_file_global ;
	print $FH_cmd $cmd_file_global_content ;
	close( $FH_cmd );
	open( my $FH_cmd_wrapper, ">", $cmd_file_wrapper ) or die "Cannot create ".$cmd_file_wrapper ;
	print $FH_cmd_wrapper 'head -$SGE_TASK_ID '.$cmd_file_global.' | tail -1 | $SHELL'."\n" ;
	print $FH_cmd_wrapper 'exit $?' ;
	close( $FH_cmd_wrapper );

	# Submit
	my $ressources_opt = SGEScheduler::_get_ressources_opt( $mem, $virtual_mem, $cpu );
	my $nb_command = scalar(@commands) ;
	my $cmd_log = `qsub -V -sync y -t 1-$nb_command $ressources_opt $name_opt -e $cmd_file_global.e'\$JOB_ID'.'\$TASK_ID' -o $cmd_file_global.o'\$JOB_ID'.'\$TASK_ID' $cmd_file_wrapper` ;
	my $job_id = undef ;
	if( $cmd_log =~ /Your job-array (\d+)/ ) { #"Your job-array 5023372.1-2:1 ("test.qarray") has been submitted"
		$job_id = $1 ;
	} else {
		die "Unable to retrieve job_id for command '".$cmd_file_global."'." ;			
	}
		
	# Check status
	my $username = getpwuid($<);
	my $qacct_log = `qacct -o $username -j $job_id 2> /dev/null` ;
	if( $? ){ # Retry after wait
		sleep(30);
		$qacct_log = `qacct -o $username -j $job_id` ;
	}
	my @tasks_logs = split(/\n={5,}\n/, $qacct_log) ;
	foreach my $current_task_log (@tasks_logs) {
		chomp( $current_task_log );
		if( $current_task_log ne "" ) {
			my $failed = undef ;
			if( $current_task_log =~ /\n\s*failed\s+(\d+)\s*\n/ ){
				$failed = int($1) ;
			}
			my $exit_status = undef ;
			if( $current_task_log =~ /\n\s*exit_status\s+(\d+)\s*\n/ ){
				$exit_status = int($1) ;
			}
			if( !defined($failed) || !defined($exit_status) ){
				die "Error with exit status parsing for '".$cmd_file_global."' (job ID: ".$job_id.")" ;	
			}
			if( $failed != 0 || $exit_status != 0 ){
				die "Error in command '".$cmd_file_global."' (job ID: ".$job_id.") see ".$cmd_file_global.".e*" ;	
			}
		}
	}

	# Delete tmp files
	unlink($cmd_file_wrapper);
	unlink($cmd_file_global);
	unlink(glob $cmd_file_global.".e*");
	unlink(glob $cmd_file_global.".o*");
}


=head2 procedure serial_submit

 Usage        : $scheduler->serial_submit( $cmd_set_name, @commands )
 Function     : The implementation of this abstract method is used to submit
                all the commands lines in serial execution mode.
 Args         : [str] The name for the set of command.
                [array] The list of executed commands lines.

=cut
sub serial_submit {
	my ($self, $cmdSet_name, @commands) = @_ ;
	
	# Packs creation
	my @commands_pack = ();
	my $cmd_file_content = "#!/bin/bash\n\n" ;
	my $virtual_mem = "" ;
	my $mem = "" ;
	my $cpu = "" ;
	my $command_idx = 1 ;
	foreach my $current_command (@commands) {
		if( $cmd_file_content ne "#!/bin/bash\n\n" && 
		    (($current_command->{'virtual_mem'} ne $virtual_mem) || ($current_command->{'mem'} ne $mem) || ($current_command->{'cpu'} ne $cpu))
		){
			my %previous_pack = ( 'cmd' => $cmd_file_content,
			                      'cpu' => $cpu,
			                      'mem' => $mem,
			                      'virtual_mem' => $virtual_mem
			); 
			push( @commands_pack, \%previous_pack);
			$cmd_file_content = "#!/bin/bash\n\n" ;
			$command_idx = 1 ;
		}
		$virtual_mem = $current_command->{'virtual_mem'} ;
		$mem = $current_command->{'mem'} ;
		$cpu = $current_command->{'cpu'} ;
		$cmd_file_content .= $current_command->{'cmd'}." ;\nif [ \$? -ne 0 ]; then\n\techo \"The command ".$command_idx." in ##COMMAND_FILE## has failed\" >&2\n\texit 1;\nfi\n\n" ;
		$command_idx++ ;
	}
	if( $cmd_file_content ne "#!/bin/bash\n\n" ){
		my %previous_pack = ( 'cmd' => $cmd_file_content,
		                      'cpu' => $cpu,
		                      'mem' => $mem,
		                      'virtual_mem' => $virtual_mem
		); 
		push( @commands_pack, \%previous_pack);
	}
	
	# Packs execution
	my $pack_idx = 0 ;
	foreach my $current_pack (@commands_pack) {
		my $cmd_file = $self->{'working_dir'}."/tmpCmd_".time."_".int(rand(10000))."_".$pack_idx.".sh" ;
		my $name_opt = "" ;
		if( defined($cmdSet_name) ){
			$cmd_file = $self->{'working_dir'}."/".$cmdSet_name."_".time."_".int(rand(10000)).".sh" ;
			$name_opt = "-N ".$cmdSet_name ;
		}
		$$current_pack{'cmd'} =~ s/##COMMAND_FILE##/$cmd_file/g ;
		
		# Write command file
		open( my $FH_cmd, ">", $cmd_file ) or die "Cannot create ".$cmd_file ;
		print $FH_cmd $$current_pack{'cmd'} ;
		close( $FH_cmd );
		
		# Submit
		my $ressources_opt = SGEScheduler::_get_ressources_opt( $$current_pack{'mem'}, $$current_pack{'virtual_mem'}, $$current_pack{'cpu'} );
		my $cmd_log = `qsub -V -sync y $ressources_opt $name_opt -e $cmd_file.e -o $cmd_file.o $cmd_file` ;
		my $job_id = undef ;
		if( $cmd_log =~ /Your job (\d+)/ ) { #"Your job 5023349 ("test.qsub") has been submitted"
			$job_id = $1 ;
		} else {
			die "Unable to retrieve job_id for command '".$cmd_file."'." ;			
		}
		
		# Check status
		my $username = getpwuid($<);
		my $qacct_log = `qacct -o $username -j $job_id 2> /dev/null` ;
		if( $? ){ # Retry after wait
			sleep(30);
			$qacct_log = `qacct -o $username -j $job_id` ;
		}
		my $failed = undef ;
		if( $qacct_log =~ /\n\s*failed\s+(\d+)\s*\n/ ){
			$failed = int($1) ;
		}
		my $exit_status = undef ;
		if( $qacct_log =~ /\n\s*exit_status\s+(\d+)\s*\n/ ){
			$exit_status = int($1) ;
		}
		if( !defined($failed) || !defined($exit_status) ){
			die "Error with exit status parsing for '".$cmd_file."' (job ID: ".$job_id.")" ;
		}
		if( $failed != 0 || $exit_status != 0 ){
			die "Error in command '".$cmd_file."' (job ID: ".$job_id.") see ".$cmd_file.".e" ;
		}

		# Delete tmp files
		unlink($cmd_file);
		unlink($cmd_file.".e");
		unlink($cmd_file.".o");
		
		$pack_idx++ ;
	}
}


1;