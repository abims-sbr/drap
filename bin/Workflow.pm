# $Id$

=pod
=head1 NAME

Workflow - Workflow management.

=head1 DESCRIPTION

 Creates, submits and follows a workflow of commands sets.

=head1 USAGE

my $workflow = new Workflow( $binpath."/submit.pl", undef, $out_dir );
$workflow->add_component( 'inclusion', get_inclusion_cmdSet() );
$workflow->add_component( 'chimera', get_chimeraSelf_cmdSet() );
$workflow->add_component( 'orf', get_ORF_cmdSet() );
$workflow->submit() ;

=head1 VERSION

 1.0.1

=head1 AUTHORS

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package Workflow ;
use strict ;
use warnings ;
use Cwd 'getcwd' ;
use File::Spec ;
use File::Path qw(make_path remove_tree) ;
use POSIX 'strftime' ;
use SchedulerFactory ;



=head2 function new

 Usage        : $wf = new Workflow( $submit_bin_path[, $json_components, 
                                    $out_dir] )
 Function     : Creates, and returns a new workflow.
 Returns      : [Workflow] The workflow.
 Args         : [str] The path to the JSON CmdSet submitter: submit.pl.
                [array ref] The list of named components. A component is a 
                            CmdSet see as a mono block by the workflow. The
                            rerun, the status and the folders tree are based on
                            these blocks as most smallest units.
                            Example: 
                              \[
                            	{"name": "filtering", "cmd": <CmdSet>},
                            	{"name": "alignment", "cmd": <CmdSet>}
                               ]
                [str] The path to the output directory.

=cut
sub new {
	my ($class, $submit_bin_path, $json_components, $out_dir) = @_ ;
	if( !defined($out_dir) ) { $out_dir = getcwd ; }
	if( !defined($json_components) ) { my @empty_cmpts = () ; $json_components = \@empty_cmpts ; }
	
	my @components = ();
	my $self = {
		'out_dir'  => File::Spec->rel2abs($out_dir),
		'json_components' => \@components,
		'submit_bin_path' => $submit_bin_path
	};
	bless( $self, $class );

	# Add components
	for( my $cpt_idx=0 ; $cpt_idx < scalar(@$json_components) ; $cpt_idx++ ){
		$self->add_component( $$json_components['name'], $$json_components['cmd'] );
	}

	return $self ;
}


=head2 procedure add_component

 Usage        : $wf->add_component( $command, [$component_name], $env )
 Function     : Add the component to the workflow.
 Args         : [str] The CmdSet to add as component.
                [str] The component name used in status, job submission and 
                      folders tree.
                [hash] the commands to set the env execution

=cut
sub add_component {
	my ( $self, $command, $component_name, $env ) = @_ ;
	my $component_name_clean = $component_name ;
	$component_name_clean =~ tr/ ./__/ ;
	my $component_dir = File::Spec->rel2abs( $self->{out_dir}."/".$component_name_clean );
	my $component_json_path = $self->{out_dir}."/".(scalar(@{$self->{json_components}}))."_".$component_name_clean.".json" ;
	my %component = ( 'path'=> $component_json_path, 'status' => "to_do", 'name' => $component_name );
	
	if( -e $component_dir && -e $component_dir."/achieved" ){
		$component{'status'} = "already_processed" ;
		push( @{$self->{json_components}}, \%component );
	} else {
		push( @{$self->{json_components}}, \%component );
		
		# Create/clean out dir
		if( -e $component_dir ){
			remove_tree( $component_dir );
		}
		make_path($component_dir);
		
		my @cmds;
		if (exists $env->{$component_name.'_env'}) {
			# Add the env settings command
			push (@cmds, new Cmd($env->{$component_name.'_env'}));
			#my @cmd_with_env = ( new Cmd($env->{$component_name.'_env'}), $command ); 
			#$submission_cmd = new CmdSet( \@cmd_with_env, "serial", "wf_cmpt_".$component_name_clean, SchedulerFactory->instantiate($component_dir) );
		}	

		# Add the achieved file creation
		push(@cmds, $command, new Cmd("touch ##STEP_DIR##/achieved") ); 
		my $submission_cmd = new CmdSet( \@cmds, "serial", "wf_cmpt_".$component_name_clean, SchedulerFactory->instantiate($component_dir) );
		
		# Write command
		my $component_json_str = JSON->new->utf8->allow_nonref->pretty->convert_blessed->encode( $submission_cmd );
		$component_json_str =~ s/##STEP_DIR##/$component_dir/g ;
		open( my $COMPONENT_SCRIPT_FH, ">", $component_json_path ) or die "Can't create file ".$component_json_path."\n" ;
		print $COMPONENT_SCRIPT_FH $component_json_str ;
		close $COMPONENT_SCRIPT_FH ;
	}
}


=head2 procedure submit

 Usage        : $wf->submit( [$mode], [$scheduler_type] )
 Function     : Execute the components commands lines.
 Args         : [str] The execution mode ("serial" or "parallel").
                [str] The scheduler type (see SchedulerFactory).

=cut
sub submit {
	my ( $self, $mode, $scheduler_type ) = @_ ;
	if( !defined($mode) ) { $mode = "serial" ; }
	if( !defined($scheduler_type) ) { $scheduler_type = "local" ; }
	
	# Set workflow
	print "Workflow construction:\n" ;
	my $workflow_cmd_set = new CmdSet(
		undef,
		$mode,
		"workflow",
		SchedulerFactory->instantiate($self->{out_dir}, $scheduler_type)
	);
	for( my $cmpt_idx=0 ; $cmpt_idx < scalar(@{$self->{json_components}}) ; $cmpt_idx++ ){
		my $component = ${$self->{json_components}}[$cmpt_idx] ;
		my $component_name = defined($$component{'name'}) || $$component{'name'} eq "" ? $$component{'name'} : "step_" + $cmpt_idx ;
		if( $$component{'status'} eq "already_processed" ){
			print "\t".$component_name."\t[Already processed]\n" ;
		} else {
			print "\t".$component_name."\t[To do]\n" ;
			$workflow_cmd_set->add_cmd( new Cmd("perl ".$self->{submit_bin_path}." ".$$component{'path'}) );
		}
	}
	
	# Execute workflow
	my $status_msg = "Workflow failed." ;
	my $exit_status = 1 ;
	eval {
		print "Workflow execution:\n"
			."\tStarted at ".strftime("%d/%m/%Y %H:%M:%S", gmtime)."\n" ;
    	$workflow_cmd_set->submit() ;
    	$status_msg = "Workflow complete." ;
    	$exit_status = 0 ;
    	print "\tEnded at ".strftime("%d/%m/%Y %H:%M:%S", gmtime)."\n" ;
    	1;
	} or do {
		print STDERR $@ ;
		print "\tError at ".strftime("%d/%m/%Y %H:%M:%S", gmtime)."\n" ;
	};
	
	# Write workflow status
	open( my $WORKFLOW_STATUS_FH, ">", $self->{out_dir}."/workflow_status" );
	print $WORKFLOW_STATUS_FH $status_msg."\n" ;
	close( $WORKFLOW_STATUS_FH );
	
	# Exit
	exit $exit_status ;
}


1;