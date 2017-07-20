# $Id$

=pod
=head1 NAME

CmdSet - Commands lines representation.

=head1 DESCRIPTION

 Stores and submits a list of commands lines.

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

package CmdSet ;
use strict ;
use warnings ;
use Cmd ;
use SchedulerFactory ;



=head2 function new

 Usage        : $cmd = new CmdSet([$commands], [$mode], [$name], [$scheduler], $array)
 Function     : Creates, and returns a new commands set. The array attribute supports three
                values: empty string = the CmdSet does not contain any array ;
                'has' = the CmdSet contains one or more array CmdSet ; 'is' =  the CmdSet 
                is an array  CmdSet. An array CmdSet is a CmdSet with an undetermined number
                of atomic jobs.
 Returns      : [CmdSet] The commands set.
 Args         : [array ref] The reference to the array of Cmd.
                [str] The submission mode ("serial" or "parallel").
                [str] The commands set name.
                [AbstractScheduler] The scheduler used to submit commands.
                [int] The max CPU required for subcommand(s) execution.
                [str] The CmdSet array status.

=cut
sub new {
	my ($class, $commands, $mode, $name, $scheduler, $array) = @_ ;
	if( !defined($commands) ) { my @list = (); $commands = \@list; }
	if( !defined($mode) ) { $mode = "serial"; }
	if( !defined($scheduler) ) { $scheduler = SchedulerFactory->instantiate() ; }
	if( !defined($array) ) { $array = '' ; }
	
	my $self = {
		'commands'  => $commands,
		'mode'      => $mode,
		'scheduler' => $scheduler,
		'array'     => $array,
		'_name'     => "",
		'_max_cpu'  => 1,
	};
	bless( $self, $class );

	$self->set_name( $name );
	$self->set_max_cpu();
	if (@{$commands}) {
		foreach my $cmd ( @{$commands} ){
			if( UNIVERSAL::isa($cmd, 'CmdSet') ){
				if ( $cmd->{array} eq 'has' || $cmd->{array} eq 'is' ) { $self->{array} = 'has' ; }
			} else {
				if ( ! UNIVERSAL::isa($cmd, 'Cmd') ){
					die "The command added to the commands set must be a Cmd or a cmdSet." ;
				}
			}
		}
	}
	return $self ;
}


=head2 procedure set_name

 Usage        : $cmdSet->set_name( $name )
 Function     : Change the commands set name.
 Args         : [str] The commands set new name.

=cut
sub set_name {
	my ($self, $name) = @_ ;
	if( defined($name) && ($name =~ /\s/) ){
		die "'".$name."' is an invalid name for a command set. A name must not contain any space." ;
	}
	$self->{_name} = $name ;
}


=head2 procedure set_max_cpu

 Usage        : $cmdSet->set_max_cpu( $cpu )
 Function     : Compute the commands set _max_cpu. The max CPU is the number of CPU 
                required for subcommand(s) execution.
 Args         : [str] The commands set new _max_cpu.

=cut
sub set_max_cpu {
	my ($self) = shift ;
	if (@{$self->{commands}}) {
		if ( $self->{mode} eq 'parallel' ){
			# _max_cpu for parallel execution is n cmd * n cpu/cmd
			if( UNIVERSAL::isa($self->{commands}->[0], 'CmdSet') ){
				$self->{_max_cpu} = $self->{commands}->[0]->{_max_cpu} * scalar(@{$self->{commands}});
			} else {
				$self->{_max_cpu} = $self->{commands}->[0]->{cpu} * scalar(@{$self->{commands}});
			}
		} else {
			foreach my $cmd ( @{$self->{commands}} ){
				if( UNIVERSAL::isa($cmd, 'CmdSet') ){
					if ( $cmd->{_max_cpu} > $self->{_max_cpu} ) { $self->{_max_cpu} = $cmd->{_max_cpu} ; }
				} else {
					if ( $cmd->{cpu} > $self->{_max_cpu} ) { $self->{_max_cpu} = $cmd->{cpu} ; }
				}
			}
		}
	}
}


=head2 procedure add_cmd

 Usage        : $cmdSet->add_cmd( $scmd )
 Function     : Add a command to the end of the commands set.
 Args         : [Cmd] The commands to add.

=cut
sub add_cmd {
	my ($self, $cmd) = @_ ;
	if( UNIVERSAL::isa($cmd, 'CmdSet') ){
		if ( $cmd->{array} eq 'has' || $cmd->{array} eq 'is' ) { $self->{array} = 'has' ; }
	} else {
		if ( ! UNIVERSAL::isa($cmd, 'Cmd') ) {
			die "The command added to the commands set must be a Cmd or a CmdSet." ;
		}
	}
	push( @{$self->{commands}}, $cmd );
	$self->set_max_cpu();
}


=head2 procedure submit

 Usage        : $cmdSet->submit()
 Function     : Submit the commad set on scheduler.
 Args         : None

=cut
sub submit() {
	my ($self) = @_ ;
	
	foreach my $cmd ( @{$self->{commands}} ){
		if( UNIVERSAL::isa($cmd, 'CmdSet') ){
			die "You cannot submit a CmdSet if it at least one of the sub-commands is a CmdSet." ;
		}
	}
	
	if( $self->{mode} eq "serial" ){
		$self->{scheduler}->serial_submit( $self->{_name}, @{$self->{commands}} );
	} else {
		if( UNIVERSAL::isa($self->{scheduler}, 'LocalScheduler') ){
			$self->{scheduler}->parallel_submit( $self->{_name}, $self->{_max_cpu}, $self->{array}, @{$self->{commands}} );
		} else {
			$self->{scheduler}->parallel_submit( $self->{_name}, @{$self->{commands}} );
		}
	}
}


=head2 procedure from_hash

 Usage        : $cmdSet->from_hash(%hash_representation)
 Function     : Set the commands set object from an hash representation of the
                commands set (see get_hash and load_json).
 Returns      : [str] The hash representation.
 Args         : None

=cut
sub from_hash {
	my ($self, $hash) = @_ ;
	
	# Mode
	$self->{mode} = $$hash{mode} ;
	
	# Name
	$self->set_name( $$hash{name} );
	
	# Max CPU
	$self->{_max_cpu} = $$hash{max_cpu} ;
	
	# Array
	$self->{array} = $$hash{array} ;
	
	# Scheduler
	if( defined($$hash{scheduler}) ){
		$self->{scheduler} = SchedulerFactory->instantiate( $$hash{scheduler}{working_dir}, $$hash{scheduler}{type}, $$hash{scheduler}{cpu} );
	}
	
	# Commands
	foreach my $command ( @{$$hash{commands}} ){
		if( defined($$command{mode}) ){ # command is a command set
			my $sub_cmd_set = new CmdSet();
			$sub_cmd_set->from_hash( $command );
			push( @{$self->{commands}}, $sub_cmd_set );
		} else {
			my $sub_cmd = new Cmd();
			$sub_cmd->from_hash( $command );
			push( @{$self->{commands}}, $sub_cmd );
		}
	}
}


=head2 function get_hash

 Usage        : %hash_representation = $cmdSet->get_hash()
 Function     : Returns an hash representation of the commands set.
 Returns      : [hash] The hash representation.
 Args         : None

=cut
sub get_hash {
	my ($self) = @_ ;
	
	my @new_commands = ();
	foreach my $command (@{$self->{commands}}){
		push( @new_commands, $command->get_hash() );
	}
	
	my %hash = (
		'commands'  => \@new_commands,
		'mode'      => $self->{mode},
		'name'      => $self->{_name},
		'scheduler' => $self->{scheduler},
		'max_cpu'   => $self->{_max_cpu},
		'array'     => $self->{array}
	);
	return \%hash ;
}


=head2 function TO_JSON

 Usage        : $json_representation = $cmdSet->TO_JSON()
 Function     : Returns a json representation of the commands set. This is the
                method used by JSON::encode.
 Returns      : [hash] The json representation.
 Args         : None

=cut
sub TO_JSON {
	my ($self) = @_ ;
	return $self->get_hash() ;
}


1;