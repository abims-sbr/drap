# $Id$

=pod
=head1 NAME

AbstractScheduler - Generic class for schedulers.

=head1 DESCRIPTION

 This class allows to submit commands. Each extended class implement the 
 submission on a specific scheduler (local, SGE, SLURM, ...).

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

package AbstractScheduler ;
use strict ;
use warnings ;
use Cmd ;
use Cwd ;


=head2 function new

 Function     : Creates, and returns a new scheduler.
 Returns      : [AbstractScheduler] The scheduler.
 Args         : [str] The path to the working directory.

=cut
sub new {
	my ($class, $working_dir) = @_ ;
	if( !defined($working_dir) ) { $working_dir = getcwd ; }
	
	my $self = {
		'working_dir' => $working_dir
	};
	bless( $self, $class );

	return $self ;
}


=head2 function _get_type

 Usage        : $type = $cheduler->_get_type()
 Function     : Returns the type of scheduler.
 Returns      : [str] The type of scheduler.
 Args         : None

=cut
sub _get_type {
	my ($self) = @_ ;
	return( substr(lc(ref($self)), 0, -9) ); # "LocalScheduler" become "local" 
}


=head2 function get_hash

 Usage        : %hash_representation = $cheduler->get_hash()
 Function     : Returns an hash representation of the scheduler.
 Returns      : [hash] The hash representation.
 Args         : None

=cut
sub get_hash {
	my ($self) = @_ ;
	
	my %hash = (
		'working_dir' => $self->{working_dir},
		'type'        => $self->_get_type(),
		'cpu'         => $self->{cpu}
	);
	return \%hash ;
}


=head2 function TO_JSON

 Usage        : $json_representation = $scheduler->TO_JSON()
 Function     : Returns a json representation of the scheduler. This is the 
                method used by JSON::encode.
 Returns      : [hash] The json representation.
 Args         : None

=cut
sub TO_JSON {
	my ($self) = @_ ;
	return $self->get_hash() ;
}


=head2 procedure parallel_submit

 Function     : The implementation of this abstract method is used to submit
                all the commands lines in parallel execution mode.

=cut
sub parallel_submit { die "Not implemented yet." ; }


=head2 procedure serial_submit

 Function     : The implementation of this abstract method is used to submit
                all the commands lines in serial execution mode.

=cut
sub serial_submit { die "Not implemented yet." ; }

1;