# $Id$

=pod
=head1 NAME

Cmd - Command line representation.

=head1 DESCRIPTION

 The command line object representation.

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

package Cmd ;
use strict ;
use warnings ;



=head2 function new

 Usage        : $cmd = new Cmd( $cmd_line, [$cpu], [$mem], [$virtual_mem] )
 Function     : Creates, and returns a new command.
 Returns      : [Cmd] The command line.
 Args         : [str] The number of CPUs used by the command.
                [str] The memory used by the command.
                [str] The virtual memory used by the command.

=cut
sub new {
	my ($class, $cmd, $cpu, $mem, $virtual_mem) = @_ ;
	if( !defined($cpu) || $cpu eq '') { $cpu = 1 ; }
	if( !defined($mem) ) { $mem = "" ; }
	if( !defined($virtual_mem) ) { $virtual_mem = "" ; }
	
	my $self = {
		'cmd' => $cmd,
		'cpu' => int($cpu),
		'mem' => $mem,
		'virtual_mem' => $virtual_mem
	};
	bless( $self, $class );

	return $self ;
}


=head2 function get_hash

 Usage        : %hash_representation = $cmd->get_hash()
 Function     : Returns an hash representation of the command.
 Returns      : [hash] The hash representation.
 Args         : None

=cut
sub get_hash {
	my ($self) = @_ ;
	my %hash = (
		'cmd' => $self->{cmd},
		'cpu' => $self->{cpu},
		'mem' => $self->{mem},
		'virtual_mem' => $self->{virtual_mem}
	);
	return \%hash ;
}


=head2 procedure from_hash

 Usage        : $cmd->from_hash(%hash_representation)
 Function     : Set the command object from an hash representation of the
                command (see get_hash and load_json).
 Returns      : [str] The hash representation.
 Args         : None

=cut
sub from_hash {
	my ($self, $hash) = @_ ;
	
	foreach my $attr ( keys %$hash ){
		$self->{$attr} = $$hash{$attr};
	}
}


=head2 function TO_JSON

 Usage        : $json_representation = $cmd->TO_JSON()
 Function     : Returns a json representation of the command. This is the 
                method used by JSON::encode.
 Returns      : [hash] The json representation.
 Args         : None

=cut
sub TO_JSON {
	my ($self) = @_ ;
	return $self->get_hash() ;
}


1;