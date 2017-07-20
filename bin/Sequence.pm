# $Id$

=pod
=head1 NAME

Sequence - Storage for biological sequence

=head1 DESCRIPTION

 A sequence can has id, sequence, description and quality. 

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

package Sequence ;
use strict ;
use warnings ;


=head2 function new

 Usage        : $seq = new Sequence($id, $sequence[, $description, $quality])
 Function     : Creates, and returns a new sequence.
 Returns      : [Sequence] The sequence.
 Args         : [str] The sequence ID.
                [str] The sequence string (example: "ATGCATG").
                [str] The sequence description.
                [str] The sequence quality (example: "I9I#9II").

=cut
sub new {
	my ($class, $id, $sequence, $description, $quality) = @_ ;
	
	my $self = {
		'id'   => $id,
		'desc' => $description,
		'seq'  => $sequence,
		'qual' => $quality
	};
	bless( $self, $class );

	return $self ;
}

=head2 function length

 Usage        : $length = $record->length()
 Function     : Returns the sequence length.
 Returns      : [int] The length.
 Args         : none

=cut
sub length {
	my ($self) = @_ ;
	return length($self->{'seq'}) ;
}

1;