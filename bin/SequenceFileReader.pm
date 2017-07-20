# $Id$

=pod
=head1 NAME

SequenceFileReader - Sequence file factory.

=head1 USAGE

my $r = SequenceFileReader->instantiate( "test.fastq.gz" );
while( my $record = $r->next() ){
	...
}
$r->close();

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

package SequenceFileReader ;
use strict ;
use warnings ;
use FileIO ;
use FastaIO ;
use FastqIO ;


=head2 function instantiate

 Usage        : $file_obj = SequenceFileReader->instantiate( $filepath )
 Function     : Returns the appropriate sequence file reader.
 Returns      : [FastaIO or FastqIO] The sequence reader.
 Args         : [str] The file path.

=cut
sub instantiate {
	my ($class, $filepath) = @_ ;
	my $FH = new FileIO( $filepath );
	my $first_line = $FH->next() ;
	$FH->close();
	
	if( !defined($first_line) ){
		return new FastqIO( $filepath );
	} elsif( $first_line =~ '@' ){
		return new FastqIO( $filepath );
	} elsif( $first_line =~ '>' ){
		return new FastaIO( $filepath );
	} else {
		die "The file ".$filepath." does not have a valid format for 'SequenceFileReader'." ;
	}
}


1;