# $Id$

=pod
=head1 NAME

FastaIO - Manage I/O file in fasta format.

=head1 USAGE

# Write Fasta
my $w = new FastaIO( "test.fasta.gz", ">" );
$w->write( new Sequence("id_1", "ATGC", "desc 1") );
$w->write( new Sequence("id_2", "ATGC") );
$w->write( new Sequence("id_3", "ATGC", "desc 3") );
$w->close();

# Read Fasta
my $r = new FastaIO( "test.fasta.gz" );
while( my $record = $r->next() ){
	...
}
$r->close();

=head1 VERSION

 1.0.2

=head1 AUTHORS

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package FastaIO ;
use strict ;
use warnings ;
use FileIO ;
use Sequence ;

our @ISA = ("FileIO");


=head2 function next

 Usage        : $seq_record = $Fasta_obj->next()
 Function     : Returns the next sequence in the fasta file.
 Returns      : [Sequence] The next sequence.
 Args         : none

=cut
sub next {
	my ($self) = @_ ;
	my $FH = $self->{'file_handle'} ;
	
	if( eof($FH) ){
		return undef ;
	} else {
		my $next_line = <$FH> ;
		# First line of file
		if( !defined($self->{'next_header'}) ){
			chomp($next_line);
			$self->{'next_header'} = $next_line ;
			$next_line = <$FH> ;
		}
		# Sequence seq
		my $seq_str = "" ;
		while( defined($next_line) && !($next_line =~ '^>') ){
			chomp($next_line);
			$seq_str .= $next_line ;
			$next_line = <$FH> ;
		}
		# Sequence attributes	
		my @fields = split( /\s/, $self->{'next_header'} );
		my $seq_id = shift( @fields );
		$seq_id =~ s/^>// ;
		my $seq_desc = undef ;
		if( scalar(@fields) > 0 ){
			$seq_desc = join( " ", @fields );
		}
		# Store next seq header
		if( !eof($FH) ){
			chomp($next_line);
			$self->{'next_header'} = $next_line ;
		}
		
		return new Sequence( $seq_id, $seq_str, $seq_desc );
	}
}

=head2 function getAll

 Usage        : $seq_record = $Fasta_obj->getAll()
 Function     : Returns all sequences in the fasta file.
 Returns      : [hash] The hash of sequences.
 Args         : none

=cut

sub getAll {
	my ($self) = @_ ;
	my %sequences;
	my $FH = $self->{'file_handle'} ;

	while (	my $next_line = <$FH> ) {
		# First line of file
		if( !defined($self->{'next_header'}) ){
			chomp($next_line);
			$self->{'next_header'} = $next_line ;
			$next_line = <$FH> ;
		}
		# Sequence seq
		my $seq_str = "" ;
		while( defined($next_line) && !($next_line =~ '^>') ){
			chomp($next_line);
			$seq_str .= $next_line ;
			$next_line = <$FH> ;
		}
		# Sequence attributes	
		my @fields = split( /\s/, $self->{'next_header'} );
		my $seq_id = shift( @fields );
		$seq_id =~ s/^>// ;
		my $seq_desc = undef ;
		if( scalar(@fields) > 0 ){
			$seq_desc = join( " ", @fields );
		}
		# Store next seq header
		if( !eof($FH) ){
			chomp($next_line);
			$self->{'next_header'} = $next_line ;
		}
		$sequences{$seq_id} = new Sequence( $seq_id, $seq_str, $seq_desc );
	}
	return %sequences;
}

=head2 procedure write

 Usage        : $Fasta_obj->write( $seq_record )
 Function     : Writes the sequence in the fasta file.
 Args         : [Sequence] The sequence written.

=cut
sub write {
	my ($self, $sequence_record) = @_ ;
	my $FH = $self->{'file_handle'} ;
	print $FH FastaIO->seqToFastaLine($sequence_record)."\n" ;
}


=head2 function seqToFastaLine

 Usage        : $FastaIO->seqToFastaLine( $seq_record )
 Function     : Returns the sequence in fasta format.
 Returns      : [str] The sequence in fasta format.
 Args         : [Sequence] The sequence processed.

=cut
sub seqToFastaLine {
	my ($class, $sequence) = @_ ;
	my $string = ">".$sequence->{'id'} ;
	if( defined($sequence->{'desc'}) ){
		$string .= " ".$sequence->{'desc'} ;
	}
	$string .= "\n".$sequence->{'seq'} ;
	return $string ;
}


1;