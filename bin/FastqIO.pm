# $Id$

=pod
=head1 NAME

FastqIO  - Manage I/O file in fastq format.

=head1 DESCRIPTION

# Write Fastq
my $w = new FastqIO( "test.fastq.gz", ">" );
$w->write( new Sequence("id_1", "ATGC", "desc 1", "####") );
$w->write( new Sequence("id_2", "ATGC", undef, "####") );
$w->write( new Sequence("id_3", "ATGC", "desc 3", "####") );
$w->close();

# Read Fastq
my $r = new FastqIO( "test.fastq.gz" );
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

package FastqIO ;
use strict ;
use warnings ;
use FileIO ;
use Sequence ;

our @ISA = ("FileIO");


=head2 function next

 Usage        : $seq_record = $Fastq_obj->next()
 Function     : Returns the next sequence in the fastq file.
 Returns      : [Sequence] The next sequence.
 Args         : none

=cut
sub next {
	my ($self) = @_ ;
	my $FH = $self->{'file_handle'} ;
	
	if( eof($FH) ){
		return undef ;
	} else {
		my $seq_header = <$FH> ;
		chomp($seq_header);
		my @fields = split( /\s/, $seq_header );
		my $seq_id = shift( @fields );
		$seq_id =~ s/^@// ;
		my $seq_desc = undef ;
		if( scalar(@fields) > 0 ){
			$seq_desc = join( " ", @fields );
		}
		my $seq_str = <$FH> ;
		chomp($seq_str);
		my $seq_sep = <$FH> ;
		chomp($seq_sep);
		my $seq_qual = <$FH> ;
		chomp($seq_qual);
		
		return new Sequence( $seq_id, $seq_str, $seq_desc, $seq_qual );
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
	
	while ( my $seq_header = <$FH> ) {
		chomp($seq_header);
		my @fields = split( /\s/, $seq_header );
		my $seq_id = shift( @fields );
		$seq_id =~ s/^@// ;
		my $seq_desc = undef ;
		if( scalar(@fields) > 0 ){
			$seq_desc = join( " ", @fields );
		}
		my $seq_str = <$FH> ;
		chomp($seq_str);
		my $seq_sep = <$FH> ;
		chomp($seq_sep);
		my $seq_qual = <$FH> ;
		chomp($seq_qual);
		
		$sequences{$seq_id} = new Sequence( $seq_id, $seq_str, $seq_desc, $seq_qual );
	}
	return %sequences;
}

=head2 procedure write

 Usage        : $Fastq_obj->write( $seq_record )
 Function     : Writes the sequence in the fastq file.
 Args         : [Sequence] The sequence written.

=cut

sub write {
	my ($self, $sequence_record) = @_ ;
	my $FH = $self->{'file_handle'} ;
	print $FH FastqIO->seqToFastqLine($sequence_record)."\n" ;
}


=head2 function seqToFastqLine

 Usage        : $FastqIO->seqToFastqLine( $seq_record )
 Function     : Returns the sequence in fastq format.
 Returns      : [str] The sequence in fastq format.
 Args         : [Sequence] The sequence processed.

=cut
sub seqToFastqLine {
	my ($class, $sequence) = @_ ;
	my $string = "@".$sequence->{'id'} ;
	if( defined($sequence->{'desc'}) ){
		$string .= " ".$sequence->{'desc'} ;
	}
	$string .= "\n".$sequence->{'seq'} ;
	$string .= "\n+" ;
	$string .= "\n".$sequence->{'qual'} ;
	return $string ;
}


1;