# $Id$

=pod
=head1 NAME

FileIO - Manage I/O file in standard and gziped.

=head1 USAGE

# Write file
my $w = new FileIO( "test.txt", ">" );
$w->write( "test1" );
$w->write( "test2" );
$w->write( "test3" );
$w->write( "test4" );
$w->close();

# Write gzip file
my $w_gz = new FileIO( "test.txt.gz", ">" );
$w_gz->write( "test1" );
$w_gz->write( "test2" );
$w_gz->write( "test3" );
$w_gz->write( "test4" );
$w_gz->close();

# Append gzip file
my $w2_gz = new FileIO( "test2.txt.gz", ">" );
$w2_gz->write( "test1" );
$w2_gz->write( "test2" );
$w2_gz->close();
my $a_gz = new FileIO( "test2.txt.gz", ">>" );
$a_gz->write( "test3" );
$a_gz->write( "test4" );
$a_gz->close();

# Read file
my $r = new FileIO( "test.txt" );
while( my $line = $r->next() ){
	print $line ;
}
$r->close();

# Read gzip file
my $r_gz = new FileIO( "test.txt.gz" );
while( my $line = $r_gz->next() ){
	print $line ;
}
$r_gz->close();

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

package FileIO ;
use strict ;
use warnings ;
use IO::Zlib ;


=head2 function new

 Usage        : $file = new FileIO($filepath[, $mode])
 Function     : Creates, and returns a new file object.
 Returns      : [FileIO] The file object.
 Args         : [str] The file path.
                [str] Mode to open the file ('<', '>', '>>').

=cut
sub new {
	my ($class, $filepath, $mode) = @_ ;
	if( !defined($mode) ){ $mode = "<" }
	
	my $self = {
		'filepath' => $filepath,
		'mode'     => $mode
	};
	bless( $self, $class );
	
	if( $filepath =~ /\.gz$/ ){ # read or write zipped
		my $zlib_mode = "w" ;
		if( $mode eq ">>" ){ $zlib_mode = "a" }
		elsif( $mode eq "<" ){ $zlib_mode = "r" }
		$self->{'file_handle'} = IO::Zlib->new($filepath, $zlib_mode) or die "Unable to open '".$filepath."': ".$! ;
	} elsif ($filepath =~ /^std(in|out)$/i ) { # read or write STDIO
		die "Unable to append to STDOUT\n" if ($mode eq ">>") ;
		my $action = $filepath =~ /stdin/i ? 'read from' : 'write to';
		$mode .= '-';
		open( my $FH, $mode ) or die "Unable to $action '".$filepath."': ".$! ;
		$self->{'file_handle'} = $FH ;
	} else { # read or write unzipped
		open( my $FH, $mode, $filepath ) or die "Unable to open '".$filepath."': ".$! ;
		$self->{'file_handle'} = $FH ;
	}

	return $self ;
}

sub DESTROY {
	my ($self) = @_ ;
	$self->close();
}


=head2 procedure close

 Usage        : $file_obj->close()
 Function     : Closes the file handle.
 Args         : none

=cut
sub close {
	my ($self) = @_ ;
	if( defined($self->{'file_handle'}) ){
		close( $self->{'file_handle'} );
	}
}


=head2 function next

 Usage        : $line = $file_obj->next()
 Function     : Returns the next line in the file.
 Returns      : [str] The next line.
 Args         : none

=cut
sub next {
	my ($self) = @_ ;
	my $FH = $self->{'file_handle'} ;
	return( <$FH> );
}


=head2 procedure write

 Usage        : $file_obj->write( $line )
 Function     : Writes the line in the fastq file.
 Args         : [str] The line written.

=cut
sub write {
	my ($self, $text) = @_ ;
	my $FH = $self->{'file_handle'} ;
	print $FH $text."\n" ;
}


1;