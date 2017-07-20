# $Id$

=pod
=head1 NAME

ConfigFile - Configuration file

=head1 DESCRIPTION

 ConfigFile is used to create an object from a configuration file.
 
   Example of config file content:
   ............................................................................
   # commentary
   [SECTION1] # commentary
   value1
   value2
                  
   [SECTION2]
   key1=value3 # commentary
   key2=value4
   ............................................................................
   
   Resulting ConfigFile object attributes:
   ............................................................................
   'SECTION1': ['value1', 'value2'],
   'SECTION2': {
      'key1': 'value3',
      'key2': 'value4'
   }
   ............................................................................

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

package ConfigFile ;
use strict ;
use warnings ;

sub new {
	my ($class, $config_path) = @_ ;
	
	my $self = {
		'_config_path' => $config_path
	};
	bless( $self, $class );

	$self->_parse();

	return $self ;
}


=head2 process _parse

 Usage        : $config->_parse()
 Function     : Loads the attributes from the configuration file.
 Args         : None.

=cut
sub _parse {
	my ( $self ) = @_ ;
	
	my $current_section = undef ;
	open(my $FH_config, "<", $self->{_config_path}) or die "The file '".$self->{_config_path}."' cannot be read." ;
	foreach my $line ( <$FH_config> ){
		chomp( $line );
		if( $line !~ /^\#/ ){ # Skip commentary lines
			if( $line =~ /^\[([^\]]+)\]/ ){ # Section
				$current_section = $1 ;
			} elsif( $line =~ /^([^\#]+?)=([^\#]+)/ ){ # The section is an hash: 'Key=value'
				my $key = $1 ;
				my $value = $2 ;
				$key =~ s/^\s+|\s+$//g ;
				$value =~ s/^\s+|\s+$//g ;
				if( $value ne "" ){
					if( !defined($self->{$current_section}) ){
						my %hash = ();
						$self->{$current_section} = \%hash ;
					}
					$self->{$current_section}{$key} = $value ;
				}
			} elsif( $line ne "" ) { # The section is an array: 'value'
				if( $line =~ /^([^\#]+)/ ){
					my $value = $1 ;
					$value =~ s/^\s+|\s+$//g ;
					if( !defined($self->{$current_section}) ){
						my @array = ();
						$self->{$current_section} = \@array ;
					}
					push( @{$self->{$current_section}}, $value );
				}
			}
		}
	}
	close( $FH_config );
}

1;