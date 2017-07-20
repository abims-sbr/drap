#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 transdecoderMetrics.pl

=head1 SYNOPSIS

 transdecoderMetrics.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms transdecoder log in JSON report.
 
=head1 OPTIONS

=over 8

=item B<INPUT>

 The log file from transdecoder.

=back

=head1 VERSION

 0.1

=head1 AUTHORS

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict ;
use warnings ;
use List::Util ;
use Pod::Usage ;
use Getopt::Long ;
use JSON ;


MAIN:
{
	# Parameters
	my $file = scalar(@ARGV) > 0 ? $ARGV[0] : undef ;
	my $help = 0 ;
	GetOptions(
		'h|help'  => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Process
	my %nb_ORF_distribution = parse_log( $file );
	my %metrics = (
		'nb_ORF_distribution' => \%nb_ORF_distribution
	);
	
	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=head2 function parse_log
#
# Usage        : %metrics = parse_log( $log )
# Function     : Returns the number of contigs by number of ORF.
# Returns      : [hash] The number of contigs by number of ORF.
# Args         : [str] transdecoder log file.
#
#=cut
sub parse_log {
	my ($step_file) = @_ ;
	my %contigs_by_nb_ORF = () ;
	
	# Get data
	open(my $FH_input, "<", $step_file) or die "Can't open ".$step_file." for read: ".$! ;
	my $header = <$FH_input> ;
	while (my $line = <$FH_input>) {
		chomp($line);
		my ( $contig_id, $nb_ORF, $coords ) = split( /\s+/, $line );
		if( exists($contigs_by_nb_ORF{$nb_ORF}) ){
			$contigs_by_nb_ORF{$nb_ORF} += 1 ;
		} else {
			$contigs_by_nb_ORF{$nb_ORF} = 1 ;
		}
	}
	close( $FH_input );

	# Complete data
	my $max_nb_orf = List::Util::max(keys(%contigs_by_nb_ORF)) ;
	for( my $nb_ORF = 0 ; $nb_ORF <= $max_nb_orf ; $nb_ORF++ ){
		if( !exists($contigs_by_nb_ORF{$nb_ORF}) ){
			$contigs_by_nb_ORF{$nb_ORF} = 0 ;
		}
	}
	
	return %contigs_by_nb_ORF ;
}