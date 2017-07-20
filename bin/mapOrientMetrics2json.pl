#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 mapOrientMetrics2json.pl

=head1 SYNOPSIS

 mapOrientMetrics2json.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms orientation metrics in JSON report.

=head1 OPTIONS

=over 8

=item B<INPUT>

 Output file from maporientMetrics.pl.

=back

=head1 VERSION

 0.1

=head1 AUTHOR

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict ;
use warnings ;
use Pod::Usage ;
use Getopt::Long ;
use JSON ;


MAIN:
{
	# Parameters
	my $file = scalar(@ARGV) > 0 ? $ARGV[0] : undef ;
	my $help = 0 ;
	GetOptions(
		'h|help'   => \$help,
	);
	
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Process
	my %metrics = parse_orient_metrics( $file );
	
	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=head2 function parse_orient_metrics
#
# Usage        : %metrics = parse_orient_metrics( $log )
# Function     : Returns the assembly metrics from mapOrientMetrics.pl 
#                output.
# Returns      : [hash] The metrics.
# Args         : [str] mapOrientMetrics.pl output file.
#
#=cut
sub parse_orient_metrics {
	my ( $input ) = @_ ;
	my %hash = ();
	my %metrics = ();
	$metrics{"contigs_by_ratio"} = \%hash ;
	
	my $section = "" ;
	my @keys = () ;
	open(my $FH, "<", $input);
	while( my $line = <$FH> ) {
		chomp( $line );
		if( !($line =~ /^#/) ) {
			my( $contig, $nb_forward, $nb_reverse ) = split( /\t/, $line );
			my $orientation_ratio = int( sprintf("%.0f", (int($nb_forward)*100)/(int($nb_forward) + int($nb_reverse))) );
			$metrics{"contigs_by_ratio"}{$orientation_ratio} += 1 ;
		}
	}
	close( $FH );

	return %metrics ;
}