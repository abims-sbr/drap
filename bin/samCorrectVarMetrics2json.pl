#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 samCorrectVarMetrics2json.pl

=head1 SYNOPSIS

 samCorrectVarMetrics2json.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms samCorrectVariation log in JSON report.
 
=head1 OPTIONS

=over 8

=item B<INPUT>

 The log file from samCorrectVariation.

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
	my %metrics = parse_log( $file );

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=head2 function parse_log
#
# Usage        : %metrics = parse_log( $log )
# Function     : Returns indel correction metrics.
# Returns      : [hash] Metrics corrected_contigs, errors_in_reads, 
#                insertion_case, deletion_case, substitution_case,
#                contig_discrodant.
# Args         : [str] samCorrectVariation log file.
#
#=cut
sub parse_log {
	my ($step_file) = @_ ;
	my $nb_corrected_contigs = 0 ;
	my $nb_errors_in_reads = 0 ;
	my $nb_contig_discordant = 0 ;
	my %nb_case = (
		"remove"  => 0,
		"add"     => 0,
		"replace" => 0,
	);
	
	# Get data
	my %corrected_contigs = () ;
	open(my $FH_input, "<", $step_file) or die "Can't open ".$step_file." for read: ".$! ;
	while (my $line = <$FH_input>) {
		chomp($line);
		if( $line =~ /^(.+): (add|remove|replace) [A-Z]+ .+ position \d+ \((\d+)\/(\d+)\)$/ ){ # <CONTIG_NAME>: remove <NUCLEOTID(s)> at position <POSITION> (<NB_READS_WITH_REMOVE>/<NB_READS_TOTAL>)
			$corrected_contigs{ $1 } = 1 ;
			$nb_case{$2} += 1 ;
			$nb_errors_in_reads += $4 - $3 ;
			if( $4 == $3 ){
				$nb_contig_discordant += 1 ;
			}
		}
	}
	close( $FH_input );
	$nb_corrected_contigs = scalar(keys %corrected_contigs) ;
	
	return (
		"corrected_contigs" => $nb_corrected_contigs,
		"errors_in_reads" => $nb_errors_in_reads,
		"insertion_case" => $nb_case{"remove"},
		"deletion_case" => $nb_case{"add"},
		"substitution_case" => $nb_case{"replace"},
		"contig_discrodant" => $nb_contig_discordant,
	);
}