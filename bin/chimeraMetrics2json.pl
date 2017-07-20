#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 chimeraMetrics2json.pl

=head1 SYNOPSIS

 chimeraMetrics2json.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms self-chimera log in JSON report.
 
=head1 OPTIONS

=over 8

=item B<INPUT>

 The log file from self-chimera.

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
		'h|help' => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Process
	my ($nb_removed, $nb_trim, $nb_nt_lost) = parse_log( $file );
	my %metrics = (
		'nb_chimera_removed' => $nb_removed,
		'nb_chimera_trim'    => $nb_trim,
		'nb_nt_lost'         => $nb_nt_lost
	);

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=head2 function parse_log
#
# Usage        : ($nb_removed, $nb_trim, $nb_nt_lost) = parse_log( $log )
# Function     : Returns the number of contigs removed, trim, and the total 
#                number of nt lost.
# Returns      : [array] nb_ctig_removed, nb_ctig_trim, nb_nt_lost.
# Args         : [str] self-chimera log file.
#
#=cut
sub parse_log {
	my ($log_file) = @_ ;
	my $nb_removed = 0 ;
	my $nb_trim = 0 ;
	my $nb_nt_lost = 0 ;
	
	# Get data
	open(my $FH_input, "<", $log_file) or die "Can't open ".$log_file." for read: ".$! ;
	my $first_section = 1 ;
	my $header = <$FH_input> ;
	while (my $line = <$FH_input>) {
		chomp($line);
		if( $line =~ /^\s*$/ ){
			$first_section = 0 ;
		}
		if( $first_section ){
			my ( $contig_id, $contig_length, $chimera_type, $keep_pos ) = split( /\s+/, $line );
			if( $keep_pos eq 'discarded' ){
				$nb_removed++ ;
				$nb_nt_lost += $contig_length ;
			} else {
				$nb_trim++ ;
				my ($start, $stop) = split(/-/, $keep_pos);
				$nb_nt_lost += $contig_length - ($stop - $start +1) ;
			}
		}
	}
	close( $FH_input );

	return ($nb_removed, $nb_trim, $nb_nt_lost) ;
}