#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 expressMetrics.pl

=head1 SYNOPSIS

 expressMetrics.pl --title INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms express log in JSON report.
 
=head1 OPTIONS

=over 8

=item B<INPUT>

 The log file from express.

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
		'h|help'    => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Process
	my %fpkm_stat = parse_log( $file );
	my %metrics = (
		'fpkm_distribution' => \%fpkm_stat
	);
	
	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=head2 function parse_log
#
# Usage        : %metrics = parse_log( $log )
# Function     : Returns the FPKM distribution metrics.
# Returns      : [hash] The FPKM min, lower quartile, mean, median, upper 
#                quartile, max.
# Args         : [str] express log file.
#
#=cut
sub parse_log {
	my ($step_file) = @_ ;
	my @sorted_fpkm = () ;
	my $nb_seq = 0 ;
	my $fpkm_sum = 0 ;
	
	# Get data
	my @fpkm = () ;
	open(my $FH_input, "<", $step_file) or die "Can't open ".$step_file." for read: ".$! ;
	my $header = <$FH_input> ;
	while (my $line = <$FH_input>) {
		chomp($line);
		$nb_seq++ ;
		my @line_fields = split(/\t/, $line) ;
		my $current_fpkm = sprintf( "%.2f", $line_fields[10] )*1 ;
		push( @fpkm, $current_fpkm );
		$fpkm_sum += $current_fpkm ;
	}
	close( $FH_input );
	@sorted_fpkm = sort {$a <=> $b} @fpkm ;
	undef( @fpkm );

	# FPKM stat
	return (
		'min'            => $sorted_fpkm[0],
		'lower_quartile' => $sorted_fpkm[int(($nb_seq/4) + 0.5) -1],
		'mean'           => sprintf("%.2f", $fpkm_sum/$nb_seq)*1,
		'median'         => median( \@sorted_fpkm ),
		'upper_quartile' => $sorted_fpkm[int((3*($nb_seq/4)) + 0.5) -1],
		'max'            => $sorted_fpkm[$nb_seq - 1]
	);
}

#=head2 function median
#
# Usage        : $med = median( $ref_sorted_list )
# Function     : Returns the median.
# Returns      : [int] The median value.
# Args         : [Array ref] Ascending values.
#
#=cut
sub median
{
    my ( $sorted_list ) = @_ ;
    my $nb_seq = scalar( @$sorted_list );
	if( $nb_seq%2 ) {
        return $$sorted_list[($nb_seq/2) -1] ;
	} else {
		return sprintf("%.2f", ($$sorted_list[($nb_seq/2) -1] + $$sorted_list[($nb_seq/2)])/2)*1 ;
	}
}