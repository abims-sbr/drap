#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 alignmentMetrics2json.pl

=head1 SYNOPSIS

 alignmentMetrics2json.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms flagstat output in JSON report.

=head1 OPTIONS

=over 8

=item B<INPUT>

 The output file from flagstat.

=back

=head1 VERSION

 0.3

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
	my %metrics = parse_log( $file );

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
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

#=head2 function parse_log
#
# Usage        : %metrics = parse_log( $log )
# Function     : Returns the alignment metrics (mapped, paired, ...).
# Returns      : [hash] The number alignment metrics.
# Args         : [str] flagstat log file.
#
#=cut
sub parse_log {
	my ($step_file) = @_ ;
	my %aln_metrics = () ;
	
	# Get data
	my $nb_lines = 0 ;
	my $nb_secondary = 0 ;
	my $nb_reads = 0 ;
	open(my $FH_input, "<", $step_file) or die "Can't open ".$step_file." for read: ".$! ;
	while (my $line = <$FH_input>) {
		chomp($line);
		if( $line =~ /^(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)\s*$/ ){
			$nb_lines = int($1) + int($2) ;
		} elsif( $line =~ /^(\d+) \+ (\d+) secondary\s*$/ ){
			$nb_secondary = int($1) + int($2) ;
			$nb_reads = $nb_lines - $nb_secondary;
		} elsif( $line =~ /^(\d+) \+ (\d+) mapped \(/ ){
			$aln_metrics{'mapped'} = int($1) + int($2) - $nb_secondary ;
		} elsif( $line =~ /^(\d+) \+ (\d+) paired in sequencing\s*$/ ){
			$aln_metrics{'paired'} = int($1) + int($2) ;
		} elsif( $line =~ /^(\d+) \+ (\d+) read1\s*$/ ){
			if( int($1) != 0 && $nb_reads != (int($1) + int($2) *2) ){ # paired-end AND single
				$aln_metrics{'nb_read_1'} = $nb_reads - (int($1) + int($2)) ;
			} elsif( int($1) != 0 ){ # paired-end only
				$aln_metrics{'nb_read_1'} = int($1) + int($2) ;
			} else { # single only
				$aln_metrics{'nb_read_1'} = $nb_reads ;
			}
		} elsif( $line =~ /^(\d+) \+ (\d+) read2\s*$/ ){
			$aln_metrics{'nb_read_2'} = int($1) + int($2) ;
		} elsif( $line =~ /^(\d+) \+ (\d+) properly paired \(/ ){
			$aln_metrics{'properly_paired'} = int($1) + int($2) ;
		} elsif( $line =~ /^(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)\s*$/ ){
			$aln_metrics{'mate_on_other_chr'} = int($1) + int($2) ;
		}
	}
	close( $FH_input );
	
	return %aln_metrics ;
}