#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 buscoMetrics2json.pl

=head1 SYNOPSIS

 buscoMetrics2json.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms busco summary in JSON report.
 
=head1 OPTIONS

=over 8

=item B<INPUT>

 The summary file from busco.

=back

=head1 VERSION

 0.1.0

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
		'h|help' => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Global score
	my %metrics = ();
	my ($complete, $complete_single_copy, $complete_duplicated, $fragmented, $missing, $total, $notation, $lineage) = parse_summary( $file );
	$metrics{'complete'} = $complete ;
	$metrics{'complete_single_copy'} = $complete_single_copy ;
	$metrics{'complete_duplicated'} = $complete_duplicated ;
	$metrics{'fragmented'} = $fragmented ;
	$metrics{'missing'} = $missing ;
	$metrics{'total'} = $total ;
	$metrics{'notation'} = $notation ;
	$metrics{'lineage'} = $lineage ;

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


=head2 function parse_suummary

 Usage        : ($complete_single_copy, $complete_duplicated, $fragmented,
                $missing, $total, $notation) = parse_log( $log )
 Function     : Returns complete_single_copy, complete_duplicated, fragmented
                missing, total and BUSCO notation.
 Returns      : [array] complete_single_copy, complete_duplicated, fragmented
                missing, total and BUSCO notation.
 Args         : [str] busco log file.

=cut
sub parse_summary {
	my ($summary_file) = @_ ;
	my $complete = undef ;
	my $complete_single_copy = undef ;
	my $complete_duplicated = undef ;
	my $fragmented = undef ;
	my $missing = undef ;
	my $total = undef ;
	my $notation = undef ;
	my $lineage = undef ;
		
	# Get metrics
	open(my $FH_input, "<", $summary_file) or die "Can't open ".$summary_file." for read: ".$! ;	
	while (my $line = <$FH_input>) {
		chomp($line);
		# Content example:
		# # BUSCO version is: 2.0 
		# # The lineage dataset is: vertebrata_odb9 (Creation date: 2016-10-21, number of species: 65, number of BUSCOs: 2586)
		# # To reproduce this run: ...
		# #
		# # Summarized benchmarking in BUSCO notation for ...
		# # BUSCO was run in mode: tran
		# 
		# C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:2586
		# 
		# 0	Complete BUSCOs (C)
		# 0	Complete and single-copy BUSCOs (S)
		# 0	Complete and duplicated BUSCOs (D)
		# 0	Fragmented BUSCOs (F)
		# 2586	Missing BUSCOs (M)
		# 2586	Total BUSCO groups searched

		if( $line =~ /The lineage dataset is: (\w+) \(.+\)$/ ){
			$lineage = $1 ;
		} elsif( $line =~ /\s*(C:.+?%\[S:.+?%,D:.+?%\],F:.+?%,M:.+?%,n:\d+)$/ ){
			$notation = $1 ;
		} elsif( $line =~ /\s*(\d+)\s+Complete BUSCOs \(C\)$/ ){
			$complete = $1*1 ;
		} elsif( $line =~ /\s*(\d+)\s+Complete and single-copy BUSCOs \(S\)$/ ){
			$complete_single_copy = $1*1 ;
		} elsif( $line =~ /\s*(\d+)\s+Complete and duplicated BUSCOs \(D\)$/ ){
			$complete_duplicated = $1*1 ;
		} elsif( $line =~ /\s*(\d+)\s+Fragmented BUSCOs \(F\)$/ ){
			$fragmented = $1*1 ;
		} elsif( $line =~ /\s*(\d+)\s+Missing BUSCOs \(M\)$/ ){
			$missing = $1*1 ;
		} elsif( $line =~ /\s*(\d+)\s+Total BUSCO groups searched$/ ){
			$total = $1*1 ;
		}
	}
	close( $FH_input );

	return ($complete, $complete_single_copy, $complete_duplicated, $fragmented, $missing, $total, $notation, $lineage);
}
