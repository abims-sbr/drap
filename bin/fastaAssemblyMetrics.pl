#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 fastaAssemblyMetrics.pl

=head1 SYNOPSIS

 fastaAssemblyMetrics.pl --input FASTA > TSV

=head1 DESCRIPTION

 This program generates descriptive statistics on contigs.

=head1 OPTIONS

=over 8

=item B<-i/--input>

 The contigs sequences file.

=back

=head1 VERSION

 0.2.3

=head1 AUTHORS

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=head1 KEYWORDS

 assembly
 
=cut

use strict ;
use warnings ;
use Pod::Usage ;
use Getopt::Long ;
use FindBin ;
use lib ("$FindBin::Bin") ;
use SequenceFileReader ;


MAIN:
{
	# Parameters
	my $input = undef ;
	my $help = 0 ;
	GetOptions(
		'i|input=s' => \$input,
		'h|help'    => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|KEYWORDS|OPTIONS|VERSION"
	) if ($help) ;
	pod2usage("$0: '--input' is required.") if !defined($input) ;
	pod2usage("$0: '".$input." doesn't exist.") if !(-e $input) ;
	pod2usage("$0: '".$input." is not readable.") if !(-r $input) ;
	
	# Get data
	my $nb_seq = 0 ;
	my $lengths_sum = 0 ;
	my @lengths = () ;
	my $FH_input = SequenceFileReader->instantiate($input);
	while( my $record = $FH_input->next() ){
		$nb_seq++ ;
		push( @lengths, $record->length() );
		$lengths_sum += $record->length() ;
	}
	my @sorted_lengths = sort {$a <=> $b} @lengths ;
	undef( @lengths );
	
	# Lengths stat
	my $min = 0 ;
	my $lower_quartile = 0 ;
	my $mean = 0 ;
	my $median = 0 ;
	my $upper_quartile = 0 ;
	my $max = 0 ;
	my $N50 = 0 ;
	my $L50 = 0 ;
	if( $nb_seq > 0 ){
		$min = $sorted_lengths[0] ;
		$lower_quartile = $sorted_lengths[int(($nb_seq/4) + 0.5) -1] ;
		$mean = $lengths_sum/$nb_seq ;
		$median = median( \@sorted_lengths );
		$upper_quartile = $sorted_lengths[int((3*($nb_seq/4)) + 0.5) -1] ;
		$max = $sorted_lengths[$nb_seq - 1] ;
		$N50 = N50( \@sorted_lengths, $lengths_sum );
		$L50 = L50( \@sorted_lengths, $lengths_sum );
	}
	
	# Write output
	print "##Global\n" ;
	print "#Nb_seq\tN50\tL50\n" ;
	print $nb_seq."\t".$N50."\t".$L50."\n" ;
	print "\n" ;
	print "##Length\n" ;
	print "#Sum\tMin\tLower_quartile\tMean\tMedian\tUpper_quartile\tMax\n" ;
	print $lengths_sum."\t".$min."\t".$lower_quartile."\t".$mean."\t".$median."\t".$upper_quartile."\t".$max."\n" ;
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

#=head2 function N50
#
# Usage        : $n50 = N50( $ref_sorted_lengths, $lengths_sum )
# Function     : Returns the N50.
# Returns      : [int] The N50 value.
# Args         : [Array ref] Ascending lengths.
#                [int] Lengths sum.
#
#=cut
sub N50 {
	my ( $sorted_lengths, $lengths_sum ) = @_ ;
	my $idx = scalar(@$sorted_lengths) -1 ;
	my $current_sum = 0 ;
	while( $current_sum < $lengths_sum/2 ) {
		$current_sum += $$sorted_lengths[$idx] ;
		$idx-- ;
	}
	return $$sorted_lengths[$idx] ;
}

#=head2 function L50
#
# Usage        : $l50 = L50( $ref_sorted_lengths, $lengths_sum )
# Function     : Returns the L50.
# Returns      : [int] The L50 value.
# Args         : [Array ref] Ascending lengths.
#                [int] Lengths sum.
#
#=cut
sub L50 {
	my ( $sorted_lengths, $lengths_sum ) = @_ ;
	my $idx = scalar(@$sorted_lengths) -1 ;
	my $current_sum = 0 ;
	while( $current_sum < ($lengths_sum/2) ) {
		$current_sum += $$sorted_lengths[$idx] ;
		$idx-- ;
	}
	return( scalar(@$sorted_lengths) - $idx );
}