#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 transrateMetrics2json.pl

=head1 SYNOPSIS

 transrateMetrics2json.pl INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms transrate log in JSON report.
 
=head1 OPTIONS

=over 8

=item B<INPUT>

 The log file from transrate STDOUT.

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
	my ($score, $optimal_score, $optimal_cutoff, $contigs_log_path) = parse_log( $file );
	$metrics{'score'} = $score ;
	$metrics{'optimal_score'} = $optimal_score ;
	$metrics{'optimal_cutoff'} = $optimal_cutoff ;
	# Score by contig
	if( !defined($contigs_log_path) ){
		die "The file containing contig metrics for each contig cannot be retrieved from ".$file."\n" ;
	}
	$metrics{'contigs_scores'} = parse_contigs_log( $contigs_log_path );

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


=head2 function parse_log

 Usage        : ($score, $optimal_score, $optimal_cutoff, 
                 $contigs_scores_path) = parse_log( $log )
 Function     : Returns assembly score, assembly optimal score, assembly 
                optimal cutoff, and path to the file containing the contig 
                metrics for each contig.
 Returns      : [array] assembly score, assembly optimal score, assembly 
                optimal cutoff, and path to the file containing the contig 
                metrics for each contig.
 Args         : [str] transrate log file.

=cut
sub parse_log {
	my ($log_file) = @_ ;
	my $score = undef ;
	my $optimal_score = undef ;
	my $optimal_cutoff = undef ;
	my $contigs_scores_path = undef ;
	
	# Get global scores
	open(my $FH_input, "<", $log_file) or die "Can't open ".$log_file." for read: ".$! ;	
	while (my $line = <$FH_input>) {
		chomp($line);
		# Content example:
		#[ INFO] 2015-09-23 11:02:04 : TRANSRATE ASSEMBLY SCORE     0.0077
		#[ INFO] 2015-09-23 11:02:04 : -----------------------------------
		#[ INFO] 2015-09-23 11:02:04 : TRANSRATE OPTIMAL SCORE      0.0438
		#[ INFO] 2015-09-23 11:02:04 : TRANSRATE OPTIMAL CUTOFF     0.2825
		#[ INFO] 2015-09-23 11:02:04 : good contigs                     25
		#[ INFO] 2015-09-23 11:02:04 : p good contigs                 0.28
		#[ INFO] 2015-09-23 11:02:04 : Writing contig metrics for each contig to path/to/contigs.csv
		#[ INFO] 2015-09-23 11:02:04 : Writing analysis results to assemblies.csv
		if( $line =~ /TRANSRATE ASSEMBLY SCORE\s+(.+)$/ ){
			$score = $1*1 ;
		} elsif( $line =~ /TRANSRATE OPTIMAL SCORE\s+(.+)$/ ){
			$optimal_score = $1*1 ;	
		} elsif( $line =~ /TRANSRATE OPTIMAL CUTOFF\s+(.+)$/ ){
			$optimal_cutoff = $1*1 ;
		} elsif( $line =~ /Writing contig metrics for each contig to (.+)$/ ){
			$contigs_scores_path = $1 ;
		}
	}
	close( $FH_input );

	return ($score, $optimal_score, $optimal_cutoff, $contigs_scores_path);
}

=head2 function parse_contigs_log

 Usage        : %scores_disribution = parse_contigs_log( $contigs_log )
 Function     : Returns the contigs score distribution (quartiles).
 Returns      : [hash ref] The distribution.
 Args         : [str] contigs log file.

=cut
sub parse_contigs_log {
	my ( $contigs_log_path ) = @_ ;
	my $nb_contig = 0 ;
	my $scores_sum = 0 ;
	my @scores = () ;
	my $score_idx = undef ;
	
	open(my $FH_input, "<", $contigs_log_path) or die "Can't open ".$contigs_log_path." for read: ".$! ;	
	while (my $line = <$FH_input>) {
		chomp($line);
		# Content example:
		#contig_name,length,prop_gc,gc_skew,at_skew,cpg_count,cpg_ratio,orf_length,linguistic_complexity_6,has_crb,reference_coverage,hits,in_bridges,p_good,p_bases_covered,p_seq_true,score,p_not_segmented,eff_length,eff_count,tpm,coverage
		#oases_splA_CL18Contig1_1,271,0.413284,-0.107143,-0.144654,20,1.748387,89,0.059814,false,NA,NA,0,0.1875,0.96679,0.944361,0.141565,0.826962,271,32.0079,5810.72,220.75
		#oases_splA_CL1Contig1_1,890,0.385393,0.311953,0.722121,80,2.681733,190,0.091064,false,NA,NA,0,0.761821,1.0,0.900818,0.335377,0.488701,890,590.809,32658.7,1240.7
		my @fields = split( ",", $line );
		if( !defined($score_idx) ){
			($score_idx) = grep { $fields[$_] eq "score" } 0..$#fields ;
		} else {
			$nb_contig++ ;
			my $current_score = $fields[$score_idx]*1 ;
			$scores_sum += $current_score ;
			push( @scores, $current_score );
		}
	}
	close( $FH_input );
	my @sorted_scores = sort {$a <=> $b} @scores ;
	undef( @scores );


	my %scores_distrib = (
		'min'            => $sorted_scores[0],
		'lower_quartile' => $sorted_scores[int(($nb_contig/4) + 0.5) -1],
		'mean'           => sprintf("%.6f",  ($scores_sum/$nb_contig))*1,
		'median'         => median( \@sorted_scores ),
		'upper_quartile' => $sorted_scores[int((3*($nb_contig/4)) + 0.5) -1],
		'max'            => $sorted_scores[$nb_contig - 1]
	);
	return \%scores_distrib ;
}

=head2 function median

 Usage        : $med = median( $ref_sorted_list )
 Function     : Returns the median.
 Returns      : [float] The median value.
 Args         : [Array ref] Ascending values.

=cut
sub median
{
    my ( $sorted_list ) = @_ ;
    my $nb_seq = scalar( @$sorted_list );
	if( $nb_seq%2 ) {
        return $$sorted_list[($nb_seq/2) -1] ;
	} else {
		return sprintf("%.6f", ($$sorted_list[($nb_seq/2) -1] + $$sorted_list[($nb_seq/2)])/2)*1 ;
	}
}