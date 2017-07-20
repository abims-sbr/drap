#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 coverageMetrics2json.pl

=head1 SYNOPSIS

 coverageMetrics2json.pl BAM INPUT [--help] > JSON_OUT

=head1 DESCRIPTION

 Transforms samtools depth output in JSON report.
 
=head1 OPTIONS

=over 8

=item B<BAM>

 The ordered and filtered BAM.

=back

=over 8

=item B<INPUT>

 The output file from samtools depth.

=back

=head1 VERSION

 0.1.0

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
use File::Basename;
use List::Util ;
use Pod::Usage ;
use Getopt::Long ;
use JSON ;


MAIN:
{
	# Parameters
	my $sam_depth_out = scalar(@ARGV) > 0 ? $ARGV[0] : undef ;
	my $bam = scalar(@ARGV) > 1 ? $ARGV[1] : undef ;
	my @analysed_coverage = ( 60, 70, 80, 90 );
	my $help = 0 ;
	GetOptions(
		'h|help'  => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($sam_depth_out) ;
	pod2usage("$0: '".$sam_depth_out." doesn't exist.") if !(-e $sam_depth_out) ;
	pod2usage("$0: '".$sam_depth_out." is not readable.") if !(-r $sam_depth_out) ;
	pod2usage("$0: 'BAM' is required.") if !defined($bam) ;
	pod2usage("$0: '".$bam." doesn't exist.") if !(-e $bam) ;
	pod2usage("$0: '".$bam." is not readable.") if !(-r $bam) ;
	
	# Process
	my $contigs_sizes = get_ref_region_lengths( $bam, $sam_depth_out );
	my %metrics = parse_log( $sam_depth_out, \@analysed_coverage, $contigs_sizes );

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=head2 function get_ref_region_lengths
#
# Usage        : %lengths = get_ref_region_lengths( $bam, $out_file )
# Function     : Returns the contigs lengths.
# Returns      : [hash] contig length by contig name.
#                Example : {
#                  "contig_A001": 10,
#                  "contig_A002": 15
#                }
# Args         : [str] path to BAM file.
#                [ref] path to the program output file (used for temp file).
#
#=cut
sub get_ref_region_lengths {
	my ($bam, $out_file) = @_ ;
	my %contigs_sizes = () ;
	
	# Tmp file
	my ($out_filename, $out_dir, $out_ext) = fileparse($out_file, qr/\.[^.]*/);
	my $tmp_bam_header = "tmp_".time."_".int(rand(10000)) ;
	if( $out_dir ne "" ){
		$tmp_bam_header = $out_dir."/".$tmp_bam_header ;
	}
	
	# Retrieve lengths
	`samtools view -H $bam > $tmp_bam_header` ;
		$? and die "Cannot parse '".$bam."'\n" ;
	open( my $FH_BAM_HEADER, "<", $tmp_bam_header );
	while( my $line = <$FH_BAM_HEADER> ){ # @SQ    SN:ORENI.1.1    LN:1745
		chomp($line);
		if( $line =~ /^\@SQ\s+SN\:([^\s]+)\s+LN\:(\d+)$/ ){
			$contigs_sizes{$1} = int($2) ;
		}
	}
	close( $FH_BAM_HEADER );
	unlink( $tmp_bam_header );
	
	return \%contigs_sizes ;
}


#=head2 function parse_log
#
# Usage        : %metrics = parse_log( $depth_out, $analysed_coverage, $contigs_lengths )
# Function     : Returns the coverage metrics.
# Returns      : [hash] "depth_partial_coverage" - The number of contig by depth for each analysed coverage.
#                "depth_histogram" - The number of contigs nt by depth.
#                Example : {
#                  "depth_partial_coverage": {
#                      "60": {          
#                          "1": 10,       # It exists 10 contigs with 60% of the contig covered by at list 1 read
#                          "20": 2        # It exists 2 contigs with 60% of the contig covered by at list 20 reads
#                      }
#                  },
#                  "depth_histogram": {
#                      "1": 100,          # 100 nt in contigs are supported by 1 read.
#                      "20": 2000         # 2000 nt in contigs are supported by 20 reads.
#                  }
#                }
# Args         : [str] samtools depth log file.
#                [ref] list of coverage evaluated.
#                [ref] hash of contig length by contig name.
#
#=cut
sub parse_log {
	my ($step_file, $analysed_coverage, $contigs_lengths) = @_ ;
	my %coverage_hist = () ;
	my %depth_partial_coverage = () ;
	foreach my $current_analysed_cov (@$analysed_coverage){
		my %empty_hash = ();
		$depth_partial_coverage{$current_analysed_cov} = \%empty_hash ;
	}
	
	# Get data
	my $prev_contig_name = "" ;
	my @coverage = ();
	my $prev_contig_pos = 0 ;
	open(my $FH_input, "<", $step_file) or die "Can't open ".$step_file." for read: ".$! ;
	while (my $line = <$FH_input>) {
		chomp($line);
		my ($contig_name, $contig_pos, $depth) = split( /\t/, $line );
		if( $prev_contig_name ne "" and $prev_contig_name ne $contig_name ){
			# Complete missing positions
			while( $prev_contig_pos != $$contigs_lengths{$prev_contig_name} ){
				push( @coverage, 0 );
				$prev_contig_pos++ ;
			}
			# Store metrics
			my @rev_sort_coverage = sort {$b <=> $a} @coverage ; # sort numerically descending
			undef(@coverage);
			my $nb_pos = scalar(@rev_sort_coverage);
			foreach my $current_analysed_cov (@$analysed_coverage){
				my $nb_pos = int( (($nb_pos*$current_analysed_cov)/100) + 0.5 );
				$depth_partial_coverage{$current_analysed_cov}{$rev_sort_coverage[$nb_pos -1]} += 1 ;
			}
			foreach my $current_depth (@rev_sort_coverage){
				$coverage_hist{$current_depth} += 1 ;
			}
			# Next contig
			@coverage = ();
			$prev_contig_pos = 0 ;
		}
		while( $prev_contig_pos != ($contig_pos - 1) ){ # complete missing positions
			push( @coverage, 0 );
			$prev_contig_pos++ ;
		}
		push( @coverage, $depth );
		$prev_contig_pos = $contig_pos ;
		$prev_contig_name = $contig_name ;
	}
	if( $prev_contig_name ne "" ){
		# Complete missing positions
		while( $prev_contig_pos != $$contigs_lengths{$prev_contig_name} ){
			push( @coverage, 0 );
			$prev_contig_pos++ ;
		}
		# Store metrics
		my @rev_sort_coverage = sort {$b <=> $a} @coverage ; # sort numerically descending
		undef(@coverage);
		my $nb_pos = scalar(@rev_sort_coverage);
		foreach my $current_analysed_cov (@$analysed_coverage){
			my $nb_pos = int( (($nb_pos*$current_analysed_cov)/100) + 0.5 );
			$depth_partial_coverage{$current_analysed_cov}{$rev_sort_coverage[$nb_pos -1]} += 1 ;
		}
		foreach my $current_depth (@rev_sort_coverage){
			$coverage_hist{$current_depth} += 1 ;
		}
	}
	close( $FH_input );
	
	return (
		'depth_partial_coverage' => \%depth_partial_coverage,
		'depth_histogram'        => \%coverage_hist
	);
}