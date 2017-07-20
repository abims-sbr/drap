#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 blat2json.pl

=head1 SYNOPSIS

 blat2json.pl INPUT [--maxoverlap INT --minidentity FLOAT --mincoverage FLOAT 
                     --help] > JSON_OUT
=head1 DESCRIPTION

 Transforms endBlat.pl output in JSON report.
 
=head1 OPTIONS

=over 8

=item <INPUT>

 The output file from endBlat.pl.

=back

=item maxoverlap INT

 Several proteins can have the same match area on one contig. These proteins 
 as merged in count when the match on the contig is totaly included 
 on the other. Or when the non-overlapped region concern $max_out_overlap 
 nucleotids. [Defaul: 10]

=back

=item minidentity FLOAT

 Minimum identity ratio to use match in process [Default 0.8].

=back

=item mincoverage FLOAT

 Minimum query coverage ratio to use match in process [Default 0.8].

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
	my $max_out_overlap = 10 ;
	my $min_coverage = 0.8 ;
	my $min_identity = 0.8 ;
	my $help = 0 ;
	GetOptions(
		'o|maxoverlap=i'  => \$max_out_overlap,
		'i|minidentity=f' => \$min_identity,
		'c|mincoverage=f' => \$min_coverage,
		'h|help'           => \$help,
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'INPUT' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Process
	my %metrics = parse_out( $file, $max_out_overlap, $min_coverage, $min_identity );

	print to_json( \%metrics, { utf8 => 1, pretty => 1 } );
}


#=item function parse_out
#
# Usage        : %metrics = parse_out( $filepath, [$max_out_overlap, $min_coverage, $min_identity] )
# Function     : Returns metrics on protein distribution.
# Returns      : [hash] Metrics nb_diff_proteins, nb_contig_with_prot,
#                nb_prot_by_contig.
# Args         : [str] endBlat output file.
#                [int] Several proteins can have the same match area on 
#                      one contig. These proteins as merged in count 
#                      when the match on the contig is totaly included 
#                      on the other. Or when the non-overlapped region 
#                      concern $max_out_overlap nucleotids. [Defaul: 10]
#                [int] Minimum query coverage ratio to use match in process 
#                      [default: 0.8]
#                [int] Minimum identity ratio to use match in process 
#                      [default: 0.8]
#
#=cut
sub parse_out {
	my ($step_file, $max_out_overlap, $min_coverage, $min_identity) = @_ ;
	if( !defined($max_out_overlap) ){
		$max_out_overlap = 10 ;
	}
	if( !defined($min_coverage) ){
		$min_coverage = 0.8 ;
	}
	if( !defined($min_identity) ){
		$min_identity = 0.8 ;
	}

	# Get data
	my %contigs_matches = () ;
	my %proteins = () ;
	open(my $FH_input, "<", $step_file) or die "Can't open ".$step_file." for read: ".$! ;
	my $header = <$FH_input> ; # qName	score	qStart	qEnd	qStrand	qLength	%qCoverage	%identity	%similarity	tName	tStart	tEnd	tStrand	tLength	%tCoverage	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert
	while (my $line = <$FH_input>) {
		chomp($line);
		my @fields = split("\t", $line);
		my $prot_name = $fields[0] ;
		my $score = $fields[1] ;
		my $query_coverage = substr($fields[6], 0, -1);
		my $identity = substr($fields[7], 0, -1);
		my $contig_strand = substr($fields[4], -1) ;
		my $contig_name = $fields[9] ;
		my $contig_start = $fields[10] ;
		my $contig_end = $fields[11] ;
		if( $query_coverage >= ($min_coverage*100) && $identity >= ($min_identity*100) ){
			if( !exists($contigs_matches{$contig_name}) ) {
				my @init = ();
				$contigs_matches{$contig_name} = \@init ;
			}
			my %current_match = (
				'start'  => $contig_start,
				'end'    => $contig_end,
				'strand' => $contig_strand,
				#~ 'query'  => $prot_name
			);
			$proteins{$prot_name} = 1 ;
			push( @{$contigs_matches{$contig_name}}, \%current_match );
		}
	}
	close( $FH_input );

	# Process count
	my $nb_diff_proteins = scalar(keys %proteins) ;
	my $nb_contig_with_prot = 0 ;
	my %nb_contig_by_nb_prot = ();
	#~ my %nb_prot_without_overlap = ();
	foreach my $contig_name (keys %contigs_matches) {
		$nb_contig_with_prot++ ;
		# Remove match included in other match
		my @matches = @{$contigs_matches{$contig_name}} ;
		for( my $current_idx = 0 ; $current_idx < scalar(@matches) ; $current_idx++ ){
			for( my $compared_idx = 0 ; $compared_idx < scalar(@matches) ; $compared_idx++ ){
				if( ($current_idx != $compared_idx) && ($matches[$current_idx]{'strand'} eq $matches[$compared_idx]{'strand'}) ){
					if( ($matches[$current_idx]{'start'}-$max_out_overlap <= $matches[$compared_idx]{'start'}) && ($matches[$current_idx]{'end'}+$max_out_overlap >= $matches[$compared_idx]{'end'}) ){ # The current match is downstream
						$matches[$current_idx]{'start'} = $matches[$current_idx]{'start'} < $matches[$compared_idx]{'start'} ? $matches[$current_idx]{'start'} : $matches[$compared_idx]{'start'} ;
						$matches[$current_idx]{'end'} = $matches[$current_idx]{'end'} > $matches[$compared_idx]{'end'} ? $matches[$current_idx]{'end'} : $matches[$compared_idx]{'end'} ;
						splice(@matches, $compared_idx, 1);
						if( $compared_idx < $current_idx ){
							$current_idx -= 1 ;
						}
						$compared_idx -= 1 ;
					}
				}
			}
		}
		my $nb_prot = scalar(@matches);
		if( !exists($nb_contig_by_nb_prot{$nb_prot}) ){
			$nb_contig_by_nb_prot{$nb_prot} = 0 ;
		}
		$nb_contig_by_nb_prot{$nb_prot}++ ;
	}
	
	return (
		"nb_diff_proteins"     => $nb_diff_proteins,
		"nb_contig_with_prot"  => $nb_contig_with_prot,
		"nb_contig_by_nb_prot" => \%nb_contig_by_nb_prot
	);
}