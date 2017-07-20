#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 fastaMetrics2tsv.pl

=head1 SYNOPSIS

 fastaMetrics2tsv.pl FASTA_FILE_NAME FASTA_METRICS [--help] > TSV_OUT

=head1 DESCRIPTION

 Extract nb sequences from contig metrics.

=head1 OPTIONS

=over 8

=item B<INPUT>

 Fasta file name.
 Output file from fastaAssemblyMetrics.pl.

=back

=head1 VERSION

 0.1

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
use Pod::Usage ;
use Getopt::Long ;


MAIN:
{
	# Parameters
	my $name = $ARGV[0]||undef ;
	my $file = $ARGV[1]||undef ;
	my $help = 0 ;
	GetOptions(
		'h|help'   => \$help,
	);
	
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: 'FILE NAME' is required.") if !defined($name) ;
	pod2usage("$0: 'MECTRICS FILE' is required.") if !defined($file) ;
	pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
	pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;

	# Process
	my %metrics = parse_fasta_metrics( $file );
	
	print join("\t", $name, $metrics{nb_seq}, 0)."\n";
}


#=head2 function parse_fasta_metrics
#
# Usage        : %metrics = parse_fasta_metrics( $log )
# Function     : Returns the assembly metrics from fastaAssemblyMetrics.pl 
#                output.
# Returns      : [hash] The metrics.
# Args         : [str] fastaAssemblyMetrics.pl output file.
#
#=cut
sub parse_fasta_metrics {
	my ( $input ) = @_ ;
	my %metrics = ();
	
	my $section = "" ;
	my @keys = () ;
	open(my $FH, "<", $input);
	while( my $line = <$FH> ) {
		chomp( $line );
		if( $line =~ /^##/ ) {
			$section = substr($line, 2);
		} elsif( $line =~ /^#/ ) {
			@keys = split(/\t/, substr($line, 1));
		} elsif( scalar(@keys) > 0 ) {
			my @values = split(/\t/, $line);
			for( my $idx = 0 ; $idx < scalar(@keys) ; $idx++ ){
				if( $section eq "Length" ) {
					if( !exists($metrics{'length_distribution'}) ) {
						my %hash = ();
						$metrics{'length_distribution'} = \%hash;
					}
					$metrics{'length_distribution'}{lc($keys[$idx])} = sprintf("%.2f", $values[$idx])*1 ;
				} else {
					$metrics{lc($keys[$idx])} = sprintf("%.2f",  $values[$idx])*1 ;
				}
			}
			@keys = ();
		}
	}
	close( $FH );

	return %metrics ;
}