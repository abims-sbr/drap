#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 mapOrientMetrics.pl

=head1 SYNOPSIS

 mapOrientMetrics.pl BAM_INPUT [--help] > TSV_METRICS

=head1 DESCRIPTION

 Prints by contigs the number of forward and the number of reverse reads. Only 
 single or first read in pair is porcessed.

=head1 OPTIONS

=over 8

=item B<BAM_INPUT>

 Alignment file in BAM format.

=back

=head1 VERSION

 0.1.1

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
use Pod::Usage ;
use Getopt::Long ;


sub isAligned {
	my ($flag) = @_ ;
	return !(int($flag) & 4) ;
}

sub isPaired {
	my ($flag) = @_ ;
	return (int($flag) & 1) ;
}

sub isR1 {
	my ($flag) = @_ ;
	if( isPaired($flag) ){
		return (int($flag) & 64) ;
	} else{
		return 1 ;
	}
}

sub isFirstInAlignment {
	my ($flag) = @_ ;
	return !(int($flag) & 256) ;
}

sub isSupplementaryAlignment {
	my ($flag) = @_ ;
	return (int($flag) & 2048) ;
}

sub isForward {
	my ($flag) = @_ ;
	return !(int($flag) & 16) ;
}

MAIN:
{
	my %orientation = ();
	
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
	
	# Load data
	open(my $FH_BAM, "samtools view ".$file." |") or die "Unable to open '".$file."'.\n" ;
	while( my $line = <$FH_BAM> ){
		chomp( $line );
		my @line_fields = split( /\t/, $line );
		if( isAligned($line_fields[1]) && isR1($line_fields[1]) && isFirstInAlignment($line_fields[1]) && !isSupplementaryAlignment($line_fields[1]) ){
			my $contig = $line_fields[2] ;
			my $orientation = "forward" ;
			if( !isForward($line_fields[1]) ){
				$orientation = "reverse" ;
			}
		
			if( !defined($orientation{$contig}) ){
				my %contig_orientations = (
					'forward' => 0,
					'reverse' => 0
				);
				$orientation{$contig} = \%contig_orientations ;
			}
			
			$orientation{$contig}{$orientation} += 1 ;
		}
	}
	close( $FH_BAM );
	
	# Output
	print "#Contig\tforward\treverse\n" ;
	foreach my $contig ( sort(keys(%orientation)) ){
		print $contig."\t"
			.$orientation{$contig}{'forward'}."\t"
			.$orientation{$contig}{'reverse'}."\n" ;
	}
}