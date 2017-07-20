#!/usr/bin/perl -w

=head1 Description

	read a FASTA file as STDIN
	if parameter 1 is a file, consider it as a list of accession numbers
	otherwise consider all paramters as accession numbers of sequences
	to be extrated from the FASTA file

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2017 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict;
use warnings;
use FindBin ;
use lib ("$FindBin::Bin");
use FastaIO;

($#ARGV >= 0) || &usage("wrong number of parameters");

my @acc;
if (open(FIN, $ARGV[0])) { # file or command providing a list of AC
	while (<FIN>)	{
		chomp;
		push(@acc, $_);
	}
	close(FIN);
} else { # AC to be extracted are given as parameters
	@acc = @ARGV;
}

my $inputFastaFh = new FastaIO("STDIN");
my %fastaSeqs = $inputFastaFh->getAll();
my $outputFastaFh = new FastaIO("STDOUT", ">");
foreach my $acc (@acc) {
	$outputFastaFh->write($fastaSeqs{$acc});
}

sub usage( $ ) {
		print STDERR "$_[0]\n";
		system("pod2text $0");
		exit(-1);
}
