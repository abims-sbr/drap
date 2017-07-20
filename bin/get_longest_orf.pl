#!/usr/bin/perl

use strict;
use Pod::Usage;
use Getopt::Long;

my ($fa, $fi, $na, $aa, $stats, $help, $man);

GetOptions(
	'na' => \$na,
	'aa' => \$aa,
	'stats' => \$stats,
	'find=i' => \$fi,
	'f=s' => \$fa,
	'help' => \$help,
	'man' => \$man
);

pod2usage(-exitstatus => 0, -verbose => 99, -sections => "NAME|SYNOPSIS|OPTIONS|DESCRIPTION") if ($help);
pod2usage(-verbose => 2) if ($man);
(-e $fa) || pod2usage("-f is not a file: $fa");
($na || $aa || $stats || defined($fi)) || pod2usage("Please specify output type");

my $find = defined($fi) ? $fi : $na ? 2 : 0;

open(GETORF, "cat $fa | getorf -auto -filter -find $find |") or die($!);
my ($length, $contig, $header, $longest, $seq, $flag);
if ($stats) {
	while (<GETORF>) {
		if (/^>(.+)_\d+ \[(\d+) - (\d+)\].*$/) {
			$length = abs($3-$2+1);
			if ($1 ne $contig) {
				printf("%s\t%d\n",$header,$longest) if ($header);
				($contig,$header,$longest) = ($1,"$1\t$2\t$3",$length);
			}
			else {
				($contig,$header,$longest) = ($1,"$1\t$2\t$3",$length) if ($length > $longest);
			}
		}
	}
	printf("%s\t%d\n",$header,$longest);
}
else {
	while (<GETORF>) {
		if (/^>(.+)_\d+ \[(\d+) - (\d+)\].*$/) {
			$length = abs($3-$2+1);
			if ($1 ne $contig) {
				printf(">%s\n%s",$header,$seq) if ($header);
				($contig,$header,$longest,$flag,$seq) = ($1,"$1#$2-$3",$length,1,'');
			}
			else {
				$length > $longest ? ($contig,$header,$longest,$flag,$seq) = ($1,"$1#$2-$3",$length,1,'') : $flag = 0;
			}
		}
		else {
			$seq .= $_ if ($flag);
		}
	}
	printf(">%s\n%s",$header,$seq);
}
close GETORF;

=head1 NAME

get_longest_orf.pl

=head1 SYNOPSIS

get_longest_orf.pl [-h|options] -f file.fa

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-na>

Write fasta format nucleic acids longest ORFs.

=item B<-aa>

Write fasta format amino acids longest ORFs.

=item B<-stats>

Write tsv format position and length of longest ORFs.

=item B<-find>

Find argument given to the EMBOSS getorf command. See getorf -h for more information.
Overwrite -na or -aa argument.

=item B<-f>

Input fasta file.

=back

=head1 DESCRIPTION
  
Read a fasta file with multiple entries.
Find the longest ORF (region that is free of STOP codons if option -find not defined) with the getorf EMBOSS tool and write output to STDOUT.
In ouput fasta format (-na or -aa), sequence names are concatenated with #<orf_start>-<orf_stop>.
Remove it and keep original names piping output in [sed -e 's/\(>.*\)#.*/\1/'].

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=head1 VERSION

1

=head1 DATE

2013

=head1 KEYWORDS

ORF frame fasta longest

=cut

