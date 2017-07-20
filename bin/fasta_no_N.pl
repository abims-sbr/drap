#!/usr/bin/perl -w

=head1 Description

  Read a fasta file with multiple entries as STDIN
  Print records without N as STDOUT
  With -ex option, do not remove records with N at contig extremities but trim them

  cat file1 | fasta_no_N.pl -ex 

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict;

&usage() if (-t STDIN);
my $ex = $ARGV[0] && $ARGV[0] eq '-ex' ? 1 : 0;
$/="\n>";

while (my $entry=<STDIN>) {
	chomp $entry;
	$entry=~s/^>//;
	$entry=~s/(\S+).*?\n//;
	my $name=$1;
	$entry=~s/\n//g;
	next if ($entry =~ /N/ && not $ex);
	$entry=~s/^N+//;
	$entry=~s/N+$//;
	printf (">%s\n%s\n",$name,$entry) unless ($entry =~ /N/);
}

sub usage {
	system("pod2text $0");
	exit(1);
}
