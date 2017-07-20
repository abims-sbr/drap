#!/usr/bin/perl -w

=head1 Description

  Read a fasta file with multiple entries (parameter 1)
  Print entry name and entry length as STDOUT

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict;

#------------------------------------------------------------

$/="\n>";

while (my $entry=<STDIN>) {
	chomp $entry;
	$entry=~s/^>//;
	$entry=~s/(\S+).*?\n//;
	my $name=$1;
	$entry=~s/\n//g;
	printf ("%s\t%d\n",$name,length($entry));
}
