#!/usr/bin/perl

use Pod::Usage;
use Getopt::Long;

my($opt_help, $opt_man);

GetOptions(
	'help' => \$opt_help,
	'man'  => \$opt_man,
)
or pod2usage( "Try '$0 --help' for more information.");

pod2usage( -verbose => 1 ) if $opt_help;
pod2usage( -verbose => 2 ) if $opt_man;

pod2usage( -msg => 'Wrong number of params') unless ($ARGV[0] && $ARGV[1]);
pod2usage( -msg => 'Only param 1 OR 2 could be -') if ($ARGV[0] eq '-' && $ARGV[1] eq '-');

$first = 1 if ($ARGV[1] eq '-');
$last  = 1 if ($ARGV[0] eq '-');
$ARGV[0]--;
$ARGV[1]--;
while (<STDIN>) {
	chomp;
	@t = split(/\t/,$_);
	if ($first) {
		$tmp = splice(@t, $ARGV[0], 1);
		unshift(@t, $tmp);
	}
	elsif ($last) {
		$tmp = splice(@t, $ARGV[1], 1);
		push(@t, $tmp);
	}
	else {
		$tmp = $t[$ARGV[1]];
		$t[$ARGV[1]] = $t[$ARGV[0]];
		$t[$ARGV[0]] = $tmp;
	}
	print join("\t", @t)."\n";
}

=pod

=head1 NAME

 interchange_cols.pl

=head1 SYNOPSIS

 interchange_cols.pl <column to interchange|-> <-|column to interchange> < file.tsv

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

 Reads a tsv file as STDIN. 
 Print to STDOUT the same tsv file with columns given as first and second params interchanged.

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

tsv column interchange

=head1 EXAMPLE

 cat file.tsv | interchange_cols.pl 2 5 # interchanges columns 2 and 5
 cat file.tsv | interchange_cols.pl - 3 # moves column 3 at the last position
 cat file.tsv | interchange_cols.pl 4 - # moves column 4 at the first position

=cut
