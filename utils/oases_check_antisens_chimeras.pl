#!/usr/bin/perl -w

use strict;

if (-t STDIN) {
	print <<HELP;
cat exonerate.out.tsv | awk -v cov=90 -v id=90 '\$7>=cov&&\$8>=id' | oases_check_antisens_chimeras.pl
HELP
	exit;
}

my (@tmp, %query, %target, %q, $nb_chimera);
while (<STDIN>) {
	chomp;
	@tmp = split("\t", $_);
	my $strand = $tmp[13];
	my ($start, $stop) = $strand eq '-' ? ($tmp[11], $tmp[10]) : ($tmp[10], $tmp[11]);
	if (exists $q{$tmp[0]} && $q{$tmp[0]} eq $tmp[9]) {
		if (($target{$tmp[9]}->{$tmp[0]}->{'start'} >= $stop || $target{$tmp[9]}->{$tmp[0]}->{'stop'} <= $start) && ($strand ne $target{$tmp[9]}->{$tmp[0]}->{'strand'})) {
			print $target{$tmp[9]}->{$tmp[0]}->{'hit'}."\n";
			print "$_\n";
		}
	}
	$q{$tmp[0]} = $tmp[9];
	($target{$tmp[9]}->{$tmp[0]}->{'start'}, $target{$tmp[9]}->{$tmp[0]}->{'stop'}, $target{$tmp[9]}->{$tmp[0]}->{'strand'}, $target{$tmp[9]}->{$tmp[0]}->{'hit'}) = ($start, $stop, $strand, $_);
}
