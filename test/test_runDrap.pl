#!/usr/bin/perl 

# Usage: 
# perl test_runDrap.p
# perl test_runDrap.pl --local
# Launch runDrap with combinations of options in @switches

use strict;
use warnings;
use Math::Combinatorics;

my $local = '';
my $local_dir = '_';
if ($ARGV[0] && $ARGV[0] eq '--local') {
	$local = '--local';
	$local_dir = '_local_';
}
my @switches = ('-2 sample_R2.fastq.gz', '--dbg trinity', '--no_trim', '--no_norm');
my $test = 1;
my $base = "runDrap $local -1 sample_R1.fastq.gz";

&next_cmd;
foreach my $n ( 1 .. @switches ) {
	my $combinat = Math::Combinatorics->new(count => $n, data => [@switches]);
	while (my @combination = $combinat->next_combination) {
		next_cmd(@combination);
	}
}

sub next_cmd {
	my $cmd = join(' ', $base, @_, "-o drap".$local_dir.$test++);
	print "$cmd\n";
	system($cmd);
}
