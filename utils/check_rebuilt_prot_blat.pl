#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use List::Util qw(min max sum);
if (-t STDIN) {
	print <<HELP;
sed 1d blat.best.tsv | awk -v cov=90 -v id=90 '\$7>=cov&&\$8>=id' | check_rebuilt_prot_blat.pl <-nostat> <-chim> <-fusion> <-prot> <-hit> <-gene> <-f ensembl.pep.all.fa>
 -nostat: don't print statistics line
 -chim: print putative chimeric contigs
 -fusion: print putative fusion contigs
 -prot: print proteins matching putative chimeric or fusion contigs
 -hit: print putative chimeric or fusion contigs hits
 -gene: print genes matching contigs
HELP
	exit;
}

my ($nostat, $chim, $prot, $hit, $fusion, $gene, $ensembl);
GetOptions('nostat' => \$nostat, 'chim' => \$chim, 'prot' => \$prot, 'hit' => \$hit, 'fusion' => \$fusion, 'gene' => \$gene, 'f=s' => \$ensembl);

my %pos = (
	'query'   => 0,
	'target'  => 9,
	'start'   => 10,
	'stop'    => 11,
	'strand'  => 12,
	'tlenght' => 13
);

my (%ensembl, %e, @t1, @t2);
open(IN, "$ensembl") or die "Unable to open reference file $ensembl\n";
while (<IN>) {
	if (/^>/) {
		s/^>//;
		my @t1 = split(/ /,$_);
		my @t2 = split(/:/,$t1[2]);
		($ensembl{$t1[0]}->{chr}, $ensembl{$t1[0]}->{start}, $ensembl{$t1[0]}->{stop}, $ensembl{$t1[0]}->{strand}) = splice(@t2,2,4);
		if ($t1[3] =~ s/gene://) {
			$ensembl{$t1[0]}->{gene} = $t1[3];
			$e{$t1[3]} ='';
		}
	}
}
my $totalProts = keys %ensembl;
my $totalGenes = keys %e;

my (@tmp, %query, %target, %q, %chimeras, %fusions, %tlength, %gene);
while (<STDIN>) {
	chomp;
	@tmp = split("\t", $_);
	$q{$tmp[$pos{query}]} = $tmp[$pos{target}];
	$query{$tmp[$pos{query}]}++;
	die ("Unknow gene for protein $tmp[$pos{query}]\n") unless ($ensembl{$tmp[$pos{query}]}->{gene});
	$gene{$ensembl{$tmp[$pos{query}]}->{gene}}++;
	my ($start, $stop) = ($tmp[$pos{start}], $tmp[$pos{stop}]);
	my $hitLength = $stop - $start + 1;
	$tlength{$tmp[$pos{target}]} = $tmp[$pos{tlenght}];
	push(@{$target{$tmp[$pos{target}]}}, \%{{'query'=>$tmp[$pos{query}], 'start'=>$start, 'stop'=>$stop, 'hitLength'=>$hitLength}});
}

foreach my $t (keys %target) {
	next if (scalar(@{$target{$t}}) < 2);
	# draw hit footprints on target sequence;
	my @seq = (0) x $tlength{$t};
	foreach my $hit (@{$target{$t}}) {
		map { $seq[$_] = 1 } $hit->{start}-1..$hit->{stop}-1;
	}
	my $seq = join('',@seq);
	# identify distinct matching areas
	my @matchingAreaEnds;
	while ($seq =~ m/1+/g) { 
		push(@matchingAreaEnds, pos($seq));
	}
	# assign a matching area to each hit
	foreach my $hit (@{$target{$t}}) {
		my $i;
		for ($i = 0 ; $i <= $#matchingAreaEnds ; $i++) {
			last if $matchingAreaEnds[$i] >= $hit->{stop};
		}
		$hit->{area} = $i+1;
	}
	my @footprints = split(/1+/,$seq);
	if ($#footprints > 1) {
		my $isChim = 0;
		my $isFusion = 1;
		for (my $i = 0 ; $i < scalar(@{$target{$t}}) ; $i++) {
			for (my $j = $i+1 ; $j < scalar(@{$target{$t}}) ; $j++) {
				# proteins from two genes matching to two distinct areas -> chimera
				if ($ensembl{$target{$t}->[$i]->{query}}->{gene} ne $ensembl{$target{$t}->[$j]->{query}}->{gene} && $target{$t}->[$i]->{area} ne $target{$t}->[$j]->{area}) {
					my $twoChr = $ensembl{$target{$t}->[$i]->{query}}->{chr} ne $ensembl{$target{$t}->[$j]->{query}}->{chr};
					my $distance = $ensembl{$target{$t}->[$i]->{query}}->{start} < $ensembl{$target{$t}->[$j]->{query}}->{start} ? $ensembl{$target{$t}->[$j]->{query}}->{start}-$ensembl{$target{$t}->[$i]->{query}}->{stop} : $ensembl{$target{$t}->[$i]->{query}}->{start}-$ensembl{$target{$t}->[$j]->{query}}->{stop};
					# proteins from genes from different chr or genes distance > 10kb matching juxtaposed areas -> not fusion
					$isFusion = 0 if ($twoChr || ($distance > 10000 && abs($target{$t}->[$i]->{area}-$target{$t}->[$j]->{area}) == 1));
					$isChim = 1;
					last if ($isChim && not $isFusion);
				}
			}
			last if ($isChim && not $isFusion);
		}	
		$chimeras{$t} = 1 if ($isChim && not $isFusion);
		$fusions{$t} = 1 if ($isFusion && $isChim);
	}
}

my @protInChimeras;
foreach my $target (keys %chimeras) {
	print "$target\n" if ($chim && not ($prot || $hit));
	if ($chim && $hit) {
		foreach my $h (sort { $a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop} } @{$target{$target}}) {;
			my $queryLocation = join(':',$ensembl{$h->{query}}->{chr}, $ensembl{$h->{query}}->{start}, $ensembl{$h->{query}}->{stop}, $ensembl{$h->{query}}->{strand});
			print join("\t", $h->{query}, $ensembl{$h->{query}}->{gene}, $queryLocation, $target, $h->{start}, $h->{stop}, $tlength{$target})."\n";
		}
	}
	map { push_if_new($_->{query}, \@protInChimeras) } @{$target{$target}};
}

my @protInFusions;
foreach my $target (keys %fusions) {
	print "$target\n" if ($fusion && not ($prot || $hit));
	if ($fusion && $hit) {
		foreach my $h (sort { $a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop} } @{$target{$target}}) {;
			my $queryLocation = join(':',$ensembl{$h->{query}}->{chr}, $ensembl{$h->{query}}->{start}, $ensembl{$h->{query}}->{stop}, $ensembl{$h->{query}}->{strand});
			print join("\t", $h->{query}, $ensembl{$h->{query}}->{gene}, $queryLocation, $target, $h->{start}, $h->{stop}, $tlength{$target})."\n";
		}
	}
	map { push_if_new($_->{query}, \@protInFusions) } @{$target{$target}};
}

if ($prot && $chim) {
	print join("\n", @protInChimeras)."\n";
}

if ($prot && $fusion) {
	print join("\n", @protInFusions)."\n";
}

if ($hit && not ($chim || $fusion)) {
	foreach my $target (keys %target) {
		next if (exists $chimeras{$target} || exists $fusions{$target});
		foreach my $h (sort { $a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop} } @{$target{$target}}) {;
			my $queryLocation = join(':',$ensembl{$h->{query}}->{chr}, $ensembl{$h->{query}}->{start}, $ensembl{$h->{query}}->{stop}, $ensembl{$h->{query}}->{strand});
			print join("\t", $h->{query}, $ensembl{$h->{query}}->{gene}, $queryLocation, $target, $h->{start}, $h->{stop}, $tlength{$target})."\n";
		}
	}
}

if ($gene) {
	print join("\n", keys %gene)."\n";
}

my @qstats = get_stats(\@{[values %query]});
my @tvalues = map { scalar(@{$target{$_}}) } keys %target;
my @tstats = get_stats(\@tvalues);
my $prot_with_max_hit = shift(@{[sort { $query{$b} <=> $query{$a} } keys %query]});
my $contig_with_max_hit = shift(@{[sort { scalar(@{$target{$b}}) <=> scalar(@{$target{$a}}) } keys %target]});
my $numProts = scalar(keys %query);
my $numGenes = scalar(keys %gene);
unless ($nostat) { 
	#printf ("Proteins: %d ; Genes: %d ; Contigs: %d ; Chimeras/fusions/all: %d/%d/%d ; Prot in chimeras/fusions/all: %d/%d/%d; Max hits/prot: %d ; Max hits/contig: %d\n", scalar(keys %query), scalar(keys %gene), scalar(keys %target), scalar(keys %chimeras), scalar(keys %fusions), scalar(keys %fusions)+scalar(keys %chimeras), scalar(@protInChimeras), scalar(@protInFusions), scalar(@protInFusions)+scalar(@protInChimeras), $query{$prot_with_max_hit}, scalar(@{$target{$contig_with_max_hit}}));
	printf ("Proteins: %d (%.0f%%); Genes: %d (%.0f%%); Contigs: %d ; Chimeras/fusions/all: %d/%d/%d ; Prot in chimeras/fusions/all: %d/%d/%d; Hits by prot min/median/mean/max/moreThanOne: %s ; Hits by contig: %s\n", $numProts, $numProts / $totalProts * 100 , $numGenes, $numGenes / $totalGenes * 100, scalar(keys %target), scalar(keys %chimeras), scalar(keys %fusions), scalar(keys %fusions)+scalar(keys %chimeras), scalar(@protInChimeras), scalar(@protInFusions), scalar(@protInFusions)+scalar(@protInChimeras), join('/', @qstats), join('/', @tstats));
}

sub push_if_new {
	my ($first, $second) = @_;
	foreach my $item (@$second) {
		return 0 if ($item eq $first);
	}
	push(@$second, $first);
	return 1;
}

sub get_stats {
	my $ref = shift;
	# return min, median, mean, max, moreThanOne
	my $median;
	my $mid = int @$ref/2;
	my @sorted_values = sort { $a <=> $b } @$ref;
	if (@$ref % 2) {
	    $median = $sorted_values[$mid];
	} else {
	    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
	}
	my $mean = sprintf("%.2f",sum(@$ref)/scalar(@$ref));
	my $min = min @$ref;
	my $max = max @$ref;
	my $moreThanOne = 0;
	map { $moreThanOne++ if $_ > 1 } @$ref;
	return ($min, $median, $mean, $max, $moreThanOne);
}
