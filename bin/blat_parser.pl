#!/usr/bin/perl

use strict;
use Pod::Usage;
use Getopt::Long;

my ($man, $help, $focus, $nb_best);
GetOptions(
	'f=s'  => \$focus,
	'b=i'  => \$nb_best,
	'help' => \$help,
	'man'  => \$man,
);

pod2usage(1) if ($help);
pod2usage(-exitstatus => 0, -verbose => 2) if ($man);
pod2usage(-exitstatus => 0, -verbose => 1) if (-t STDIN);
$focus ||= 'q';
$nb_best ||= 1;

my (%option, $line);

$option{'m'} = 1;
$option{'p'} = pslIsProtein($line) || 0;

my %h;
my $f = $focus eq 'q' ? 9 : 13;
while ($line = <STDIN>) {
	next unless ($line =~ /^\d+/);
	my @t = split(/\t/,$line);
	my $pid = sprintf("%.2f",get_pid(@t));
	$t[11] += 1;
	$t[15] += 1;
	my $Qcov = sprintf("%.2f",(abs($t[12]-$t[11])-$t[5]+1)/$t[10]*100);
	my $Tcov = sprintf("%.2f",(abs($t[16]-$t[15])-$t[7]+1)/$t[14]*100);
	my @strand = split('',$t[8]);
	my @list = ($t[9],get_score(@t),$t[11],$t[12],$strand[0],$t[10],$Qcov,$pid,'.',$t[13],$t[15],$t[16],$strand[1]||'.',$t[14],$Tcov,$t[4],$t[5],$t[6],$t[7]);
	if (exists($h{$t[$f]})) {
		scalar(@{$h{$t[$f]}}) == $nb_best ? push_if_better(\@{$h{$t[$f]}}, \@list) : push_and_sort(\@{$h{$t[$f]}},\@list);
	}
	else {
		push(@{$h{$t[$f]}},\@list);
	}
}

#print join("\t",qw(qName score qStart qEnd qSize qCoverage identity strand tName tStart tEnd tSize tCoverage qNumInsert qBaseInsert tNumInsert tBaseInsert))."\n";
print join("\t",qw(qName score qStart qEnd qStrand qLength %qCoverage %identity %similarity tName tStart tEnd tStrand tLength %tCoverage qNumInsert qBaseInsert tNumInsert tBaseInsert))."\n";

foreach my $seq (sort keys %h) {
	foreach my $match (@{$h{$seq}}) {
		print join("\t",@$match)."\n";
	}
}

sub push_if_better {
	my ($first, $second) = @_;
	return unless (${$$first[-1]}[1] < $$second[1]);
	$$first[-1] = $second;
	my @t = sort { $$b[1] <=> $$a[1] } @$first;
	@$first = @t;
	return;
}

sub push_and_sort {
	my ($first, $second) = @_;
	push(@$first, $second);
	my @t = sort { $$b[1] <=> $$a[1] } @$first;
	@$first = @t;
	return;
}

sub get_pid {
	my @line = @_;
	my $pid = (100.0 - (&pslCalcMilliBad(@line) * 0.1));
	return $pid;
}

sub get_score {
	my @line = @_;
	my $sizeMul;
	if ($option{'p'}) {
			$sizeMul = 3;
	} else {
			$sizeMul = 1;
	}
	my $score = $sizeMul * ($line[0] + ($line[2] >> 1)) - $sizeMul * $line[1] - $line[4] - $line[6];
	return $score;
}

sub pslCalcMilliBad {
	my @cols = @_;

	# sizeNul depens of dna/Prot
	my $sizeMul;
	if ($option{'p'}) {
			$sizeMul = 3;
	} else {
			$sizeMul = 1;
	}

	# cols[0]  matches
	# cols[1]  misMatches
	# cols[2]  repMaches
	# cols[4]  qNumInsert
	# cols[6]  tNumInsert
	# cols[11] qStart
	# cols[12] qEnd
	# cols[15] tStart
	# cols[16] tEnd

	my $qAliSize = $sizeMul * ($cols[12] - $cols[11]);
	my $tAliSize = $cols[16] - $cols[15];

	# I want the minimum of qAliSize and tAliSize
	my $aliSize;
	$qAliSize < $tAliSize ? $aliSize = $qAliSize : $aliSize = $tAliSize;

	# return 0 is AliSize == 0
	return 0 if ($aliSize <= 0);

	# size diff
	my $sizeDiff = $qAliSize - $tAliSize;
	if ($sizeDiff < 0) {
			if ($option{'m'}) {
					$sizeDiff = 0;
			} else {
					$sizeDiff = -($sizeDiff);
			}
	}

	# insert Factor
	my $insertFactor = $cols[4];
	$insertFactor += $cols[6] unless ($option{'m'});
	my $milliBad = (1000 * ($cols[1]*$sizeMul + $insertFactor + &round(3*log( 1 + $sizeDiff)))) / ($sizeMul * ($cols[0] + $cols[2] + $cols[1]));
	return $milliBad;
}

sub round {
	my $number = shift;
	return int($number + .5);
}

sub pslIsProtein {
	my $l = shift;
	my @cols = split(/\t/,$l);
	
	# cols[8]  strand
	# cols[14] tSize
	# cols[15] tStart
	# cols[16] tEnd
	# cols[17] blockCount
	# cols[18] blockSizes
	# cols[20] tStarts
	
	my $lastBlock = $cols[17] - 1;
	my $strand = substr($cols[8],1,1) || '';
	my @blockSizes = split(/,/, $cols[18]);
	my @tStarts = split(/,/, $cols[20]);
	
	return (
		($strand eq '+' && $cols[16] == $tStarts[$lastBlock]+3*$blockSizes[$lastBlock])
		|| 
		($strand eq '-' && $cols[15] == $cols[14]-($tStarts[$lastBlock]+3*$blockSizes[$lastBlock]))
	);
}

=head1 NAME

 blat_parser.pl

=head1 SYNOPSIS

 Usage: cat blat.psl | blat_parser.pl [options]

=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
 Print a brief help message and exits.
    
=item B<-man>
    
 Prints the manual page and exits.
    
=item B<-f>
    
 Focus best matches extraction on query (q) or target (t) [q]

=item B<-b>
    
 Max number of best matches to extract [1]

=back

=head1 DESCRIPTION

 Read as STDIN blat results in psl format.
 Compute score and identity percent as done on UCSC Web-based Blat (http://genome.ucsc.edu/FAQ/FAQblat.html#blat4).
 Compute qCoverage and tCoverage which are respectively the query and target fractions covered by the match.
 Print on STDOUT the max_nb_matches best matches (matches are ordered by score).
 Please refer to psl format page (http://genome.ucsc.edu/FAQ/FAQformat.html#format2) to get definitions of output fields.

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=head1 VERSION

 2

=head1 DATE

 2013

=head1 KEYWORDS

 blat psl parser

=cut

