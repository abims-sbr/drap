#!/usr/bin/perl

use strict;
use Pod::Usage;
use Getopt::Long;

my ($fa, $bed, $log, $keep, $overlap, $orf_length_file, $help, $man);

GetOptions(
	'f=s' => \$fa,
	'b=s' => \$bed,
	'l|log=s' => \$log,
	'keep' => \$keep,
	'overlap' => \$overlap,
	'orf-length=s' => \$orf_length_file,
	'help' => \$help,
	'man' => \$man
);

pod2usage(-exitstatus => 0, -verbose => 99, -sections => "NAME|SYNOPSIS|OPTIONS|DESCRIPTION") if ($help);
pod2usage(-verbose => 2) if ($man);
(-e $fa || $fa =~ /stdin|^-/) || pod2usage("-f is not a file: $fa");
(-e $bed) || pod2usage("-b is not a file: $bed");
qx(which bedtools) || die("Unable to find bedtools\n");
# log file
my $fhl= *STDERR;
if($log) { open(LOG, ">$log") or die "Can't open file $log"; $fhl= *LOG; }

my (%orf, @t, $tmp_start, $prev_ctg, $prev_end);
open(BED, "sed 1d $bed ".q(| cut -f1,3,7,8 | awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$2}' | sort -k1,1 -k2,2g | bedtools merge -c 4 -o distinct -i - |)) or die "Can't process file $bed\n";
# open gives following fields: contig name, ORF start, ORF end, contig length
while (<BED>) {
	chomp(@t = split(/\t/, $_));
	# ORFs start position extracted from the bed file are 0-based; conversion to 1-based
	$t[1]++;
	if (exists $orf{$t[0]}) {
		if ($overlap) {
			# cut with overlap between previous end +1 and current start -1
			$tmp_start = $orf{$t[0]}->[-1]->{'end'}+1;
			$orf{$t[0]}->[-1]->{'end'} = $t[1]-1; # extend previous slice in 3' until current start -1
			push (@{$orf{$t[0]}}, \%{{ 'start' => $tmp_start, 'end' => $t[2] }}); # extend current slice in 5' until previous end +1
		}
		else {
			# cut in the middle of two consecutives ORFs
			$tmp_start = int(($t[1]-$orf{$t[0]}->[-1]->{'end'})/2)+$orf{$t[0]}->[-1]->{'end'};
			$orf{$t[0]}->[-1]->{'end'} = $tmp_start-1; 
			push (@{$orf{$t[0]}}, \%{{ 'start' => $tmp_start, 'end' => $t[2] }}); 
		}
		$orf{$t[0]}->[-1]->{'orf_size'} = $t[2]-$t[1]+1;
	}
	else {
		push (@{$orf{$t[0]}}, \%{{ 'start' => 1, 'end' => $t[2], 'orf_size' => $t[2]-$t[1]+1 }}); # push first slice
		$orf{$prev_ctg}->[-1]->{'end'} = $prev_end if (exists $orf{$prev_ctg}); # extend last slice of previous contig to contig end
	}
	$prev_ctg = $t[0];
	$prev_end = $t[3];
}
$orf{$prev_ctg}->[-1]->{'end'} = $prev_end; # extend last slice of last contig to contig end
close BED;

my ($fh, $fh_length);
open($fh_length, ">$orf_length_file") or die "Can't open file $orf_length_file" if($orf_length_file); 
if($fa =~ /stdin|^-/) { $fh = *STDIN }
else { open(F, $fa) or die "Can't open $fa\n"; $fh = *F }
local $/ = '>';
my ($name, $seq, $c, $l);
print $fhl join("\t", qw(Contig #ORFs slice_coords))."\n";
while (<$fh>) {
	chomp;
	($name, $seq) = $_ =~ /(\S+).*?\n(.+)/s or next;
	$seq =~ s/\n//g;
	$c = 0;
	if (exists $orf{$name}) {
		foreach my $slice (@{$orf{$name}}) {
			print_fasta($name.'_'.++$c, $slice->{'orf_size'}, substr($seq, $slice->{'start'}-1, $slice->{'end'}-$slice->{'start'}+1), $fh_length);
		}
		printf $fhl ("%s\t%d\t%s\n", $name, $c, join(',', map { $_->{'start'}.'-'.$_->{'end'} } @{$orf{$name}}));
	}
	elsif ($keep) {
		print_fasta($name.'_1', 0, $seq, $fh_length);
		printf $fhl ("%s\t%d\t%s\n", $name, $c, '1-'.length($seq));
	}
}
close $fh;
close $fhl;

sub print_fasta {
	my ($name, $s, $seq, $fh) = @_;
	if ($fh) {
		print ">$name\n";
		print $fh "$name\t$s\n";
	} else {
		print ">$name [ORF_SIZE=$s]\n";
	}
	while ($seq) {
		print substr($seq, 0, 60, '')."\n";
	}
}

=head1 NAME

transdecoder_bed_to_false_fusion_cutter.pl

=head1 SYNOPSIS

transdecoder_bed_to_false_fusion_cutter.pl [-h|options] -f contigs.fa|stdin -b file.transdecoder.bed

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<--keep>

Keep contigs without ORF predicted by transdecoder (only cut false fusion).

=item B<--overlap>

Cut false fusion contigs producing overlapping ends.

=item B<--orf-length>

ORF length output file. If not provided, ORF length will be printed to fasta header line.

=item B<-f>

Input fasta file.

=item B<-b>

Input transdecoder bed file.

=item B<-l, --log>

Log file. If log file not provided, print log messages on STDERR.

=back

=head1 DESCRIPTION
  
Cut false fusion contigs defined as contigs with multiple non overlapping ORFs predicted by transdecoder.
Overlapping ORFs predicted by transdecoder are merged with bedtools merge before processing.
In ouput fasta format, contig names are concatenated with "_<slice_number>".
Two strategies to cut the false fusion contigs are available.

The default behavior cuts in the middle of two consecutives ORFs
producing front-end sequences:

If x stretches define ORF positions:
-------xxxxxxxxxx----------xxxxxxxxxxxxxxxxxxxxxxx--------

Default behavior will give:
-------xxxxxxxxxx-----
                      -----xxxxxxxxxxxxxxxxxxxxxxx--------

The overlap strategy produces overlapping ends.
First base is the first base next to the previous ORF.
Terminal base is the last base before the consecutive ORF.

Overlap strategy will give:
-------xxxxxxxxxx----------
                 ----------xxxxxxxxxxxxxxxxxxxxxxx--------

By default, contigs without ORFs predicted by transdecoder are discarded.

Dependancy: bedtools.

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=head1 VERSION

1

=head1 DATE

2014

=head1 KEYWORDS

ORF fasta chimera cut

=cut

