#!/usr/bin/perl -w

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO::Writer::HSPTableWriter;
use Bio::Search::Tiling::MapTiling;
use Bio::Search::Hit::GenericHit; 

my $self_id_cutoff = 96;
my $self_cov_cutoff = 60;
my $self_global_cutoff = 80;
my $debug = 0;
my $log;
my $man;

GetOptions(
	'i=i' => \$self_id_cutoff,
	'c=i' => \$self_cov_cutoff,
	'g=i' => \$self_global_cutoff,
	'l=s' => \$log,
	'man' => \$man,
);

pod2usage(-exitstatus => 0, -verbose => 2) if ($man);
pod2usage(-exitstatus => 0, -verbose => 1) if (-t);

# log file
my $logfh= *STDERR;
if($log) { open(LOG, ">$log}") or die "Can't open file $log"; $logfh= *LOG; }

my $seqI = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta');
my $seqO = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'fasta') or die($!);
my $writer = Bio::SearchIO::Writer::HSPTableWriter->new(-columns => [qw(query_name query_length score start_query end_query strand_query start_hit end_hit strand_hit)]);
my $all;
my $others;
my $type;
print $logfh ("# Contigs considered as putative chimera (seqID length chimType newFragment):\n");
while (my $seq = $seqI->next_seq) {
	my $seq_id = $seq->display_id;
	(my $name = $seq_id) =~ s|/.+||;
	cut_N_strech($seq);
	my $length = $seq->length;
	my $split = $length;
	my $factory = Bio::Tools::Run::StandAloneBlast->new(
		-program => 'blastn',
		-S => '3', # align on both strand
		-F => 'F', # don't filter query seq
		-m => 'T'  # use megablast
	);
	my $report = $factory->bl2seq($seq, $seq);
	my $result = $report->next_result;
	my $hit = $result->next_hit if ($result);
	if ($hit) {
		my @hsps = $hit->hsps;
		shift(@hsps); # remove first self complete hit
		if (@hsps && $hsps[0]->percent_identity >= $self_id_cutoff) {
			my @new_hsps = @hsps;
			my $new_hit = Bio::Search::Hit::GenericHit->new(
				-name => $hit->name,
				-algorithm => $hit->algorithm,
				-hsps => \@new_hsps,
				-hsp_factory  => $hit->hsp_factory
			); # build a new hit object without the first self complete hit	to compute tiled hsps global length
			my $tiling = Bio::Search::Tiling::MapTiling->new($new_hit);
			if ($debug) {
				my @tmp = split(/\n/, $writer->to_string($result,[qw(query_name score start_query end_query strand_query start_hit end_hit strand_hit)]));
				splice(@tmp, 0, 3);
				push(@tmp, sprintf ("-> Duplication rate: %.0f%%\n\n", $tiling->length('query')/$length*100));
				$all .= join("\n", @tmp);
			}
			
			# We process three cases:
			#
			# A one block trans self match described in only one hsp (processed if hsp length > contig length/100*$self_cov_cutoff)
			# k31_Locus_117_Transcript_2/2_Confidence_1.000_Length_2677   k31_Locus_117_Transcript_2/2_Confidence_1.000_Length_2677   99.36   2677   17  0   1     2677  2677  1     0.0  5172
			#
			# A two blocks trans self match described in two hsps (processed if hsp length*2 > contig length/100*$self_cov_cutoff)
			# k25_Locus_55_Transcript_2/2_Confidence_1.000_Length_5939    k25_Locus_55_Transcript_2/2_Confidence_1.000_Length_5939    100.00  2953   0   0   2987  5939  2953  1     0.0  5854
			# k25_Locus_55_Transcript_2/2_Confidence_1.000_Length_5939    k25_Locus_55_Transcript_2/2_Confidence_1.000_Length_5939    100.00  2953   0   0   1     2953  5939  2987  0.0  5854
			#
			# A two blocks cis self match described in two hsps (processed if hsp length*2 > contig length/100*$self_cov_cutoff)
			# k31_Locus_5_Transcript_2/2_Confidence_1.000_Length_6005     k31_Locus_5_Transcript_2/2_Confidence_1.000_Length_6005     100.00  2640   0   0   3366  6005  1     2640  0.0  5233
			# k31_Locus_5_Transcript_2/2_Confidence_1.000_Length_6005     k31_Locus_5_Transcript_2/2_Confidence_1.000_Length_6005     100.00  2640   0   0   1     2640  3366  6005  0.0  5233
			#
			if (
				$hsps[0]->start('query') == $hsps[0]->start('hit') && 
				$hsps[0]->end('query') == $hsps[0]->end('hit') && 
				$hsps[0]->strand('query') != $hsps[0]->strand('hit')
			) {
				if ($hsps[0]->length('query') >= $length/100*$self_cov_cutoff) {
					$split = $hsps[0]->length('query')/2+$hsps[0]->start('query');
					$type = 'trans_self_one_block';
				}
				goto SPLIT;
			}
			if ($hsps[0]->strand('query') != $hsps[0]->strand('hit') && non_overlapping($hsps[0])) {
				if ($hsps[0]->length('query')*2 >= $length/100*$self_cov_cutoff) {
					$split = $hsps[0]->start('query') < $hsps[0]->start('hit') ? $hsps[0]->start('hit') : $hsps[0]->start('query');
					$type = 'trans_self_two_blocks';
				}
				goto SPLIT;
			}
			if ($hsps[0]->strand('query') == $hsps[0]->strand('hit') && non_overlapping($hsps[0])) {
				if ($hsps[0]->length('query')*2 >= $length/100*$self_cov_cutoff) {
					$split = $hsps[0]->start('query') < $hsps[0]->start('hit') ? $hsps[0]->start('hit') : $hsps[0]->start('query');
					$type = 'cis_self';
				}
			}
			SPLIT:
			if ($split < $length) {
				my $newFrag;
				if ($split >= $length/2) {
					$split = int(--$split); # split include the first base of the duplicated block
					$seq->seq($seq->subseq(1, $split));
					$newFrag = "1-$split";
				}
				else {
					$split = int($split);
					$seq->seq($seq->subseq($split, $length));
					$newFrag = "$split-$length";
				}
				print $logfh join("\t",$seq_id, $length, $type, $newFrag)."\n";
			}
			elsif ($tiling->length('query') >= $length/100*$self_global_cutoff) {
				# Contigs with repeated blocks are discarded if #bases repeated greater than self_global_cutoff
				printf $logfh ("%s\t%d\trepeated_blocks_(%.0f%%)\tdiscarded\n", $seq_id, $length, $tiling->length('query')/$length*100);
				next;
			}
			else {
				if ($hsps[0]->length('query')*2 >= $length/100*$self_cov_cutoff) {
					my @tmp = split(/\n/, $writer->to_string($result,[qw(query_name score query_length start_query end_query strand_query start_hit end_hit strand_hit)]));
					splice(@tmp, 2, 1);
					push(@tmp, sprintf ("-> Duplication rate: %.0f%%\n\n", $tiling->length('query')/$length*100));
					$others .= join("\n", @tmp);
				}
			}
		}
	}
	$seqO->write_seq($seq);
}
if ($others) {
	print $logfh ("\n# Overlapping hsps under global cutoff ($self_global_cutoff%):\n");
	print $logfh ($others);
}

if ($all) {
	print $logfh ("\n# All processed matches:\n");
	print $logfh ($all);
}

sub cut_N_strech {
	my $seq = shift;
	my $gaplen = 14;
	my $gap = ('N') x $gaplen;
	for (my $gapb = index($seq->seq, $gap); $gapb >= 0; ) {
		my $reflen = $seq->length;
		my $gape = $gapb + $gaplen;
		$gape++ while ($gape < $reflen && substr($seq->seq, $gape, 1) eq "N"); 
		my $gapw = $gape - $gapb;
		my $keep = 6 + ($gapw % 3);
		my $facut = substr($seq->seq, 0, $gapb).substr("NNNNNNNNN", 0, $keep).substr($seq->seq, $gape);
		$seq->seq($facut);
		$gapb = index($seq->seq, $gap);
	}
}

sub non_overlapping {
	my $hsp = shift;
	return 0 if ($hsp->start('query') < $hsp->start('hit') && $hsp->end('query') > $hsp->start('hit'));
	return 0 if ($hsp->start('query') > $hsp->start('hit') && $hsp->end('hit') > $hsp->start('query'));
	return 1;
}
	
=head1 NAME

self_chimeras_filter.pl

=head1 SYNOPSIS

cat transcripts.fa | self_chimeras_filter.pl [options]

=head1 OPTIONS

=over 8

=item B<-man>

Print the man page and exit.

=item B<-i>

identity cutoff: only matches with identity greater or equal than -i will be processed [96]

=item B<-c>

coverage cutoff: the longest self match have to cover at least -c percent of the contig length to consider contig as a chimera [60]

=item B<-g>

global cutoff: all self matches have to cover at least -g percent of the contig length to consider contig as a chimera [80]

=item B<-log>

Print log messages to this file.

=back

=head1 DESCRIPTION
  
  Read a fasta file as STDIN.
  For each contig, compress gaps greater than 14 n's in a row (avoid bl2seq crash).
  Perform a bl2seq alignment for each contig against itself.
  Considering only self matches greater or equal than identity cutoff, a contig is considered as putative chimera if:
    - the longest (i.e. the first) self match covers at least -c percent of the contig length
    - or all self matches length cover at least -g percent of the contig length
  The position to split a putative chimera depends on the self match type: 
    - if the chimera is a one block match, position is the middle of the match
    - if the chimera is a two blocks match, position is the start of the second block
  Contigs with repeated blocks are discarded.
  Write all contigs free of chimeras to STDOUT. Write putative chimeras processing log to STDERR, if log file is not
  provided.

  One block trans self match example:
    # % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    99.36   2677   17  0   1     2677  2677  1     0.0    5172

  Two blocks trans self match example:
    # % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    100.00  2953   0   0   1     2953  5939  2987  0.0    5854
    100.00  2953   0   0   2987  5939  2953  1     0.0    5854
  
  Two blocks cis self match example:
    # % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    100.00  2640   0   0   1     2640  3366  6005  0.0    5233
    100.00  2640   0   0   3366  6005  1     2640  0.0    5233

  Repeated blocks example:
    # % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    97.06   34     1   0   71    104   34    1     9e-15  60.0
    97.06   34     1   0   70    103   34    1     9e-15  60.0
    97.06   34     1   0   69    102   34    1     9e-15  60.0
    97.06   34     1   0   1     34    102   69    9e-15  60.0
    97.06   34     1   0   1     34    103   70    9e-15  60.0
    97.06   34     1   0   1     34    104   71    9e-15  60.0
    100.00  27     0   0   78    104   35    9     5e-13  54.0
    100.00  27     0   0   9     35    104   78    5e-13  54.0

=head1 AUTHORS
    
 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3
    
=cut
