#!/usr/bin/perl
=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut
use strict;
use FindBin;
use lib ("$FindBin::RealBin");
use Carp 'croak';
use File::Basename;
use File::Copy;
use File::Compare;
use File::Find;
use File::Path qw(make_path remove_tree);
use AsmUtils;

my $dir  = shift;
my $step = shift;
my $opt = get_drap_config($dir) or croak("Unable to read drap config file .drap_conf.json inside $dir");

our $no_qacct = $opt->{no_qacct};
my $username = $opt->{local} ? 'anonymous' : getpwuid($<);
my $msg = {
	'preprocess' => {
		1 => "Something wrong happened during trimming step: see SRC file",
		2 => "Something wrong happened during filtering step: see SRC file",
		3 => "Something wrong happened during normalization step: see SRC file",
		4 => "Can't find strings 'sequences found' or 'Done' inside SRC file",
		5 => "Found string ERROR inside file SRC",
	},
	'dbg' => {
		1 => "Can't find string 'Exporting transcripts to...' inside file SRC",
		2 => "Can't find string 'Finished extracting transcripts' inside file SRC",
		3 => "Can't find string 'Butterfly assemblies are written to...' inside SRC file",
	},
	'merge' => {
		1 => "Can't find string 'Clean locus in...' for each kmer inside SRC file",
		2 => "Can't find string '...seqclean, without a detectable error' for each kmer inside SRC file",
		3 => "Number of contigs after chimFilter not equal to number of contigs before chimFilter minus discarded ones in SRC directory",
	},
	'post_asm' => {
		1 => "Can't find string 'transdecoder is finished.' inside SRC file",
		2 => "Cutting report all_contigs.raw.fa.transdecoder_cutter.log not compliant with number of contigs inside all_contigs.raw.fa.transdecoder_cutter.fa in SRC directory",
		3 => "Number of contigs after vecFilter not equal to number of contigs before vecFilter minus discarded ones in SRC directory",
	},
	'rmbt' => {
		1 => "Can't find bam or flagstat file in SRC directory",
		2 => "Can't find string '[main] Real time:' inside SRC file",
		3 => "Number of sequences processed in files SRC are not equal",
		4 => "Number of sequences processed in files SRC are not concordant (no. sequences on first line of flagstat should be twice as much as no. sequences in sampe file)",
		5 => "Can't find *.Log.final.out file in SRC directory",
		6 => "Number of sequences found in files SRC are not concordant (no. processed sequences in samse or Log.final.out file should be equal to no. sequences in fastq file)",
		7 => "Can't find string '[main] Real time:' twice inside SRC file",
		8 => "Number of sequences processed in samse mode not equal to number of sequences processed in aln mode in SRC file",
		9 => "Number of sequences processed in files SRC are not concordant (no. sequences on first line of flagstat should be equal to no. sequences in samse file)",
	},
	'rmbt_editing' => {
		1 => "The number of valid BAM files produced at the SRC pass of the rmbt-editing step is not consistent with input data",
		2 => "Number of contigs after samCorrectVariation not equal to number of contigs before samCorrectVariation (inside SRC files)",
	},
	'rmbt_filtering' => {
		1 => "The number of valid BAM files produced at the rmbt-filtering step is not consistent with input data",
		2 => "Number of BAM files parsed inside SRC file not equal to number of BAM files produced at rmbt-filtering step",
		3 => "Can't find string 'COMPLETED: Processed .* mapped fragments...' inside SRC file",
	},
	'reference' => {
		1 => "Temporary directory SRC exists. All tasks did not ended successfully. Check err and log files inside reference directory",
	},
	'meta_rmbt' => {
		1 => "The number of valid BAM files produced at the meta_rmbt step is not consistent with input data",
	},
	'meta_filter' => {
		1 => "Number of BAM files parsed inside SRC file not equal to number of BAM files produced at meta_rmbt step",
		2 => "Can't find string 'COMPLETED: Processed .* mapped fragments...' inside SRC file",
		3 => "Can't find string 'transdecoder is finished.' inside SRC file",
	},
	'cd-hit'    => "Can't find string 'program completed !' inside file SRC",
	'transrate' => "Can't find string 'Writing analysis results to assemblies.csv' inside file SRC",
	'open'      => "Can't open file SRC",
	'error'     => "Error file SRC is not empty",
	'empty'     => "File SRC does not exist or is empty",
	'qacct'     => "Failed or exit status of command [qacct -o $username -j SRC] not equal to 0"
};
$msg->{meta_reference} = $msg->{reference};

# Check only one step and exit
if ($step) {
	$opt->{step} = $step;
	$opt->{check_postprocess} = 1;
	unlink(glob("$dir/00-ERROR_AT_".uc($step).'_STEP'));
	step_complete($opt) ? exit 0 : err_exit($opt);
}

# Check whole assembly except postprocess step
my $success = exists($opt->{meta}) ? '00-META-ASSEMBLY_COMPLETE' : '00-ASSEMBLY_COMPLETE';
unlink(glob("$dir/$success"), glob("$dir/00-ERROR_AT_*_STEP"));
foreach $step (@{$opt->{steps}}) {
	last if ($step =~ /postprocess/);
	$opt->{step} = $step;
	print "Check $step\n";
	if ($step eq 'dbg' && $opt->{dbg} eq 'oases') {
		$opt->{task} = 1;
		foreach my $kmer (@{$opt->{kmers}}) {
			next if ($opt->{kmer_status}->{$kmer} == 1);
			$opt->{current_kmer} = $kmer;
			err_exit($opt) unless (step_complete($opt));
			$opt->{task}++;
		}
	} else {
		err_exit($opt) unless (step_complete($opt));
	}
}
process_cmd(0, "touch $dir/$success");
my $cov_filter_dir = exists($opt->{meta}) ? $opt->{dir_list}->[3] : $opt->{dir_list}->[5];
unless ($opt->{optimize}) {
	foreach (glob("$cov_filter_dir/transcripts_fpkm_*.fa $cov_filter_dir/coding_transcripts_fpkm_*.fa")) {
		symlink(basename($cov_filter_dir).'/'.basename($_), $dir.'/'.basename($_));
	}
}
print "Checking complete\n";
set_drap_config($opt->{outdir}, $opt);
exit 0;

sub err_exit {
	my $opt = shift;
	my $err_file = "$dir/00-ERROR_AT_".uc($opt->{step}).'_STEP';
	my $err_msg = get_msg($opt);
	$err_msg =~ s|SRC|$opt->{err_src}|;
	open(ERR, ">$err_file") or croak "Can't open file $err_file";
	print ERR "$err_msg\n";
	close ERR;
	set_drap_config($opt->{outdir}, $opt);
	exit 1;
}

sub get_msg {
	my $opt = shift;
	if (exists($msg->{$opt->{err_code}})) {
		return $msg->{$opt->{err_code}};
	} else {
		if ($opt->{err_code} =~ /(\w+)\.(\d+)/) {
			return $msg->{$1}->{$2};
		} else {
			return $msg->{$opt->{step}}->{$opt->{err_code}};
		} 
	}
}
