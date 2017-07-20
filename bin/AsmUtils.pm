#!/usr/bin/perl
=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut
package AsmUtils;

use strict;
no strict 'refs';
use Exporter;
use Carp 'croak';
use Pod::Usage;
use File::Find;
use File::Basename;
use JSON;
use FindBin;
use lib("$FindBin::Bin") ;
use ConfigFile ;
my $VERSION = '1.0';
our @ISA = qw(Exporter);
our @EXPORT = qw(
	step_complete process_cmd get_more_recent_file set_drap_config
	get_drap_config get_nb_frags submit_job get_cfg_param print_msg
	printOutMessage create_report write_shell clean_directories
	clean_directory clean_rmbt_directory clean_rmbt_reference
	valid_R1_R2 gunzip gzip find_in set_extended_path get_scheduler_type
	set_env_variables is_prot_fasta
);

sub step_complete {
	my $opt = shift;
	my $step = shift;

	# run checking sub
	$step ||= $opt->{step};
	my $subref = $step.'_complete';
	my $return  = $subref->($opt);

	# all subs return an array ref containing:
	# - return code 1|0
	# - err code to catch msg
	# - err src where error was found
	$opt->{err_code} = $return->[1];
	$opt->{err_src} = $return->[2];
	return $return->[0];
}

sub process_cmd {
	my $return = shift;
	my $cmd = shift;
	chomp(my $out = qx($cmd));
	return $out if ($return);
	croak "Error, command: $cmd\nDied with error code: $out" if ($out);
}

sub get_more_recent_file {
	my $dir = shift;
	my $pattern = shift;
	my @files = sort { (stat($b))[9] <=> (stat($a))[9] } glob("$dir/$pattern");
	return shift(@files);
}

sub set_drap_config {
	my $dir = shift;
	my $opt = shift;
	open(JSON, ">$dir/.drap_conf.json") or croak "Can't open file $dir/.drap_conf.json";
	print JSON to_json($opt, { utf8 => 1, pretty => 1, canonical => 1 });
	close(JSON);
}

sub get_drap_config {
	my ($dir, $retrieve, $opt) = @_;
	local $/;
	open(JSON, "$dir/.drap_conf.json") or croak "Can't open file $dir/.drap_conf.json";
	my $json = <JSON>;
	close(JSON);
	if ($opt) {
		# $opt exists, it is a restart
		my $from_json = from_json($json);
		if ($retrieve) {
			# no options passed, get previous configuration
			foreach my $key (keys %$opt) {
				next if ($key eq 'write' || $key eq 'run');
				(my $short_option_key = $key) =~ s/.+\|//;
				$short_option_key =~ s/-/_/g;
				$opt->{$key} = $from_json->{$short_option_key};
			}
		}
		# fill $opt with data which are not not options data (i.e rmbt checking...)
		foreach my $key (keys %$from_json) {
			(my $dash_key = $key) =~ s/_/-/g;
			next if (exists $opt->{$dash_key});
			$opt->{$dash_key} = $from_json->{$key};
		}
		return $opt;
	} else {
		return from_json($json);
	}
}
sub get_nb_frags {
	my $opt = shift;
	my $max = shift;
	if ($opt->{restart}){
		$opt->{nb_frags} = get_cfg_param($opt->{outdir}, 'nb_frags');
		$opt->{max_nb_frags} = get_cfg_param($opt->{outdir}, 'max_nb_frags');
	} else {
		printf("Computing the%s number of fragments to process for memory reservation... ", $max && $max > 1 ? ' max' : '');
		open(TSV, ">>$opt->{outdir}/.processed_reads.tsv") or croak "Can't open file $opt->{outdir}/.processed_reads.tsv";
		my $nb_frags;
		for (my $i = 0 ; $i < scalar(@{$opt->{R1}}) ; $i++) {
			$nb_frags = process_cmd(1, "$opt->{binpath}/fastq_nb_reads.sh $opt->{R1}->[$i] $opt->{filter}");
			$opt->{nb_frags} += $nb_frags;
			print TSV $opt->{R1}->[$i]."\t".$nb_frags."\t".$opt->{filter}."\n";
			print TSV $opt->{R2}->[$i]."\t".$nb_frags."\t".$opt->{filter}."\n" if ($opt->{paired});
		}
		close(TSV);
		$opt->{max_nb_frags} = $nb_frags unless (defined $opt->{max_nb_frags} && $opt->{max_nb_frags} > $nb_frags);
		printf("%s%d reads\n",$opt->{paired} ? '2x' : '',$opt->{nb_frags});
	}
	return $max ? $opt->{max_nb_frags} : $opt->{nb_frags};
}

sub get_nb_reads {
	my $opt = shift;
	my $fastq = shift;
	my $nb_reads;
	if (-e "$opt->{outdir}/.processed_reads.tsv") {
		$nb_reads = shift(@{[find_inside_file("$opt->{outdir}/.processed_reads.tsv", $fastq.'\t(\d+)\t'.$opt->{filter})]});
	}
	unless ($nb_reads) {
		$nb_reads = process_cmd(1, "$opt->{binpath}/fastq_nb_reads.sh $fastq $opt->{filter}");
		open(TSV, ">>$opt->{outdir}/.processed_reads.tsv") or croak "Can't open file $opt->{outdir}/.processed_reads.tsv";
		print TSV $fastq."\t".$nb_reads."\t".$opt->{filter}."\n";
		close(TSV);
	}
	return $nb_reads;
}

sub get_cfg_param {
	my $dir = shift;
	my $key = shift;
	return get_drap_config($dir)->{$key}||'';
}

sub submit_job {
	my ($cmd, $logfh) = @_;
	my $submit_msg = process_cmd(1, $cmd);
	print "$submit_msg\n";
	$submit_msg =~ /Your job(?:-array)? (\d+)(\..+)? \(".*"\) has been submitted/;
	return $1 if ($1);
	croak "Unable to get job_id from qsub return: $submit_msg";
}

sub print_msg {
	my $opt = shift;
	# display step name
	my $step = $opt->{step} =~ /reference/ ? $opt->{step}.'-'.$opt->{tool}.'#' : $opt->{step}.'#';
	(my $string = sprintf("%-25s ", $step)) =~ tr/ /./;
	$string =~ s/#/ /;
	print " $string";
	# display exit status
	if ($opt->{complete}) {
		print " complete\n";
	}
	else {
		if (exists $opt->{failed}) {
			printf(" kmers %s", join(', ',@{$opt->{failed}}));
			delete $opt->{failed};
		}
		print " failed\n";
		$opt->{display} = 0 unless ($opt->{step} =~ /reference/); # no need to print further msg if a step failed
	}
}

sub printOutMessage{
	my $message = shift;
	warn "ERROR:\n\t$message\n\n";
	pod2usage;
	exit 1;
}

sub create_report {
	my ($orig_report_folder, $exclude, $workflow_folder, $report_db_folder, $report) = @_ ;
	# Creates folder
	`rsync -r --exclude .svn --exclude $exclude $orig_report_folder $workflow_folder`;
	if ($? != 0) {
		croak "Error in report folder creation.\n" ;
	}
	# Writes master file
	open (REPORT, ">$report_db_folder/report_db.json") or croak "Can't create file $report_db_folder/report_db.json\n";
	print REPORT $report->get_json();
	close REPORT;
}

sub write_shell {
	my $opt  = shift;
	my $file = $opt->{scripts}->{$opt->{step}} || "$opt->{outdir}/$opt->{step}.sh";
	if (ref($file) eq 'ARRAY' && $#$file > 0) {
		my @lines = split("\n", $opt->{cmd});
		foreach my $i (0..$#lines) {
			open(SH, ">$file->[$i]") or croak "Can't open file $file->[$i]";
			print SH "#!/bin/csh\n";
			print SH $opt->{env}->{$opt->{step}.'_env'}."\n" if ($opt->{env}->{$opt->{step}.'_env'}); # get setup env command from drap.cfg for specific step
			print SH $lines[$i]."\n";
			close SH;
		}
	} else {
		$file = $file->[0] if (ref($file) eq 'ARRAY');
		open(SH, ">$file") or croak "Can't open file $file";
		print SH "#!/bin/csh\n" unless ($opt->{step} eq 'dbg' && $opt->{dbg} eq 'oases'); # because oases jobs will be run using the qarray command or separately in local mode
		if ($opt->{env}->{$opt->{step}.'_env'}) { # get setup env command from drap.cfg for specific step
			if ($opt->{step} eq 'dbg' && $opt->{dbg} eq 'oases') {
				map { print SH $opt->{env}->{$opt->{step}.'_env'}."; $_\n" } split("\n", $opt->{cmd});
			} else {
				print SH $opt->{env}->{$opt->{step}.'_env'}."\n";
				print SH $opt->{cmd};
			}
		} else {
			print SH $opt->{cmd};
		}
		close SH;
	}
	delete $opt->{cmd};
}

sub clean_directories {
	my ($opt, $index, $file) = @_;
	my $clean = 0;
	my @cleaned;
	my @dirs;
	my @mandatory_dirs = exists $opt->{meta} ? @{$opt->{dir_list}} : (@{$opt->{dir_list}}, glob("$opt->{dir_list}->[0]/k*"));
	foreach my $dir (@{$opt->{dir_list}}){
		if ($dir eq $opt->{dir_list}->[$index] || $clean) {
			next unless (-e $dir); # reference directories don't exist if reference not provided
			my $find_cmd = $file ? "find $dir -newer $file" : "find $dir";
			foreach (split("\n", process_cmd(1, $find_cmd))) {
				next if (find_in($_, \@mandatory_dirs)); # exclude mandatory dirs (modification date is not creation date)
				(-d) ? push(@dirs, $_) : unlink();
			}
			push(@cleaned, $dir);
			$clean = 1;
		}
	}
	map { rmdir() } sort { length($b) <=> length($a) } @dirs;
	if ($clean) {
		return sprintf(" director%s %s cleaned\n", $#cleaned > 0 ? 'ies' : 'y', join(', ', @cleaned));
	}
	return ""; # avoid the sub to return 0
}

sub clean_directory {
	my $root_dir    = shift;
	my $rm_root_dir = shift;
	my @dirs;
	find(sub { (-d) ? push(@dirs, $File::Find::name) : unlink() }, $root_dir);
	map { rmdir() unless (!$rm_root_dir && $_ eq $root_dir) } sort { length($b) <=> length($a) } @dirs;
	return " directory $root_dir cleaned\n";
}

sub clean_rmbt_directory {
	my $opt = shift;
	my $rmbt_directory = shift;
	my @deleted;
	foreach my $reference (keys %{$opt->{rmbt_status}}) {
		foreach my $alignR1 (keys %{$opt->{rmbt_status}->{$reference}}) {
			next if ($opt->{rmbt_status}->{$reference}->{$alignR1}->{status} == 1);
			foreach my $extension (qw(bam.flagstat bam sai Log.* SJ.* cmd log)) {
				unlink(glob("$rmbt_directory/*.".$opt->{rmbt_status}->{$reference}->{$alignR1}->{id}.".$extension"));
			}
			unlink(glob("$rmbt_directory/err_log_".$opt->{rmbt_status}->{$reference}->{$alignR1}->{id}."/*"));
			rmdir("$rmbt_directory/err_log_".$opt->{rmbt_status}->{$reference}->{$alignR1}->{id});
			push(@deleted, $opt->{rmbt_status}->{$reference}->{$alignR1}->{id});
			delete($opt->{rmbt_status}->{$reference}->{$alignR1});
		}
	}
	if (@deleted) {
		return sprintf(" directory %s cleaned: removed rmbt process Id%s: %s\n", $rmbt_directory, $#deleted > 0 ? 's' : '', join(', ', @deleted));
	}
	return ""; # avoid the sub to return 0
}

sub clean_rmbt_reference {
	my $opt = shift;
	my $rmbt_directory = shift;
	my $reference = shift;

	foreach my $extension (qw(cmd log)) {
		unlink(glob("$rmbt_directory/*.".$opt->{rmbt_ref_id}->{$reference}.".$extension"));
	}
	unlink(glob("$rmbt_directory/err_log_".$opt->{rmbt_ref_id}->{$reference}."/*"));
	rmdir("$rmbt_directory/err_log_".$opt->{rmbt_ref_id}->{$reference});
	if ($reference =~ m|/STAR_|) {
		unlink(glob($reference.'_'.$opt->{rmbt_ref_id}->{$reference}.'/*'));
		rmdir($reference.'_'.$opt->{rmbt_ref_id}->{$reference});
	} else {
		foreach my $extension (qw(amb ann bwt pac sa)) {
			unlink(glob("$reference.$extension"));
		}
	}
	delete($opt->{rmbt_ref_id}->{$reference});
	return " directory $rmbt_directory cleaned: removed rmbt reference: $reference\n";
}

sub valid_R1_R2 {
	my @R1 = split(//,shift);
	my @R2 = split(//,shift);
	my ($d, $i);
	return undef unless ($#R1 == $#R2);
	map { unless ($R1[$_] eq $R2[$_]) { $d++; $i=$_ } } 0..$#R1;
	return $R1[$i] == 1 && $R2[$i] == 2 && $d == 1 ? $i : undef;
}

sub gunzip {
	my $gzip = shift;
	$gzip =~ s/\.gz$//;
	return $gzip;
}

sub gzip {
	my $gzip = shift;
	$gzip .= $gzip =~ /\.gz$/ ? '' : '.gz';
	return $gzip;
}

sub find_in {
	my $string = shift;
	my $array  = shift;
	map { return 1 if ($string eq $_) } @$array;
	return 0;
}

sub is_prot_fasta {
	my $fa = shift;
	my $c  = 0;
	open (FA, $fa) or croak "Can't open file $fa";
	while (<FA>) {
		next if /^>/;
		$c++;
		if (m/[EFIJLOPQXZ]/i) {
			return 1;
		}
		last if $c > 1000;
	}
	close FA;
	return 0;
}

sub get_more_recent_dir { get_more_recent_file(@_) };

sub preprocess_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j1-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j1-*.e*');
	unless ($opt->{no_trim} || -e "$opt->{outdir}/.trim.over") {
		return [0, 1, "$opt->{outdir}/trim.log"];
	}
	return [0, 2, "$opt->{outdir}/filter.log"] unless (-e "$opt->{outdir}/.filter.over");
	unless ($opt->{no_norm}) {
		return [0, 3, $err] unless ($opt->{norm_merge_only} || -e "$opt->{outdir}/.normalize.over");
		return [0, 3, $err] unless ($opt->{pool} == 1 || $opt->{no_norm_merge} || -e "$opt->{outdir}/.merge.normalize.over");
	}
	foreach my $R (qw(R1 R2)) {
		next unless (@{$opt->{$R}});
		for (my $i = 0 ; $i < scalar(@{$opt->{$R}}) ; $i++) {
			return [0, 'empty', $opt->{"trim$R"}->[$i]] unless (!defined $opt->{"trim$R"} || -e $opt->{"trim$R"}->[$i]);
			return [0, 'empty', $opt->{"norm$R"}->[$i]] unless (!defined $opt->{"norm$R"} || -e $opt->{"norm$R"}->[$i]);
		}
		return [0, 'empty', $opt->{"merge$R"}] unless (!defined $opt->{"merge$R"} || -e $opt->{"merge$R"});
	}
	if ($opt->{dbg} eq 'oases') {
		return [0, 4, $log] if (count_inside_file($log, '^\[.+\] \d+ sequences found$') == 0);
		return [0, 4, $log] if (count_inside_file($log, '^\[.+\] Done$') == 0);
		return [0, 'empty', "$opt->{outdir}/Sequences"] unless (-s "$opt->{outdir}/Sequences");
	}
	return [0, 5, $err] if (count_inside_file($err, 'ERROR') > 0);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub dbg_complete {
	my $opt = shift;
	if ($opt->{dbg} eq 'oases') {
		my $log = get_more_recent_file("$opt->{outdir}/err_log", "j2-*.o*.$opt->{task}");
		my $err = get_more_recent_file("$opt->{outdir}/err_log", "j2-*.e*.$opt->{task}");
		return [0, 1, $log] if (count_inside_file($log, "^\\[.+\\] Exporting transcripts to $opt->{dir_list}->[0]/k$opt->{current_kmer}/transcripts.fa") == 0);
		return [0, 2, $log] if (count_inside_file($log, "^\\[.+\\] Finished extracting transcripts") == 0);
		# An empty transcripts.fa for a specific kmer could not be considerated as an error
		# return [0, 'empty', "$opt->{dir_list}->[0]/k$opt->{current_kmer}/transcripts.fa"] unless (-s "$opt->{dir_list}->[0]/k$opt->{current_kmer}/transcripts.fa");
		return [0, 'error', $err] unless (-z $err);
		return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log), $opt->{task}) == 1);
		$opt->{kmer_status}->{$opt->{current_kmer}} = 1;
	} else {
		my $log = get_more_recent_file("$opt->{outdir}/err_log", "j2-*.o*");
		my $err = get_more_recent_file("$opt->{outdir}/err_log", "j2-*.e*");
		return [0, 3, $log] if (count_inside_file($log, "^Butterfly assemblies are written to ") == 0);
		return [0, 'empty', "$opt->{dir_list}->[0]/Trinity.fasta"] unless (-s "$opt->{dir_list}->[0]/Trinity.fasta");
		return [0, 'error', $err] unless (-z $err);
		return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	}
	return [1];
}

sub merge_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j3-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j3-*.e*');
	if ($opt->{dbg} eq 'oases') {
		foreach my $kmer_dir (glob("$opt->{dir_list}->[0]/k*")) {
			my $kmer = basename($kmer_dir);
			next if (-z "$opt->{dir_list}->[0]/$kmer/transcripts.fa"); # specific kmer could generate empty transcripts.fa file
			return [0, 1, $log] if (count_inside_file($log, "^Clean locus in $kmer") == 0);
			return [0, 2, $err] if (count_inside_file($err, "$opt->{dir_list}->[0]/$kmer/seqclean, without a detectable error") == 0);
			return [0, 'empty', "$opt->{dir_list}->[0]/$kmer/transcripts.best.80\%Longest.fa"] unless (-s "$opt->{dir_list}->[0]/$kmer/transcripts.best.80\%Longest.fa");
			return [0, 'empty', "$opt->{dir_list}->[0]/$kmer/transcripts.best.80\%Longest.fa.cln"] unless (-s "$opt->{dir_list}->[0]/$kmer/transcripts.best.80\%Longest.fa.cln");
		}
	} else {
		return [0, 'error', $err] unless (-z $err);
	}
	return [0, 'empty',"$opt->{dir_list}->[1]/transcripts.fa" ] unless (-s "$opt->{dir_list}->[1]/transcripts.fa");
	return [0, 'empty', "$opt->{dir_list}->[1]/transcripts.raw.fa"] unless (-s "$opt->{dir_list}->[1]/transcripts.raw.fa");
	return [0, 'empty', "$opt->{dir_list}->[1]/transcripts.fa.chimFilter.log"] unless (-e "$opt->{dir_list}->[1]/transcripts.fa.chimFilter.log");
	my $nb_trans_raw = count_inside_file("$opt->{dir_list}->[1]/transcripts.raw.fa", '^>');
	my $nb_trans = count_inside_file("$opt->{dir_list}->[1]/transcripts.fa", '^>');
	my $nb_discarded = count_inside_file("$opt->{dir_list}->[1]/transcripts.fa.chimFilter.log", 'discarded');
	return [0, 3, $opt->{dir_list}->[1]] unless ($nb_trans == $nb_trans_raw-$nb_discarded);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub clustering_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j4-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j4-*.e*');
	return [0, 'cd-hit', "$opt->{dir_list}->[2]/all_dbg.fa.cd-hit.log"] if (!-e "$opt->{dir_list}->[2]/all_dbg.fa.cd-hit.log" || count_inside_file("$opt->{dir_list}->[2]/all_dbg.fa.cd-hit.log", '^program completed !') == 0);
	return [0, 'empty', "$opt->{dir_list}->[2]/all_dbg.fa"] unless (-s "$opt->{dir_list}->[2]/all_dbg.fa");
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub asm_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j5-*.o*');
	return [0, 'empty', "$opt->{dir_list}->[2]/all_contigs.raw.fa"] unless (-s "$opt->{dir_list}->[2]/all_contigs.raw.fa");
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub post_asm_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j6-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j6-*.e*');
	return [0, 'empty', "$opt->{dir_list}->[3]/all_contigs.fa"] unless (-s "$opt->{dir_list}->[3]/all_contigs.fa");
	return [0, 'empty', "$opt->{dir_list}->[3]/all_contigs.fa.vecFilter.log"] unless (-e "$opt->{dir_list}->[3]/all_contigs.fa.vecFilter.log");
	my $nb_ctgs_pre_filter;
	if ($opt->{type} eq 'orf') {
		return [0, 1, $err] if (count_inside_file($err, 'transdecoder is finished.') == 0);
		return [0, 'empty', "$opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.fa"] unless (-s "$opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.fa");
		my $nb_ctgs_raw = count_inside_file("$opt->{dir_list}->[2]/all_contigs.raw.fa", '^>');
		my $nb_contigs_after_cut = 0;
		open (LOG, "$opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.log") or return [0, 'open', "$opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.log"];
		my $header = <LOG>;
		map { /^.+\t(\d+)\t.+/ && $1 ? $nb_contigs_after_cut += $1 : $nb_contigs_after_cut++ } <LOG>;
		close LOG;
		$nb_ctgs_pre_filter = count_inside_file("$opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.fa", '^>');
		return [0, 2, $opt->{dir_list}->[3]] unless ($nb_contigs_after_cut == $nb_ctgs_pre_filter);
	} else {
		$nb_ctgs_pre_filter = count_inside_file("$opt->{dir_list}->[2]/all_contigs.raw.fa", '^>');
	}
	my $nb_contigs = count_inside_file("$opt->{dir_list}->[3]/all_contigs.fa", '^>');
	my $nb_discarded = count_inside_file("$opt->{dir_list}->[3]/all_contigs.fa.vecFilter.log", 'ERR: too (short|noisy)');
	return [0, 3, $opt->{dir_list}->[3]] unless ($nb_contigs == $nb_ctgs_pre_filter-$nb_discarded);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub rmbt_editing_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j7-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j7-*.e*');
	map { unlink($_) } glob("$opt->{dir_list}->[4]/*_pass.over");

	my $rmbt_checking = check_rmbt_directory($opt, $opt->{dir_list}->[4]);

	my $passes = {
		'first' => \@{["$opt->{dir_list}->[3]/all_contigs.fa", "$opt->{dir_list}->[4]/all_contigs.first_pass.raw.fa", "$opt->{dir_list}->[4]/all_contigs.first_pass.fa.cd-hit.log"]},
		'second' => \@{["$opt->{dir_list}->[4]/all_contigs.first_pass.fa", "$opt->{dir_list}->[4]/all_contigs.second_pass.raw.fa", "$opt->{dir_list}->[4]/all_contigs.second_pass.fa.cd-hit.log"]}
	};
	foreach my $pass ('first', 'second') {
		my $valid_rmbt = 0;
		foreach my $alignR1 (keys %{$opt->{rmbt_status}->{$opt->{rmbt_ref}->{"editing_$pass"}}}) {
			$valid_rmbt++ if ($opt->{rmbt_status}->{$opt->{rmbt_ref}->{"editing_$pass"}}->{$alignR1}->{status} == 1);
		}
		return [0, 'empty', "$passes->{$pass}->[0]"] unless (-e "$passes->{$pass}->[0]");
		return [0, 'empty', "$passes->{$pass}->[1]"] unless (-e "$passes->{$pass}->[1]");
		my $nb_contigs_pre = count_inside_file($passes->{$pass}->[0], '^>');
		my $nb_contigs_post = count_inside_file($passes->{$pass}->[1], '^>');
		return [0, 1, $pass] unless ($valid_rmbt == scalar(@{$opt->{alignR1}}));
		return [0, 'empty', "$opt->{dir_list}->[4]/all_contigs.${pass}_pass.raw.fa.samCorrectVariation.log"] unless (-e "$opt->{dir_list}->[4]/all_contigs.${pass}_pass.raw.fa.samCorrectVariation.log");
		return [0, 2, "$passes->{$pass}->[0] and $passes->{$pass}->[1]"] unless ($nb_contigs_pre == $nb_contigs_post);
		return [0, 'cd-hit', $passes->{$pass}->[2]] if (!-e $passes->{$pass}->[2] || count_inside_file($passes->{$pass}->[2], '^program completed !') == 0);
		return [0, 'empty', "$opt->{dir_list}->[4]/all_contigs.${pass}_pass.fa"] unless (-s "$opt->{dir_list}->[4]/all_contigs.${pass}_pass.fa");
		system("touch $opt->{dir_list}->[4]/${pass}_pass.over");
	}
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub rmbt_filtering_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j8-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j8-*.e*');
	my $rmbt_checking = check_rmbt_directory($opt, $opt->{dir_list}->[5]);
	return $rmbt_checking unless ($rmbt_checking->[0] == 1);
	my $valid_rmbt = 0;
	foreach my $alignR1 (keys %{$opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}}) {
		$valid_rmbt++ if ($opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}->{$alignR1}->{status} == 1);
	}
	return [0, 1] unless ($valid_rmbt == scalar(@{$opt->{alignR1}}));
	return [0, 2, $err] if (count_inside_file($err, 'Parsing BAM header...') != scalar(@{[glob("$opt->{dir_list}->[5]/*.bam")]}));
	return [0, 3, $err] if (count_inside_file($err, 'COMPLETED: Processed \d+ mapped fragments') == 0);
	return [0, 'empty', "$opt->{dir_list}->[5]/results.xprs"] unless (-s "$opt->{dir_list}->[5]/results.xprs");
	my $lowest_fpkm = $opt->{coverages}->[0];
	if ($opt->{type} eq 'orf') {
		my $lowest_fpkm_coding_file = "$opt->{dir_list}->[5]/coding_transcripts_fpkm_$lowest_fpkm.fa";
		return [0, 'empty', $lowest_fpkm_coding_file] unless (-s $lowest_fpkm_coding_file);
	}
	my $lowest_fpkm_file = "$opt->{dir_list}->[5]/transcripts_fpkm_$lowest_fpkm.fa";
	return [0, 'empty', $lowest_fpkm_file] unless (-s $lowest_fpkm_file);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub postprocess_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j9-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j9-*.e*');
	unless ($opt->{no_rate}) {
		return [0, 'transrate', "$opt->{outdir}/00-ASSEMBLY_RATING"] if (count_inside_file("$opt->{outdir}/00-ASSEMBLY_RATING", 'Writing analysis results to assemblies.csv') == 0);
	}
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{check_postprocess} && !$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub reference_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j10-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j10-*.e*');
	my $run_dir = $opt->{tool} eq 'exonerate' ? 6 : 7;
	my $tmp_dir = get_more_recent_dir("$opt->{dir_list}->[$run_dir]", "tmp_*");
	return [0, 1, $tmp_dir] if (-e $tmp_dir); # runBlat and runExonerate cleans tmp dir when all tasks ended successfully
	my $out_file = get_more_recent_file("$opt->{dir_list}->[$run_dir]", "*.$opt->{tool}.best.tsv");
	return [0, 'empty', $out_file] unless (-e $out_file);
	unless ($opt->{local}) {
		foreach my $run_job_id (find_inside_file($log, 'Your job (\d+) .* has been submitted')){
			return [0, 'qacct', $run_job_id] if (qacct_status($run_job_id) == 1);
		}
	}
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_merge_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j1-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j1-*.e*');
	return [0, 'empty', "$opt->{dir_list}->[0]/all_conditions_contigs.fa"] unless (-s "$opt->{dir_list}->[0]/all_conditions_contigs.fa");
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_longest_orf_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j2-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j2-*.e*');
	return [0, 'empty', "$opt->{dir_list}->[0]/all_conditions_contigs_longest_orf.faa"] unless (-s "$opt->{dir_list}->[0]/all_conditions_contigs_longest_orf.faa");
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_cluster_orf_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j3-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j3-*.e*');
	return [0, 'cd-hit', "$opt->{dir_list}->[0]/cd-hit_orfs.faa.log"] if (!-e "$opt->{dir_list}->[0]/cd-hit_orfs.faa.log" || count_inside_file("$opt->{dir_list}->[0]/cd-hit_orfs.faa.log", '^program completed !') == 0);
	return [0, 'empty', "$opt->{dir_list}->[0]/cd-hit_orfs.faa"] unless (-s "$opt->{dir_list}->[0]/cd-hit_orfs.faa");
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_longest_contig_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j4-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j4-*.e*');
	return [0, 'empty', "$opt->{dir_list}->[1]/cd-hit_contigs.fa"] unless (-s "$opt->{dir_list}->[1]/cd-hit_contigs.fa");
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_cluster_contig_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j5-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j5-*.e*');
	return [0, 'cd-hit', "$opt->{dir_list}->[1]/meta_contigs.fa.log"] if (!-e "$opt->{dir_list}->[1]/meta_contigs.fa.log" || count_inside_file("$opt->{dir_list}->[1]/meta_contigs.fa.log", '^program completed !') == 0);
	return [0, 'empty', "$opt->{dir_list}->[1]/meta_contigs.fa"] unless (-s "$opt->{dir_list}->[1]/meta_contigs.fa");
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_index_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j6-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j6-*.e*');
	if ($opt->{mapper} eq 'bwa') {
		return [0, 'empty', "$opt->{dir_list}->[1]/meta_contigs.fa.sa"] unless (-s "$opt->{dir_list}->[1]/meta_contigs.fa.sa");
	} else {
		return [0, 'empty', "$opt->{dir_list}->[2]/STAR_*/Genome"] unless (glob("$opt->{dir_list}->[2]/STAR_*/Genome"));
	}
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_rmbt_complete {
	my $opt = shift;
	my $rmbt_checking = check_rmbt_directory($opt, $opt->{dir_list}->[2]);
	return $rmbt_checking unless ($rmbt_checking->[0] == 1);
	my $valid_rmbt = 0;
	foreach my $alignR1 (keys %{$opt->{rmbt_status}->{$opt->{rmbt_ref}->{meta_rmbt}}}) {
		$valid_rmbt++ if ($opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}->{$alignR1}->{status} == 1);
	}
	return [1];
}

sub meta_filter_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j8-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j8-*.e*');
	return [0, 1, $err] if (count_inside_file($err, 'Parsing BAM header...') != scalar(@{[glob("$opt->{dir_list}->[2]/*.bam")]}));
	return [0, 2, $err] if (count_inside_file($err, 'COMPLETED: Processed \d+ mapped fragments') == 0);
	return [0, 3, $err] if (count_inside_file($err, 'transdecoder is finished.') == 0);
	return [0, 'empty', "$opt->{dir_list}->[3]/results.xprs"] unless (-s "$opt->{dir_list}->[3]/results.xprs");
	return [0, 'empty', "$opt->{dir_list}->[3]/meta_contigs.fa.transdecoder.bed"] unless (-s "$opt->{dir_list}->[3]/meta_contigs.fa.transdecoder.bed");
	my $lowest_fpkm = $opt->{coverages}->[0];
	if ($opt->{type} eq 'orf') {
		my $lowest_fpkm_coding_file = "$opt->{dir_list}->[3]/coding_transcripts_fpkm_$lowest_fpkm.fa";
		return [0, 'empty', $lowest_fpkm_coding_file] unless (-s $lowest_fpkm_coding_file);
	}
	my $lowest_fpkm_file = "$opt->{dir_list}->[3]/transcripts_fpkm_$lowest_fpkm.fa";
	return [0, 'empty', $lowest_fpkm_file] unless (-s $lowest_fpkm_file);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_postprocess_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j9-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j9-*.e*');
	return [0, 'transrate', "$opt->{outdir}/00-META-ASSEMBLY_RATING"] if (count_inside_file("$opt->{outdir}/00-META-ASSEMBLY_RATING", 'Writing analysis results to assemblies.csv') == 0);
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{check_postprocess} && !$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub meta_reference_complete {
	my $opt = shift;
	my $log = get_more_recent_file("$opt->{outdir}/err_log", 'j10-*.o*');
	my $err = get_more_recent_file("$opt->{outdir}/err_log", 'j10-*.e*');
	my $run_dir = $opt->{tool} eq 'exonerate' ? 4 : 5;
	my $tmp_dir = get_more_recent_dir("$opt->{dir_list}->[$run_dir]", "tmp_*");
	return [0, 1, $tmp_dir] if (-e $tmp_dir); # runBlat and runExonerate cleans tmp dir when all tasks ended successfully
	my $out_file = get_more_recent_file("$opt->{dir_list}->[$run_dir]", "*.$opt->{tool}.best.tsv");
	return [0, 'empty', $out_file] unless (-e $out_file);
	unless ($opt->{local}) {
		foreach my $run_job_id (find_inside_file($log, 'Your job (\d+) .* has been submitted')){
			return [0, 'qacct', $run_job_id] if (qacct_status($run_job_id) == 1);
		}
	}
	return [0, 'error', $err] unless (-z $err);
	return [0, 'qacct', get_job_id($log)] if (!$opt->{local} && qacct_status(get_job_id($log)) == 1);
	return [1];
}

sub count_inside_file {
	my $file = shift;
	my $pattern = shift;
	return scalar(find_inside_file($file, $pattern));
}

sub find_inside_file {
	my $file = shift;
	my $pattern = shift;
	my @found;
	open (IN, $file) or croak "Can't open file $file";
	map { m|$pattern| && push(@found, $1||$&) } <IN>;
	close IN;
	return @found;
}

sub check_rmbt_directory {
	my $opt = shift;
	my $rmbt_directory = shift;
	# check submissions during process
	unless ($opt->{local} || ! $opt->{restart}) {
		foreach my $log (glob("$rmbt_directory/runSam?e.*.log")) {
			foreach my $run_job_id (find_inside_file($log, 'Your job-array (\d+)\..+ .* has been submitted')){
				return [0, 'qacct', $run_job_id] if (qacct_status($run_job_id, 1) == 1);
				return [0, 'qacct', $run_job_id] if (qacct_status($run_job_id, 2) == 1);
			}
			foreach my $run_job_id (find_inside_file($log, 'Your job (\d+) .* has been submitted')){
				return [0, 'qacct', $run_job_id] if (qacct_status($run_job_id) == 1);
			}
		}
	}
	my $return;
	# need to check all rmbt processes but returns error message for the first encountered error only
	$opt->{rmbt_status} = get_cfg_param($opt->{outdir}, 'rmbt_status')||{} if ($opt->{restart} && exists $opt->{meta});
	foreach my $cmd (glob("$rmbt_directory/runSam?e.*.cmd")) {
		my ($rmbt_job_id) = $cmd =~ /\.(\d+)\.cmd$/;
		my $reference = shift(@{[find_inside_file($cmd, '-r\s+(\S+?)\s')]});
		$reference =~ s/_\d+$//; # remove id from STAR reference directory
		if (count_inside_file($cmd, '--index') > 0) {
			# index only process
			if (count_inside_file($cmd, '-m\s+star') > 0) {
				my $output = shift(@{[find_inside_file($cmd, '-o\s+(\S+?)\s')]});
				my $STAR_reference = $output.'/STAR_'.basename($reference);
				$opt->{rmbt_ref_id}->{$STAR_reference} = $rmbt_job_id;
			} else {
				$opt->{rmbt_ref_id}->{$reference} = $rmbt_job_id;
			}
		} else {
			my $rmbt_err_dir = "$rmbt_directory/err_log_$rmbt_job_id";
			my $bam = get_more_recent_file("$rmbt_directory", "*.$rmbt_job_id.bam");
			my $flagstat = $bam ? "$bam.flagstat" : '';
			if ($cmd =~ /runSampe\.\d+\.cmd$/) {
				$opt->{check_alignR1} = shift(@{[find_inside_file("$rmbt_directory/runSampe.$rmbt_job_id.cmd", '-1 (\S+?)\s')]});
				$opt->{check_alignR2} = shift(@{[find_inside_file("$rmbt_directory/runSampe.$rmbt_job_id.cmd", '-2 (\S+?)\s')]});
				$opt->{check_paired} = 1;
			} elsif ($cmd =~ /runSamse\.\d+\.cmd$/) {
				$opt->{check_alignR1} = shift(@{[find_inside_file("$rmbt_directory/runSamse.$rmbt_job_id.cmd", '-f (\S+?)\s')]});
				$opt->{check_paired} = 0;
			}
			next if ($opt->{restart} && exists $opt->{rmbt_status}->{$reference}->{$opt->{check_alignR1}} && $opt->{rmbt_status}->{$reference}->{$opt->{check_alignR1}}->{status} == 1);
			$opt->{rmbt_status}->{$reference}->{$opt->{check_alignR1}}->{id} = $rmbt_job_id;
			my $checking = check_align_process($opt, $rmbt_err_dir, $rmbt_job_id, $bam, $flagstat);
			if ($checking->[0] == 0) {
				$return ||= $checking;
				$opt->{rmbt_status}->{$reference}->{$opt->{check_alignR1}}->{status} = 0;
			} else {
				$opt->{rmbt_status}->{$reference}->{$opt->{check_alignR1}}->{status} = 1;
			}
		}
	}
	return $return||[1];
}

sub check_align_process {
	my ($opt, $rmbt_err_dir, $rmbt_job_id, $bam, $flagstat) = @_;
	return [0, 'empty', $bam||dirname($rmbt_err_dir)."/*.$rmbt_job_id.bam"] unless (-e $bam);
	return [0, 'empty', $flagstat||dirname($rmbt_err_dir)."/*.$rmbt_job_id.flagstat"] unless (-s $flagstat);
	my $nb_alignR1 = get_nb_reads($opt, $opt->{check_alignR1});
	if ($opt->{check_paired}) {
		my $nb_alignR2 = get_nb_reads($opt, $opt->{check_alignR2});
		my $nb_processed_sampe = 0;
		my $processed_sampe;
		if ($opt->{mapper} eq 'bwa') {
			my $sampe_err = get_more_recent_file($rmbt_err_dir, 'sampe.e*');
			my $aln_1_err = get_more_recent_file($rmbt_err_dir, 'aln.e*.1');
			my $aln_2_err = get_more_recent_file($rmbt_err_dir, 'aln.e*.2');
			return [0, "rmbt.2", $sampe_err] if (count_inside_file($sampe_err, '\[main\] Real time:') == 0);
			return [0, "rmbt.2", $aln_1_err] if (count_inside_file($aln_1_err, '\[main\] Real time:') == 0);
			return [0, "rmbt.2", $aln_2_err] if (count_inside_file($aln_2_err, '\[main\] Real time:') == 0);
			$nb_processed_sampe = pop(@{[find_inside_file($sampe_err, '^\[.+\] (\d+) sequences have been processed')]});
			my $nb_processed_aln_1 = pop(@{[find_inside_file($aln_1_err, '^\[.+\] (\d+) sequences have been processed')]});
			my $nb_processed_aln_2 = pop(@{[find_inside_file($aln_2_err, '^\[.+\] (\d+) sequences have been processed')]});
			my $nb_processed_flagstat = shift(@{[find_inside_file($flagstat, '^(\d+) \+ \d+ in total \(QC-passed reads \+ QC-failed reads\)')]});
			return [0, "rmbt.3", "$sampe_err and $aln_1_err"] unless ($nb_processed_sampe == $nb_processed_aln_1);
			return [0, "rmbt.3", "$sampe_err and $aln_2_err"] unless ($nb_processed_sampe == $nb_processed_aln_2);
			return [0, "rmbt.4", "$sampe_err and $flagstat"] unless ($nb_processed_sampe*2 == $nb_processed_flagstat);
			$processed_sampe = $sampe_err;
		} else {
			my $final_log = get_more_recent_file(dirname($rmbt_err_dir), "*.$rmbt_job_id.Log.final.out");
			return [0, "rmbt.5", dirname($rmbt_err_dir)] unless ($final_log);
			$nb_processed_sampe = shift(@{[find_inside_file($final_log, '^\s*Number of input reads[\s|]+(\d+)')]});
			$processed_sampe = $final_log;
		}
		return [0, "rmbt.6", "$processed_sampe and $opt->{check_alignR1}"] unless ($nb_processed_sampe == $nb_alignR1);
		return [0, "rmbt.6", "$processed_sampe and $opt->{check_alignR2}"] unless ($nb_processed_sampe == $nb_alignR2);
		unless ($opt->{local}) {
			foreach my $runSampe_job_id (find_inside_file(get_more_recent_file($rmbt_err_dir, 'runSampe.*.jid'), '(\d+)')){
				return [0, 'qacct', $runSampe_job_id] if (qacct_status($runSampe_job_id) == 1);
			}
		}
	} else {
		my ($nb_processed_samse, $nb_processed_sai, $processed_samse);
		if ($opt->{mapper} eq 'bwa') {
			my $samse_err = get_more_recent_file($rmbt_err_dir, 'samse.e*');
			return [0, "rmbt.7", $samse_err] unless (count_inside_file($samse_err, '\[main\] Real time:') == 2);
			my ($nb_processed_aln, $processed);
			open (ERR, $samse_err);
			while (<ERR>) {
				$processed = $1 if (/^\[.+\] (\d+) sequences have been processed/);
				if (/\[main\] Real time:/) {
					$nb_processed_aln ? $nb_processed_samse = $processed : $nb_processed_aln = $processed;
					$processed = 0;
				}
			}
			close ERR;
			return [0, "rmbt.8", $samse_err] unless ($nb_processed_samse == $nb_processed_aln);
			my $nb_processed_flagstat = shift(@{[find_inside_file($flagstat, '^(\d+) \+ \d+ in total \(QC-passed reads \+ QC-failed reads\)')]});
			return [0, "rmbt.9", $samse_err] unless ($nb_processed_samse == $nb_processed_flagstat);
			$processed_samse = $samse_err;
		} else {
			my $final_log = get_more_recent_file(dirname($rmbt_err_dir), "*.$rmbt_job_id.Log.final.out");
			return [0, "rmbt.5", dirname($rmbt_err_dir)] unless ($final_log);
			$nb_processed_samse = shift(@{[find_inside_file($final_log, '^\s*Number of input reads[\s|]+(\d+)')]});
			$processed_samse = $final_log;
		}
		return [0, "rmbt.6", "$processed_samse and $opt->{check_alignR1}"] unless ($nb_processed_samse == $nb_alignR1);
		unless ($opt->{local}) {
			my $runSamse_job_id = shift(@{[find_inside_file(get_more_recent_file($rmbt_err_dir, 'runSamse.*.jid'), '(\d+)')]});
			return [0, 'qacct', $runSamse_job_id] if (qacct_status($runSamse_job_id) == 1);
		}
	}
	return [1];
}

sub qacct_status {
	my $jid = shift;
	my $task = shift;
	my $username = getpwuid($<);
	my $qacct_cmd = sprintf("qacct -o $username -j %d%s", $jid, $task ? " -t $task" : '');
	printf("    Warning: no qacct checking for job %d%s\n", $jid, $task ? ".$task" : '') && return 0 if ($main::no_qacct);
	my $output = process_cmd(1, $qacct_cmd);
	while ($output =~ m/(?:failed|exit_status)\s+(\d+)/g) {
		return 1 if $1 != 0;
	}
	return 0;
}

sub get_job_id {
	my $log = shift;
	$log =~ /.+\.o(\d+)/;
	return $1;
}

sub get_scheduler_type {
	my ( $config_path ) = @_ ;
	my $scheduler_type = "local" ;
	if( -e $config_path && -r $config_path ){ # If it exist a config file
		my $config = new ConfigFile( $config_path );
		if( defined($config->{'SCHEDULER'}) && defined($config->{'SCHEDULER'}{'type'}) ){ # If the scheduler is set in config file
			$scheduler_type = $config->{'SCHEDULER'}{'type'} ;
		}
	}
	return $scheduler_type;
}

sub set_env_variables {
	my ( $config_path, $opt ) = @_ ;

	if( -e $config_path && -r $config_path ){ # If configuration file exists and if it is readable
		my $config = new ConfigFile( $config_path );
		foreach my $section ('ENV','DATABASE','SOFTWARE CONFIGURATION','SGE CONFIGURATION','SGE RESOURCES') {
			if( defined($config->{$section}) ){
				foreach my $key (keys %{$config->{$section}}) {
					while ($config->{$section}{$key} =~ /\$(\w+)/) {
						my $substitute = $1;
						my $value = $ENV{$substitute}||$config->{$section}{$substitute};
						croak "Died with error: unable to find value for variable \$$substitute in $config_path\n" unless ($value);
						$config->{$section}{$key} =~ s/\$$substitute/$value/g;
					}
					$ENV{$key} = $config->{$section}{$key};
					$opt->{env}->{$key} = $ENV{$key};
				}
			}
		}
		if( defined($config->{'SCHEDULER'}) ){
			foreach my $key (keys %{$config->{'SCHEDULER'}}) {
				$ENV{"scheduler_$key"} = $config->{'SCHEDULER'}{$key};
				$opt->{env}->{scheduler}->{$key} = $ENV{"scheduler_$key"};
			}
		}
	}
}

=head2 procedure set_extended_path

 Usage        : set_extended_path( $config_path )
 Function     : Adds in environment variable PATH the specified folders in configuration file (section "PATH") if it exists.
 Args         : [str] the configuration file path.

=cut
sub set_extended_path {
	my ( $config_path ) = @_ ;

	if( -e $config_path && -r $config_path ){ # If configuration file exists and if it is readable
		my $config = new ConfigFile( $config_path );
		if( defined($config->{'PATH'}) && scalar(@{$config->{'PATH'}}) > 0 ){ # If folders have been added to section PATH
			$ENV{PATH} = join( ':', @{$config->{'PATH'}} ).":$ENV{PATH}" ;
		}
	}
}

1;
