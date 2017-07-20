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
use Carp 'croak';
use File::Basename;
use File::Copy;
use File::Compare;
use File::Path qw(make_path remove_tree);
use POSIX qw(strftime);
use lib ("$FindBin::RealBin");
use AsmUtils;
use Report;

my $dir = shift;
my $opt = get_drap_config($dir) or croak("Unable to read drap config file .drap_conf.json inside $dir");

our $no_qacct = $opt->{no_qacct};
$opt->{display} = $opt->{restart};
my $resubmit = '';
my $clean_msg;
print "Check steps already executed:\n" if ($opt->{restart} && -f $opt->{submit_log});

# create directory tree
make_path("$dir/err_log", grep { basename($_) !~ /(exonerate|blat)/ || ($opt->{ref} && $opt->{$1}) } @{$opt->{dir_list}});

# rename existing shell scripts
foreach my $script (glob("$dir/*.sh")) {
	next if ($script =~ /\d+:\d+:\d+\.sh/); # previously renamed script
	my $mdate = strftime("%Y-%m-%d_%H:%M:%S", localtime((stat($script))[9]));
	my $name_no_ext = $dir.'/'.basename($script, '.sh');
	map { unlink() if (compare($script, $_) == 0) } glob("$name_no_ext.*-*-*_*:*:*.sh"); # remove previously renamed scripts identical to current one
	rename($script, "$name_no_ext.$mdate.sh");
}

# start report
my $orig_report_folder = $opt->{reportpath};
my $report_folder = $dir."/report";
my $report_db_folder = $report_folder."/database";
my $report = undef;
if( $opt->{restart} ){
	$report = Report::load_json( $report_db_folder."/report_db.json" );
} else {
	$report = new Report( "DRAP ".$opt->{dbg} );
}

my $nstep = 0;
# create Sequences file inside output directory
# write file 01.sh
my $preprocess_analysis = $report->get_or_create_analysis( 'Preprocess', 'Filter reads and convert to fasta.' );
$opt->{step} = $opt->{steps}->[$nstep];
$opt->{complete} = 1;
unless (step_complete($opt)) {
	$opt->{complete} = 0;
	my $local_restart = $opt->{restart};
	if ($opt->{restart}) {
		$clean_msg .= clean_directories($opt, 0);
		$opt->{kmer_status}	= map { $_, 0 } @{$opt->{kmers}} if (ref($opt->{kmer_status}) eq 'HASH');
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	$opt->{preprocessed}->{R1} = $opt->{R1}; # Path to the R1 file after current process
	$opt->{preprocessed}->{R2} = $opt->{paired} ? $opt->{R2} : ""; # Path to the R2 file after current process
	my @provided_pairs = ("R1");
	if ($opt->{paired}) {
		@provided_pairs = ("R1", "R2");
	}
	# Report reads statistics on raw data
	my $rawData_step = $preprocess_analysis->get_or_create_step('Raw data');
	unless (-f "$dir/.raw_data.over") {
		$opt->{cmd} = "find $dir -name \\*.unzip -exec \\rm -f {} \\;\n" if ($local_restart);
		my @all_fq = map {@{$opt->{preprocessed}->{$_}}} @provided_pairs;
		my @all_json = map {$report_db_folder."/".$rawData_step->get_or_create_metrics_filename( 'fastaAssembly', basename($_))} @all_fq;
		#### supprimer de @all_fq et @all_json les fichiers pour lesquels process OK
		$opt->{cmd} .= "parallel --joblog $dir/raw_data.log -j $opt->{env}->{n_cpu} -N 2 \\\n"
		              ." $opt->{binpath}/fastq_unzip.sh '{1}' $dir ';' \\\n"
		              ." $opt->{binpath}/fastaAssemblyMetrics.pl --input $dir/'{1/}'.unzip '>' $dir/'{1/.}'.metrics ';' \\\n"
		              ." $opt->{binpath}/fastaMetrics2json.pl $dir/'{1/.}'.metrics '>' '{2}' ';' \\\n";
		$opt->{cmd} .= " $opt->{binpath}/fastaMetrics2tsv.pl '{1}' $dir/'{1/.}'.metrics '>>' $dir/.processed_reads.tsv ';' \\\n" unless ($opt->{nb_frags});
		$opt->{cmd} .= " \\rm -f $dir/'{1/.}'.metrics \\\n"
		              .' ::: '.join(' ', map { $all_fq[$_], $all_json[$_] } (0 .. $#all_fq))."\n";
		$opt->{cmd} .= "if (`sed 1d $dir/raw_data.log | awk '\$7>0'` != '') then\n"
		               ."\techo 'Following raw data stats command(s) failed:'\n"
		               ."\tsed 1d $dir/raw_data.log | awk '\$7>0' | cut -f9\n"
		               ."\texit\nelse\n\t\\rm -f $dir/raw_data.log\n\ttouch $dir/.raw_data.over\nendif\n";
	}
	@{$opt->{preprocessed}->{R1}} = map { $dir.'/'.basename($_).'.unzip' } @{$opt->{R1}};
	@{$opt->{preprocessed}->{R2}} = map { $dir.'/'.basename($_).'.unzip' } @{$opt->{R2}} if ($opt->{paired});
	# Trim raw data using Trim galore
	my $min_length;
	unless ($opt->{no_trim}) {
		unless (-f "$dir/.trim.over") {
			# need to known the k-mer size used for the normalization step
			unless ($opt->{no_norm}) {
				($min_length) = $opt->{norm} =~ /(:?-k|--ksize|--KMER_SIZE)\s+(\d+)/;
				$min_length ||= $opt->{norm_source} eq 'trinity' ? 25 : 32;
			}
			$opt->{cmd} .= 'mkdir -p '.join(' ', map {"$dir/trim_".$_} (1 .. $opt->{pool}))."\n";
			$opt->{cmd} .= sprintf("parallel --joblog %s/trim.log -j %d -N %d \\\n trim_galore --no_report_file %s-o %s/trim_'{#}' %s %s ';' \\\n",
				$dir, $opt->{env}->{n_cpu}, $opt->{paired} ? 4 : 2,	$min_length ? "--length $min_length " : '', $dir, $opt->{trim},
				$opt->{paired} ? "--paired '{1}' '{2}'" : "'{1}'"
			);
			if ($opt->{paired}) {
				$opt->{cmd} .= " mv $dir/trim_'{#}'/'{1/}'_val_1.fq '{3}' ';' \\\n"
							  ." mv $dir/trim_'{#}'/'{2/}'_val_2.fq '{4}' ';' \\\n"
			} else {
				$opt->{cmd} .= " mv $dir/trim_'{#}'/'{1/}'_trimmed.fq '{2}' ';' \\\n"
			}
			$opt->{cmd} .= "\\rm -rf $dir/trim_'{#}' \\\n";
			$opt->{cmd} .= sprintf(" ::: %s\n",
				$opt->{paired} ? join(' ', map {$opt->{preprocessed}->{R1}->[$_], $opt->{preprocessed}->{R2}->[$_], gunzip($opt->{trimR1}->[$_]), gunzip($opt->{trimR2}->[$_])} (0 .. $opt->{pool}-1)) : join(" ", map {$opt->{preprocessed}->{R1}->[$_], gunzip($opt->{trimR1}->[$_])} (0 .. $opt->{pool}-1))
			);
			$opt->{cmd} .= "if (`sed 1d $dir/trim.log | awk '\$7>0'` != '') then\n"
			              ."\techo 'Following trimming command(s) failed:'\n"
			              ."\tsed 1d $dir/trim.log | awk '\$7>0' | cut -f9\n"
			              ."\texit\nelse\n\t\\rm -f $dir/trim.log\n\ttouch $dir/.trim.over\nendif\n";
		}
		@{$opt->{preprocessed}->{R1}} = map { gunzip($_) } @{$opt->{trimR1}};
		@{$opt->{preprocessed}->{R2}} = map { gunzip($_) } @{$opt->{trimR2}} if ($opt->{paired});
	}
	# Report reads statistics on trimmed data
	unless ($opt->{no_trim}) {
		my $trimGalore_step = $preprocess_analysis->get_or_create_step('Trim Galore', "Trim adapters, Trim the 3' end on quality, and remove sequences became shorter than 20bp.");
		unless (-f "$dir/.trimmed_data.over") {
			my @all_fq = map {@{$opt->{preprocessed}->{$_}}} @provided_pairs;
			my @all_json = map {$report_db_folder."/".$trimGalore_step->get_or_create_metrics_filename( 'fastaAssembly', basename($_))} @all_fq;
			$opt->{cmd} .= "parallel --joblog $dir/trimmed_data.log -j $opt->{env}->{n_cpu} -N 2 \\\n"
		                ." $opt->{binpath}/fastaAssemblyMetrics.pl --input '{1}' '>' '{1}'.metrics ';' \\\n"
		                ." $opt->{binpath}/fastaMetrics2json.pl '{1}'.metrics '>' '{2}' ';' \\\n"
		                ." $opt->{binpath}/fastaMetrics2tsv.pl '{1}'.gz '{1}'.metrics '>>' $dir/.processed_reads.tsv ';' \\\n"
		                ." \\rm -f '{1}'.metrics \\\n"
		                .' ::: '.join(' ', map { $all_fq[$_], $all_json[$_] } (0 .. $#all_fq))."\n";
			$opt->{cmd} .= "if (`sed 1d $dir/trimmed_data.log | awk '\$7>0'` != '') then\n"
		                ."\techo 'Following trimmed data stats command(s) failed:'\n"
		                ."\tsed 1d $dir/trimmed_data.log | awk '\$7>0' | cut -f9\n"
		                ."\texit\nelse\n\t\\rm -f $dir/trimmed_data.log\n\ttouch $dir/.trimmed_data.over\nendif\n";
		}
	}
	# Filter raw or trimmed data using filter_illumina
	unless (-f "$dir/.filter.over") {
		$opt->{cmd} .= sprintf("parallel --joblog %s/filter.log -j %d -N %d \\\n %s/filter_illumina %s %s-s '/' %s \\\n :::",
			$dir, $opt->{env}->{n_cpu}, $opt->{paired} ? 8 : 4,	$opt->{binpath}, $opt->{clean}, $min_length ? "-m $min_length " : '',
			$opt->{paired} ? "'{1}' '{2}' '{3}' '{4}' '{5}' '{6}' '{7}' '{8}'" : "'{1}' '{2}' '{3}' '{4}'"
		);
		for (my $i = 0 ; $i < scalar(@{$opt->{R1}}) ; $i++) {
			$opt->{cmd} .= $opt->{paired} ? " -i $opt->{preprocessed}->{R1}->[$i] -i $opt->{preprocessed}->{R2}->[$i]" : " -i $opt->{preprocessed}->{R1}->[$i]",
			$opt->{cmd} .= $opt->{paired} ? " -o $opt->{cleanR1}->[$i] -o $opt->{cleanR2}->[$i]" : " -o $opt->{cleanR1}->[$i]"
		}
		$opt->{cmd} .= "\nif (`sed 1d $dir/filter.log | awk '\$7>0'` != '') then\n"
		              ."\techo 'Following filtering command(s) failed:'\n"
		              ."\tsed 1d $dir/filter.log | awk '\$7>0' | cut -f9\n"
		              ."\texit\nelse\n\t\\rm -f $dir/filter.log $dir/*.unzip\n\ttouch $dir/.filter.over\nendif\n";
	}
	$opt->{preprocessed}->{R1} = $opt->{cleanR1};
	$opt->{preprocessed}->{R2} = $opt->{cleanR2} if ( $opt->{paired} );
	# Report reads statistics on filtered data
	my $filter_step = $preprocess_analysis->get_or_create_step('Filter Illumina');
	unless (-f "$dir/.filtered_data.over") {
		my @all_fq = map {@{$opt->{preprocessed}->{$_}}} @provided_pairs;
		my @all_json = map {$report_db_folder."/".$filter_step->get_or_create_metrics_filename( 'fastaAssembly', basename($_))} @all_fq;
		$opt->{cmd} .= "parallel --joblog $dir/filtered_data.log -j $opt->{env}->{n_cpu} -N 2 \\\n"
		              ." $opt->{binpath}/fastaAssemblyMetrics.pl --input '{1}' '>' '{1}'.metrics ';' \\\n"
		              ." $opt->{binpath}/fastaMetrics2json.pl '{1}'.metrics '>' '{2}' ';' \\\n"
		              ." \\rm -f '{1}'.metrics \\\n"
		              .' ::: '.join(' ', map { $all_fq[$_], $all_json[$_] } (0 .. $#all_fq))."\n";
		$opt->{cmd} .= "if (`sed 1d $dir/filtered_data.log | awk '\$7>0'` != '') then\n"
		              ."\techo 'Following filtered data stats command(s) failed:'\n"
		              ."\tsed 1d $dir/filtered_data.log | awk '\$7>0' | cut -f9\n"
		              ."\texit\nelse\n\t\\rm -f $dir/filtered_data.log\n\ttouch $dir/.filtered_data.over\nendif\n";
	}
	# Normalize filtered data using Trinity or Khmer
	unless ($opt->{no_norm}) {
		unless ($opt->{norm_merge_only} || -f "$dir/.normalize.over") {
			$opt->{cmd} .= "\\rm -f $dir/tmp_norm/*\n" if ($local_restart && $opt->{norm_src} eq 'trinity');
			for (my $i = 0 ; $i < scalar(@{$opt->{R1}}) ; $i++) {
				my $fastq = $opt->{paired} ? "$opt->{preprocessed}->{R1}->[$i] - $opt->{preprocessed}->{R2}->[$i]" : "$opt->{preprocessed}->{R1}->[$i]";
				if ($opt->{norm_src} eq 'trinity') {
					$opt->{cmd} .= sprintf("insilico_read_normalization.pl --seqType fq --JM %dG --CPU %d --output %s/tmp_norm %s %s%s\n",
						$opt->{norm_mem}, $opt->{env}->{n_cpu}, $dir, $opt->{norm},
						$opt->{paired} ? "--pairs_together --PARALLEL_STATS --left $opt->{preprocessed}->{R1}->[$i] --right $opt->{preprocessed}->{R2}->[$i]" : "--single $opt->{preprocessed}->{R1}->[$i]",
						$opt->{stranded} ? " --SS_lib_type $opt->{strand}" : ''
					);
					$opt->{cmd} .= "if (\$status != 0) then\n\techo Normalization failed for reads $fastq\n\texit\nendif\n";
					if ($opt->{paired}) {
						$opt->{cmd} .= sprintf("cp -L %s/tmp_norm/left.norm.fq %s\n", $dir, gunzip($opt->{normR1}->[$i]));
						$opt->{cmd} .= sprintf("cp -L %s/tmp_norm/right.norm.fq %s\n", $dir, gunzip($opt->{normR2}->[$i]));
					} else {
						$opt->{cmd} .= sprintf("cp -L %s/tmp_norm/single.norm.fq %s\n", $dir, gunzip($opt->{normR1}->[$i]));
					}
					$opt->{cmd} .= "\\rm -rf $dir/tmp_norm\n";
				} else {
					$opt->{cmd} .= sprintf("%s | normalize-by-median.py -M %s %s %s --output - - %s\n",
						$opt->{paired} ? "interleave-reads.py $opt->{preprocessed}->{R1}->[$i] $opt->{preprocessed}->{R2}->[$i]" : "cat $opt->{preprocessed}->{R1}->[$i]",
						$opt->{norm_mem}.'000000000', $opt->{paired} ? '--paired ' : '', $opt->{norm},
						$opt->{paired} ? '| split-paired-reads.py -1 '.gunzip($opt->{normR1}->[$i]).' -2 '.gunzip($opt->{normR2}->[$i]) : '> '.gunzip($opt->{normR1}->[$i])
					);
					$opt->{cmd} .= "if (\$status != 0) then\n\techo Normalization failed for reads $fastq\n\texit\nendif\n";
				}
			}
			$opt->{cmd} .= "touch $dir/.normalize.over\n";
		}
		unless ($opt->{norm_merge_only}) {
			@{$opt->{preprocessed}->{R1}} = map { gunzip($_) } @{$opt->{normR1}};
			@{$opt->{preprocessed}->{R2}} = map { gunzip($_) } @{$opt->{normR2}} if ($opt->{paired});
		}
		# Add a last normalization run after normalized fastq merge
		unless ($opt->{pool} == 1 || $opt->{no_norm_merge} || -f "$dir/.merge.normalize.over") {
			my $fastq = $opt->{paired} ? "$opt->{mergeR1} - $opt->{mergeR2}" : "$opt->{mergeR1}";
			my $mergeCleanR1 = $dir.'/'.basename(gunzip($opt->{mergeR1}),'.norm.fq').'.clean.fq';
			$opt->{cmd} .= sprintf("cat %s > %s\n", join(' ', @{$opt->{preprocessed}->{R1}}), $mergeCleanR1);
			my $mergeCleanR2 = $dir.'/'.basename(gunzip($opt->{mergeR2}),'.norm.fq').'.clean.fq' if ($opt->{paired});
			$opt->{cmd} .= sprintf("cat %s > %s\n", join(' ', @{$opt->{preprocessed}->{R2}}), $mergeCleanR2) if ($opt->{paired});
			if ($opt->{norm_src} eq 'trinity') {
				$opt->{cmd} .= sprintf("insilico_read_normalization.pl --seqType fq --JM %dG --CPU %d --output %s/tmp_norm %s %s%s\n",
					$opt->{norm_mem}, $opt->{env}->{n_cpu}, $dir, $opt->{norm},
					$opt->{paired} ? "--pairs_together --PARALLEL_STATS --left $mergeCleanR1 --right $mergeCleanR2" : "--single $mergeCleanR1",
					$opt->{stranded} ? " --SS_lib_type $opt->{strand}" : ''
				);
				$opt->{cmd} .= "if (\$status != 0) then\n\techo Normalization failed for reads $fastq\n\texit\nendif\n";
				if ($opt->{paired}) {
					$opt->{cmd} .= sprintf("cp -L %s/tmp_norm/left.norm.fq %s\n", $dir, gunzip($opt->{mergeR1}));
					$opt->{cmd} .= sprintf("cp -L %s/tmp_norm/right.norm.fq %s\n", $dir, gunzip($opt->{mergeR2}));
				} else {
					$opt->{cmd} .= sprintf("cp -L %s/tmp_norm/single.norm.fq %s\n", $dir, gunzip($opt->{mergeR1}));
				}
				$opt->{cmd} .= "\\rm -rf $dir/tmp_norm\n";
			} else {
				$opt->{cmd} .= sprintf("%s | normalize-by-median.py -M %s %s %s --output - - %s\n",
					$opt->{paired} ? "interleave-reads.py $mergeCleanR1 $mergeCleanR2" : "cat $mergeCleanR1",
					$opt->{norm_mem}.'000000000', $opt->{paired} ? '--paired ' : '', $opt->{norm},
					$opt->{paired} ? '| split-paired-reads.py -1 '.gunzip($opt->{mergeR1}).' -2 '.gunzip($opt->{mergeR2}) : '> '.gunzip($opt->{mergeR1})
				);
				$opt->{cmd} .= "if (\$status != 0) then\n\techo Normalization failed for reads $fastq\n\texit\nendif\n";
			}
			$opt->{cmd} .= "touch $dir/.merge.normalize.over\n";
		}
	}
	# Report reads statistics on normalized data	
	unless ($opt->{no_norm}) {	
		my $normalization_step = $preprocess_analysis->get_or_create_step('Normalization', 'Reduce dataset by decreasing sampling variation, discarding redundant data, and removing the majority of errors.');
		unless (-f "$dir/.normalized_data.over") {
				my @all_fq = map {@{$opt->{preprocessed}->{$_}}} @provided_pairs;
				push(@all_fq, $opt->{paired} ? (gunzip($opt->{mergeR1}), gunzip($opt->{mergeR2})) : gunzip($opt->{mergeR1})) unless ($opt->{pool} == 1 || $opt->{no_norm_merge});
				my @all_json = map {$report_db_folder."/".$normalization_step->get_or_create_metrics_filename( 'fastaAssembly', basename($_))} @all_fq;
				$opt->{cmd} .= "parallel --joblog $dir/normalized_data.log -j $opt->{env}->{n_cpu} -N 2 \\\n"
			                ." $opt->{binpath}/fastaAssemblyMetrics.pl --input '{1}' '>' '{1}'.metrics ';' \\\n"
			                ." $opt->{binpath}/fastaMetrics2json.pl '{1}'.metrics '>' '{2}' ';' \\\n"
			                ." $opt->{binpath}/fastaMetrics2tsv.pl '{1}'.gz '{1}'.metrics '>>' $dir/.processed_reads.tsv ';' \\\n"
			                ." \\rm -f '{1}'.metrics \\\n"
			                .' ::: '.join(' ', map { $all_fq[$_], $all_json[$_] } (0 .. $#all_fq))."\n";
				$opt->{cmd} .= "if (`sed 1d $dir/normalized_data.log | awk '\$7>0'` != '') then\n"
			                ."\techo 'Following normalized data stats command(s) failed:'\n"
			                ."\tsed 1d $dir/normalized_data.log | awk '\$7>0' | cut -f9\n"
			                ."\texit\nelse\n\t\\rm -f $dir/normalized_data.log\n\ttouch $dir/.normalized_data.over\nendif\n";
		}
	}
	if ($opt->{dbg} eq 'oases') {
		if ($opt->{paired}) {
			if ($opt->{no_norm} || $opt->{no_norm_merge} || $opt->{pool} == 1) {
				for (my $i = 0 ; $i < scalar(@{$opt->{R1}}) ; $i++) {
					$opt->{cmd} .= sprintf("interleave-reads.py %s %s >> %s/tmp_Sequences\n", $opt->{preprocessed}->{R1}->[$i], $opt->{preprocessed}->{R2}->[$i], $dir);
				}
			} else {
				$opt->{cmd} .= sprintf("interleave-reads.py %s %s >> %s/tmp_Sequences\n", 
					gunzip($opt->{mergeR1}), gunzip($opt->{mergeR2}), $dir
				);
			}
			$opt->{cmd} .= "cat $dir/tmp_Sequences";
		} else {
			if ($opt->{no_norm} || $opt->{no_norm_merge} || $opt->{pool} == 1) {
				$opt->{cmd} .= sprintf("cat %s", join(' ', @{$opt->{preprocessed}->{R1}}));
			} else {
				$opt->{cmd} .= sprintf("cat %s", gunzip($opt->{mergeR1}));
			}
		}
		$opt->{cmd} .= sprintf(" | velveth %s 27 %s -fastq - -noHash\n", $dir, $opt->{paired} ? '-shortPaired' : '-short');
		$opt->{cmd} .= "mv $dir/Log $opt->{dir_list}->[0]\n";
		$opt->{cmd} .= "if (-e $dir/Sequences && ! -z $dir/Sequences) then\n"
		              ."\tfind $dir -name tmp_Sequences -exec \\rm -f {} \\;\n"
		              ."endif\n";
	}
	$opt->{cmd} .= "find $dir -name \\*.clean.fq -exec \\rm -f {} \\;\n" unless ($opt->{dbg} eq 'trinity' && $opt->{no_norm});
	$opt->{cmd} .= "find $dir -name \\*.fq -print | xargs parallel -j $opt->{env}->{n_cpu} gzip '{}' ::: \n";
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 1);

# run DBG assembler
# write file 02.sh
my $dbg_analysis = $report->get_or_create_analysis('Assembly', 'Assemble sequences with De Bruijn Graph assembler.');
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
if ($opt->{dbg} eq 'oases') {
	$opt->{task} = 1;
	my $assembly_step = $dbg_analysis->get_or_create_step('Oases', 'Assemble transcriptome with following k-mers size : '.join(', ', @{$opt->{kmers}}));
	foreach my $kmer (@{$opt->{kmers}}) {
		next if $opt->{kmer_status}->{$kmer} == 1;
		$opt->{current_kmer} = $kmer;
		make_path("$opt->{dir_list}->[0]/k$kmer");
		unless (-f "$opt->{dir_list}->[0]/k$kmer/transcripts.fa" && step_complete($opt)) {
			$opt->{cmd} .= sprintf("echo VELVETH_START; velveth %s/k%s %s %s -reuse_Sequences; echo VELVETG_START; velvetg %s/k%s -read_trkg yes -min_contig_lgth 100 -cov_cutoff 4; echo OASES_START; oases %s/k%s -cov_cutoff 4; ",
				$opt->{dir_list}->[0], $kmer, $kmer, $opt->{stranded} ? '-strand_specific' : '', $opt->{dir_list}->[0], $kmer, $opt->{dir_list}->[0], $kmer
			);
			my $metrics_file = $report_db_folder."/".$assembly_step->get_or_create_metrics_filename( 'fastaAssembly', 'k-mer_'.$kmer );
			$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[0]/k$kmer/transcripts.fa > $opt->{dir_list}->[0]/k$kmer/assembler_metrics.tsv ; $opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[0]/k$kmer/assembler_metrics.tsv > $metrics_file ; \\rm -f $opt->{dir_list}->[0]/k$kmer/assembler_metrics.tsv\n";
			$clean_msg .= clean_directory("$opt->{dir_list}->[0]/k$kmer", 0) if ($opt->{restart});
			push(@{$opt->{failed}}, "k$kmer") if ($opt->{restart});
			$opt->{complete} = 0;
		}
		$opt->{task}++;
		symlink("../../Sequences", "$opt->{dir_list}->[0]/k$kmer/Sequences");
	}
} else {
	my $assembly_step = $dbg_analysis->get_or_create_step('Trinity', 'Assemble transcriptome.');
	unless (-f "$opt->{dir_list}->[0]/Trinity.fasta" && step_complete($opt)) {
		my $input_reads;
		if ($opt->{paired}) {
			if ($opt->{no_norm} || $opt->{no_norm_merge} || $opt->{pool} == 1) {
				$input_reads = sprintf("--left %s --right %s", 
					join(',', map { gzip($_) } @{$opt->{preprocessed}->{R1}}),
					join(',', map { gzip($_) } @{$opt->{preprocessed}->{R2}})
				)
			} else {
				$input_reads = sprintf("--left %s --right %s", $opt->{mergeR1}, $opt->{mergeR2});
			}
		} else {
			if ($opt->{no_norm} || $opt->{no_norm_merge} || $opt->{pool} == 1) {
				$input_reads = sprintf("--single %s", join(',', map { gzip($_) } @{$opt->{preprocessed}->{R1}}));
			} else {
				$input_reads = sprintf("--single %s", $opt->{mergeR1});
			}
		}
		if ($opt->{local}) {
			$opt->{cmd} .= sprintf("Trinity --no_cleanup --seqType fq --max_memory %dG --bflyHeapSpaceMax 4G --CPU %d --output %s --no_normalize_reads %s%s\n",
				$opt->{dbg_mem}, $opt->{env}->{n_cpu}, $opt->{dir_list}->[0], 
				$input_reads, $opt->{stranded} ? " --SS_lib_type $opt->{strand}" : ''
			);
		} else {
			$opt->{cmd} .= sprintf("Trinity --no_cleanup --seqType fq --max_memory %dG --bflyHeapSpaceMax 4G --CPU %d --no_distributed_trinity_exec --output %s --no_normalize_reads %s%s\n",
				$opt->{dbg_mem}, $opt->{env}->{n_cpu}, $opt->{dir_list}->[0], 
				$input_reads, $opt->{stranded} ? " --SS_lib_type $opt->{strand}" : ''
			);
			$opt->{cmd} .= sprintf("Trinity --no_cleanup --seqType fq --max_memory 8G --bflyHeapSpaceMax 4G --CPU 1 --grid_exec %s/%s --output %s --no_normalize_reads %s%s; ",
				$opt->{binpath}, $opt->{env}->{trinity_grid_module}, $opt->{dir_list}->[0], 
				$input_reads, $opt->{stranded} ? " --SS_lib_type $opt->{strand}" : ''
			);
		}
		my $metrics_file = $report_db_folder."/".$assembly_step->get_or_create_metrics_filename('fastaAssembly');
		$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[0]/Trinity.fasta > $opt->{dir_list}->[0]/assembler_metrics.tsv ; $opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[0]/assembler_metrics.tsv > $metrics_file ; \\rm -f $opt->{dir_list}->[0]/assembler_metrics.tsv\n";
		# no need to clean a-trinity to let Trinity restarting at the right step
		$opt->{complete} = 0;
	}
}
if (exists $opt->{cmd}) {
	if ($opt->{restart}) {
		$clean_msg .= clean_directories($opt, 1);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 2);

# merge kmers assemblies and run seqclean if oases ; remove contigs with Ns (fasta_noN) ; detect self chimeras
# write file 03.sh
my $merge_analysis = undef;
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
unless (-f "$opt->{dir_list}->[1]/transcripts.fa" && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		$clean_msg .= clean_directories($opt, 1);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	if ($opt->{dbg} eq 'oases') {
		my $nkmers = scalar(@{$opt->{kmers}});
		my @transcripts = glob("$opt->{dir_list}->[0]/k*/transcripts.1\%Locus.20\%Longest.fa");
		$opt->{cmd} = "if (`find $opt->{dir_list}->[0] -name transcripts.fa | wc -l` != $nkmers) exit\n";
		$opt->{cmd} .= "set pwd = `pwd`\n";
		unless (scalar(@transcripts) == scalar(@{$opt->{kmers}})) {
			$opt->{cmd} .= "$opt->{binpath}/oases_best_trans_chooser.sh $opt->{dir_list}->[0]\n";
		}
		$opt->{cmd} .= "find $opt->{dir_list}->[1] -name transcripts.raw.fa -exec \\rm -f {} \\;\n"
		              ."foreach f ($opt->{dir_list}->[0]/k*/transcripts.best.80\%Longest.fa)\n"
		              ."\tmkdir -p `dirname \$f`/seqclean\n\tcd `dirname \$f`/seqclean\n\tseqclean ../transcripts.best.80\%Longest.fa -o ../transcripts.best.80\%Longest.fa.cln\n";
		$opt->{cmd} .= sprintf("\tcd \$pwd\n\tset g = `dirname \$f | xargs basename`\n\tcat \$f.cln%s | sed \"s/>/>\${g}_/;s|/.*||\" >> %s/transcripts.raw.fa\nend\n",
			$opt->{discardN} ? " | $opt->{binpath}/fasta_no_N.pl -ex" : '', $opt->{dir_list}->[1]
		);
		# Report
		$merge_analysis = $report->get_or_create_analysis("Merge and clean", "Merge k-mers, clean sequence and correct self-chimera.");
		my $locusFilter_step = $merge_analysis->get_or_create_step('Locus filter', 'Filter contigs to produce one contig by locus.');
		my $seqClean_step = $merge_analysis->get_or_create_step('Seq cleaning', 'Trim ends rich in Ns, discard N-rich sequences, trim polyA tails from 3\' end or polyT from 5\' end of sequences and discard low-complexity sequences.');
		foreach my $kmer (@{$opt->{kmers}}) {
			my $locusFilter_metrics_file = $report_db_folder."/".$locusFilter_step->get_or_create_metrics_filename('fastaAssembly', 'k-mer_'.$kmer);
			$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[0]/k$kmer/transcripts.best.80\%Longest.fa > $opt->{dir_list}->[1]/locusFilter_metrics_k$kmer.tsv\n"
			              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[1]/locusFilter_metrics_k$kmer.tsv > $locusFilter_metrics_file\n"
			              ."\\rm -f $opt->{dir_list}->[1]/locusFilter_metrics_k$kmer.tsv\n";
			my $seqClean_metrics_file = $report_db_folder."/".$seqClean_step->get_or_create_metrics_filename('fastaAssembly', 'k-mer_'.$kmer);
			$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[0]/k$kmer/transcripts.best.80\%Longest.fa.cln > $opt->{dir_list}->[1]/seqclean_metrics_k$kmer.tsv\n"
			              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[1]/seqclean_metrics_k$kmer.tsv > $seqClean_metrics_file\n"
			              ."\\rm -f $opt->{dir_list}->[1]/seqclean_metrics_k$kmer.tsv\n";
		}
		my $merge_metrics_file = $report_db_folder."/".$merge_analysis->get_or_create_step('Merge', 'Merge all contigs generated by the different k-mers.')->get_or_create_metrics_filename('fastaAssembly');
		$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[1]/transcripts.raw.fa > $opt->{dir_list}->[1]/merge_metrics.tsv\n"
		              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[1]/merge_metrics.tsv > $merge_metrics_file\n"
		              ."\\rm -f $opt->{dir_list}->[1]/merge_metrics.tsv\n";
	} else {
		$opt->{cmd} = "if (! -e $opt->{dir_list}->[0]/Trinity.fasta || -z $opt->{dir_list}->[0]/Trinity.fasta) exit\n";
		$merge_analysis = $report->get_or_create_analysis('Correct chimera', 'Correct self-chimera.');
		$opt->{cmd} .= sprintf("cat $opt->{dir_list}->[0]/Trinity.fasta%s | sed 's/\\(>[^ ]*\\) .*/\\1/' > %s/transcripts.raw.fa\n", $opt->{discardN} ? " | $opt->{binpath}/fasta_no_N.pl -ex" : '', $opt->{dir_list}->[1]);
	}
	$opt->{cmd} .= sprintf("%s/submitSelfChimera.pl --max-sub 70 -i %s/transcripts.raw.fa -o %s/transcripts.fa -l %s/transcripts.fa.chimFilter.log%s\n",
		$opt->{binpath}, $opt->{dir_list}->[1], $opt->{dir_list}->[1], $opt->{dir_list}->[1], $opt->{local} ? ' --local' : ''
	);
	my $chimera_step = $merge_analysis->get_or_create_step('Correct chimera', 'Remove repeat(s) of self-chimera.');
	my $chimeraFasta_metrics_file = $report_db_folder."/".$chimera_step->get_or_create_metrics_filename('fastaAssembly');
	$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[1]/transcripts.fa > $opt->{dir_list}->[1]/chimera_fasta_metrics.tsv\n"
	              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[1]/chimera_fasta_metrics.tsv > $chimeraFasta_metrics_file\n"
	              ."\\rm -f $opt->{dir_list}->[1]/chimera_fasta_metrics.tsv\n";
	my $chimeraLog_metrics_file = $report_db_folder."/".$chimera_step->get_or_create_metrics_filename('chimeraLog');
	$opt->{cmd} .= $opt->{binpath}."/chimeraMetrics2json.pl ".$opt->{dir_list}->[1]."/transcripts.fa.chimFilter.log > ".$chimeraLog_metrics_file."\n";
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 3);

# cluster DBG contigs using cd-hit
# write file 04.sh
my $cluster_analysis = $report->get_or_create_analysis('Remove inclusion', 'Remove contigs included in other(s).');
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
unless (-f "$opt->{dir_list}->[2]/all_dbg.fa" && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		$clean_msg .= clean_directories($opt, 2);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	$opt->{cmd} = "if (! -e $opt->{dir_list}->[1]/transcripts.fa || -z $opt->{dir_list}->[1]/transcripts.fa) exit\n";
	$opt->{cmd} .= "cd-hit-est -i $opt->{dir_list}->[1]/transcripts.fa -o $opt->{dir_list}->[2]/all_dbg.fa -M 0 -d 0 -c 0.98 -T $opt->{env}->{n_cpu} > $opt->{dir_list}->[2]/all_dbg.fa.cd-hit.log\n"
	              ."cat $opt->{dir_list}->[2]/all_dbg.fa.clstr | ".q(perl -le '$/="\n>";map{s/^\d+.+>(\S+)\.\.\. \*$//m;$centroid=$1;print "$centroid\t$centroid";while(s/^\d+.+>(\S+)\.\.\..+//m){print "$1\t$centroid"}}<STDIN>')." > $opt->{dir_list}->[2]/all_dbg.fa.cd-hit.history.log\n";
	my $clustering_metrics_file = $report_db_folder."/".$cluster_analysis->get_or_create_step('cd-hit clustering')->get_or_create_metrics_filename('fastaAssembly');
	$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[2]/all_dbg.fa > $opt->{dir_list}->[2]/cd-hit_metrics.tsv\n"
	              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[2]/cd-hit_metrics.tsv > $clustering_metrics_file\n"
	              ."\\rm -f $opt->{dir_list}->[2]/cd-hit_metrics.tsv\n";
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 4);

# assemble DBG contigs using OLC assembler
# write file 05.sh
my $largeAsm_analysis = $report->get_or_create_analysis('Large scale assembly', 'Assemble contigs with Overlap Layout Consensus assembler.');
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
unless (-f "$opt->{dir_list}->[2]/all_contigs.raw.fa" && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		$clean_msg .= process_cmd(1, "$opt->{binpath}/tgicl_clean_dirs.sh $opt->{dir_list}->[2]");
		$clean_msg .= clean_directories($opt, 3);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	$opt->{cmd} = "if (! -f $opt->{dir_list}->[2]/all_dbg.fa || -z $opt->{dir_list}->[2]/all_dbg.fa) exit\n"
	             ."set pwd = `pwd`\n"
	             ."cd $opt->{dir_list}->[2]\n./runAssembly.sh\n"
	             ."cd \$pwd\n";
	open(SRC, "$opt->{binpath}/runAssembly.sh") or croak "Can't open file $opt->{binpath}/runAssembly.sh";
	open(SH, ">$opt->{dir_list}->[2]/runAssembly.sh") or croak "Can't open file $opt->{dir_list}->[2]/runAssembly.sh";
	while (<SRC>) {
		s|BINPATH|$opt->{binpath}|;
		s|CFGPATH|$opt->{cfgpath}|;
		s|TGICL_CPU|$opt->{env}->{n_cpu}|;
		print SH;
	}
	close SH;
	close SRC;
	copy("$opt->{cfgpath}/tgicl.cfg", "$opt->{dir_list}->[2]");
	chmod(0775, "$opt->{dir_list}->[2]/runAssembly.sh");
	my $largeAsm_metrics_file = $report_db_folder."/".$largeAsm_analysis->get_or_create_step('tgicl assembly')->get_or_create_metrics_filename('fastaAssembly');
	$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[2]/all_contigs.raw.fa > $opt->{dir_list}->[2]/tgicl_metrics.tsv\n"
	              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[2]/tgicl_metrics.tsv > $largeAsm_metrics_file\n"
	              ."\\rm -f $opt->{dir_list}->[2]/tgicl_metrics.tsv\n";
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 5);

# run Transdecoder and cut fusion contigs if filter type is orf ; run VecScreen filter to become TSA compliant
# write file 06.sh
my $clean_analysis = $report->get_or_create_analysis('Inner cleaning', 'Split contigs with multiple ORFs, remove vector sequences and filter by length.');
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
unless (-f "$opt->{dir_list}->[3]/all_contigs.fa" && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		$clean_msg .= clean_directories($opt, 3);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	if ($opt->{type} eq 'orf') {
		my $base = basename($opt->{dir_list}->[2]);
		$opt->{cmd} = "if (! -e $opt->{dir_list}->[2]/all_contigs.raw.fa || -z $opt->{dir_list}->[2]/all_contigs.raw.fa) exit\n"
		             ."set pwd = `pwd`\n"
		             .sprintf("cd %s\n%s -t ../%s/all_contigs.raw.fa%s\n", $opt->{dir_list}->[3], "TransDecoder.LongOrfs", $base, $opt->{stranded} ? ' -S' : '')
		             .sprintf("%s -t ../%s/all_contigs.raw.fa\n", "TransDecoder.Predict", $base)
		             ."cat ../$base/all_contigs.raw.fa | $opt->{binpath}/transdecoder_bed_to_false_fusion_cutter.pl -f stdin -b all_contigs.raw.fa.transdecoder.bed --log all_contigs.raw.fa.transdecoder_cutter.log --orf-length all_contigs.raw.fa.transdecoder_cutter.orf_length.tsv --keep > all_contigs.raw.fa.transdecoder_cutter.fa\n"
		             .'sed 1d all_contigs.raw.fa.transdecoder_cutter.log | perl -lane \'$c=$F[2]=~s/,//g;map{print"$F[0]\t$F[0]_".++$_}0..$c\' > all_contigs.raw.fa.transdecoder_cutter.history.log'."\n"
		             ."cd \$pwd\n";
		my $splitORF_step = $clean_analysis->get_or_create_step('ORF splitting', 'Find multiple ORF contigs and split them.');
		my $splitORFFasta_metrics_file = $report_db_folder."/".$splitORF_step->get_or_create_metrics_filename('fastaAssembly');
		$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.fa > $opt->{dir_list}->[3]/transdecoder_fasta_metrics.tsv\n"
		              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[3]/transdecoder_fasta_metrics.tsv > $splitORFFasta_metrics_file\n"
		              ."\\rm -f $opt->{dir_list}->[3]/transdecoder_fasta_metrics.tsv\n";
		my $splitORFLog_metrics_file = $report_db_folder."/".$splitORF_step->get_or_create_metrics_filename('splitORFLog');
		$opt->{cmd} .= "$opt->{binpath}/transdecoderMetrics2json.pl $opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.log > $splitORFLog_metrics_file\n";
	}
	$opt->{cmd} .= sprintf("%s/submitFastaTsa.pl --max-sub 70 -i %s -o %s/all_contigs.fa -l %s/all_contigs.fa.vecFilter.log --length %d%s\n",
		$opt->{binpath}, $opt->{type} eq 'orf' ? "$opt->{dir_list}->[3]/all_contigs.raw.fa.transdecoder_cutter.fa" : "$opt->{dir_list}->[2]/all_contigs.raw.fa",
		$opt->{dir_list}->[3], $opt->{dir_list}->[3], $opt->{length}, $opt->{local} ? ' --local' : ''
	);
	$opt->{cmd} .= sprintf(q(grep '^>' %s/all_contigs.fa | tr -d '>' | awk '{print $0"\\t"$0}' > %s/all_contigs.fa.vecFilter.history.log%ssed 1d %s/all_contigs.fa.vecFilter.log | grep '^#' | cut -f1 | tr -d '#' | awk '{print $0"\\t"}' >> %s/all_contigs.fa.vecFilter.history.log%ssort -o %s/all_contigs.fa.vecFilter.history.log %s/all_contigs.fa.vecFilter.history.log%s),
		$opt->{dir_list}->[3], $opt->{dir_list}->[3], "\n", $opt->{dir_list}->[3], $opt->{dir_list}->[3], "\n", $opt->{dir_list}->[3], $opt->{dir_list}->[3], "\n"
	);
	my $vectorFilter_metrics_file = $report_db_folder."/".$clean_analysis->get_or_create_step('Vector and length filter', "Cut vector sequences and filter contigs on length (>= $opt->{length} nucleotides).")->get_or_create_metrics_filename('fastaAssembly');
	$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[3]/all_contigs.fa > $opt->{dir_list}->[3]/vecscreen_metrics.tsv\n"
	              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[3]/vecscreen_metrics.tsv > $vectorFilter_metrics_file\n"
	              ."\\rm -f $opt->{dir_list}->[3]/vecscreen_metrics.tsv\n";
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 6);

# map reads back to contigs using BWA or STAR, correct variations and run again cd-hit-est
# write file 07.sh
my $aln_edit_analysis = $report->get_or_create_analysis('Alignment-Editing', 'Align reads back to contigs and edit contigs (fix insertion/deletion/substitution bugs).');
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
$opt->{rmbt_ref}->{editing_first} = $opt->{mapper} eq 'star' ? "$opt->{dir_list}->[4]/STAR_all_contigs.fa" : "$opt->{dir_list}->[3]/all_contigs.fa";
$opt->{rmbt_ref}->{editing_second} = $opt->{mapper} eq 'star' ? "$opt->{dir_list}->[4]/STAR_all_contigs.first_pass.fa" : "$opt->{dir_list}->[4]/all_contigs.first_pass.fa";
unless (-f "$opt->{dir_list}->[4]/all_contigs.second_pass.fa" && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		unlink("$opt->{dir_list}->[3]/all_contigs.fa.fai");
		unlink("$opt->{dir_list}->[4]/all_contigs.first_pass.fa.fai");
		if (! -f "$opt->{dir_list}->[4]/first_pass.over") {
			# first pass not complete - clean following rmbt process
			foreach my $rmbt_step ('editing_second', 'filtering') {
				foreach my $alignR1 (keys %{$opt->{rmbt_status}->{$opt->{rmbt_ref}->{$rmbt_step}}}) {
					$opt->{rmbt_status}->{$opt->{rmbt_ref}->{$rmbt_step}}->{$alignR1}->{status} = 0;
				}
				$clean_msg .= clean_rmbt_reference($opt, $opt->{dir_list}->[4], $opt->{rmbt_ref}->{$rmbt_step});
			}
		} elsif (! -f "$opt->{dir_list}->[4]/second_pass.over") {
			# second editing pass not complete - clean following rmbt process
			foreach my $alignR1 (keys %{$opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}}) {
				$opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}->{$alignR1}->{status} = 0;
			}
			$clean_msg .= clean_rmbt_reference($opt, $opt->{dir_list}->[4], $opt->{rmbt_ref}->{filtering});
		}
		$clean_msg .= clean_rmbt_directory($opt, $opt->{dir_list}->[4]);
		$clean_msg .= clean_directories($opt, 5);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	} else {
		# First pass
		# runSamse only for indexing reference
		$opt->{cmd} = "if (! -e $opt->{dir_list}->[3]/all_contigs.fa || -z $opt->{dir_list}->[3]/all_contigs.fa) exit\n";
		$opt->{cmd} .= sprintf ("%s/runSamse.sh -r %s/all_contigs.fa -f %s -o %s -m %s -t %d --force --index --sync%s\n",
			$opt->{binpath}, $opt->{dir_list}->[3], $opt->{alignR1}->[0], $opt->{dir_list}->[4],
			$opt->{mapper}, $opt->{env}->{n_cpu}, $opt->{local} ? ' --local' : ''
		);
	}
	unless (-f "$opt->{dir_list}->[4]/first_pass.over") {
		# runSamse or runSampe to map back reads
		$opt->{cmd} .= sprintf("set reference = %s\n", $opt->{mapper} eq 'star' ? "`find $opt->{dir_list}->[4] -name STAR_all_contigs.fa_\\*`" : "$opt->{dir_list}->[3]/all_contigs.fa");
		my $rmbt_submit = 0;
		for (my $i = 0 ; $i < scalar(@{$opt->{alignR1}}) ; $i++) {
			next if (exists $opt->{rmbt_status}->{$opt->{rmbt_ref}->{editing_first}}->{$opt->{alignR1}->[$i]} && $opt->{rmbt_status}->{$opt->{rmbt_ref}->{editing_first}}->{$opt->{alignR1}->[$i]}->{status} == 1);
			my $fastq = $opt->{paired} ? basename($opt->{alignR1}->[$i]).'/'.basename($opt->{alignR2}->[$i]) : basename($opt->{alignR1}->[$i]);
			$opt->{cmd} .= sprintf("set bam = `%s/%s -r \$reference%s %s -o %s -m %s -t %d --flagstat --bam --sort_by_pos%s%s`\n",
				$opt->{binpath}, $opt->{paired} ? 'runSampe.sh' : 'runSamse.sh', $opt->{mapper} eq 'star' ? ' -g' : '',
				$opt->{paired} ? "-1 $opt->{alignR1}->[$i] -2 $opt->{alignR2}->[$i]" : "-f $opt->{alignR1}->[$i]",
				$opt->{dir_list}->[4], $opt->{mapper}, $opt->{env}->{n_cpu}, $opt->{filter} ? ' -y' : '', $opt->{local} ? ' --local' : ''
			);
			my $flagstat_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step($opt->{mapper}." mapping - first pass")->get_or_create_metrics_filename('flagstatLog', $fastq);
			my $flagstat_metrics_basename = basename($flagstat_metrics_file, '.json');
			if ($opt->{local}) {
				$opt->{cmd} .= "$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file\n";
			} else {
				$opt->{cmd} .= "set bam_jid = `echo \$bam | sed -e 's/.*\\.\\([0-9]*\\)\\.bam/\\1/'`\n"
				              ."set hold_jid = `cat $opt->{dir_list}->[4]/err_log_\$bam_jid/runSam?e.\$bam_jid.jid | tr '\\n' ,`\n"
				              ."set metric = \"$opt->{binpath}/submitJob --name first_pass_$flagstat_metrics_basename --stdout $dir/err_log --stderr $dir/err_log --options -hold_jid \$hold_jid --binary -- '$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file'\"\n"
				              ."echo \$metric\nset metric_out = `eval \$metric`\necho \$metric_out\n"
			}
			$rmbt_submit++;
		}
		unless ($opt->{local} || !$rmbt_submit) {
			$opt->{cmd} .= "set all_jid = `cat $opt->{dir_list}->[4]/err_log_*/runSam?e.*.jid | tr '\\n' ,`\n"
			              ."set pending = \"$opt->{binpath}/submitJob --name pending_first_pass_rmbt --sync --stdout $dir/err_log --stderr $dir/err_log --options -hold_jid \$all_jid --binary -- echo first pass rmbt ended\"\n"
			              ."echo \$pending\nset pending_out = `eval \$pending`\necho \$pending_out\n"
		}
		$opt->{cmd} .= "samtools faidx $opt->{dir_list}->[3]/all_contigs.fa\n"
		              ."set all_bams = `find $opt->{dir_list}->[4] -name \\*all_contigs.\\[0-9\\]\\*.bam`\nif (\"\$all_bams\" == '') exit\n";
		$opt->{cmd} .= sprintf("%s/submitSamCorrectVar.pl --max-sub 70 --fasta %s/all_contigs.fa --bam \$all_bams --log %s/all_contigs.first_pass.raw.fa.samCorrectVariation.log --output %s/all_contigs.first_pass.raw.fa%s\n",
			$opt->{binpath}, $opt->{dir_list}->[3], $opt->{dir_list}->[4], $opt->{dir_list}->[4], $opt->{local} ? ' --local' : ''
		);
		$opt->{cmd} .= "cd-hit-est -i $opt->{dir_list}->[4]/all_contigs.first_pass.raw.fa -o $opt->{dir_list}->[4]/all_contigs.first_pass.fa -M 0 -d 0 -c 0.98 -T $opt->{env}->{n_cpu} > $opt->{dir_list}->[4]/all_contigs.first_pass.fa.cd-hit.log\n"
		              ."cat $opt->{dir_list}->[4]/all_contigs.first_pass.fa.clstr | ".q(perl -le '$/="\n>";map{s/^\d+.+>(\S+)\.\.\. \*$//m;$centroid=$1;print "$centroid\t$centroid";while(s/^\d+.+>(\S+)\.\.\..+//m){print "$1\t$centroid"}}<STDIN>')." > $opt->{dir_list}->[4]/all_contigs.first_pass.fa.cd-hit.history.log\n";
		my $correctRef_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step('Correct consensus - first pass', 'Correct insertions, deletions and substitutions on contigs. It follow the majority vote at each position of the alignment.')->get_or_create_metrics_filename('correctVariationLog');
		$opt->{cmd} .= "$opt->{binpath}/samCorrectVarMetrics2json.pl $opt->{dir_list}->[4]/all_contigs.first_pass.raw.fa.samCorrectVariation.log > $correctRef_metrics_file\n";
		my $clustering_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step('cd-hit clustering - first pass')->get_or_create_metrics_filename('fastaAssembly');
		$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[4]/all_contigs.first_pass.fa > $opt->{dir_list}->[4]/cd-hit_first_pass_metrics.tsv\n"
		              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[4]/cd-hit_first_pass_metrics.tsv > $clustering_metrics_file\n"
		              ."\\rm -f $opt->{dir_list}->[4]/cd-hit_first_pass_metrics.tsv\n";
	# Second pass
	# runSamse only for indexing reference
		$opt->{cmd} .= sprintf ("%s/runSamse.sh -r %s/all_contigs.first_pass.fa -f %s -o %s -m %s -t %d --force --index --sync%s\n",
			$opt->{binpath}, $opt->{dir_list}->[4], $opt->{alignR1}->[0], $opt->{dir_list}->[4],
			$opt->{mapper}, $opt->{env}->{n_cpu}, $opt->{local} ? ' --local' : ''
		);
	}	
	# runSamse or runSampe to map back reads
	$opt->{cmd} .= sprintf("set reference = %s\n", $opt->{mapper} eq 'star' ? "`find $opt->{dir_list}->[4] -name STAR_all_contigs.first_pass.fa_\\*`" : "$opt->{dir_list}->[4]/all_contigs.first_pass.fa");
	my $rmbt_submit = 0;
	for (my $i = 0 ; $i < scalar(@{$opt->{alignR1}}) ; $i++) {
		next if (exists $opt->{rmbt_status}->{$opt->{rmbt_ref}->{editing_second}}->{$opt->{alignR1}->[$i]} && $opt->{rmbt_status}->{$opt->{rmbt_ref}->{editing_second}}->{$opt->{alignR1}->[$i]}->{status} == 1);
		my $fastq = $opt->{paired} ? basename($opt->{alignR1}->[$i]).'/'.basename($opt->{alignR2}->[$i]) : basename($opt->{alignR1}->[$i]);
		$opt->{cmd} .= sprintf("set bam = `%s/%s -r \$reference%s %s -o %s -m %s -t %d --flagstat --bam --sort_by_pos%s%s`\n",
			$opt->{binpath}, $opt->{paired} ? 'runSampe.sh' : 'runSamse.sh', $opt->{mapper} eq 'star' ? ' -g' : '',
			$opt->{paired} ? "-1 $opt->{alignR1}->[$i] -2 $opt->{alignR2}->[$i]" : "-f $opt->{alignR1}->[$i]",
			$opt->{dir_list}->[4], $opt->{mapper}, $opt->{env}->{n_cpu}, $opt->{filter} ? ' -y' : '', $opt->{local} ? ' --local' : ''
		);
		my $flagstat_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step($opt->{mapper}." mapping - second pass")->get_or_create_metrics_filename('flagstatLog', $fastq);
		my $flagstat_metrics_basename = basename($flagstat_metrics_file, '.json');
		if ($opt->{local}) {
			$opt->{cmd} .= "$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file\n";
		} else {
			$opt->{cmd} .= "set bam_jid = `echo \$bam | sed -e 's/.*\\.\\([0-9]*\\)\\.bam/\\1/'`\n"
			              ."set hold_jid = `cat $opt->{dir_list}->[4]/err_log_\$bam_jid/runSam?e.\$bam_jid.jid | tr '\\n' ,`\n"
			              ."set metric = \"$opt->{binpath}/submitJob --name second_pass_$flagstat_metrics_basename --stdout $dir/err_log --stderr $dir/err_log --options -hold_jid \$hold_jid --binary -- '$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file'\"\n"
			              ."echo \$metric\nset metric_out = `eval \$metric`\necho \$metric_out\n"
		}
		$rmbt_submit++;
	}
	unless ($opt->{local} || !$rmbt_submit) {
		$opt->{cmd} .= "set all_jid = `cat $opt->{dir_list}->[4]/err_log_*/runSam?e.*.jid | tr '\\n' ,`\n"
		              ."set pending = \"$opt->{binpath}/submitJob --name pending_second_pass_rmbt --sync --stdout $dir/err_log --stderr $dir/err_log --options -hold_jid \$all_jid --binary -- echo second pass rmbt ended\"\n"
		              ."echo \$pending\nset pending_out = `eval \$pending`\necho \$pending_out\n"
	}
	$opt->{cmd} .= "samtools faidx $opt->{dir_list}->[4]/all_contigs.first_pass.fa\n"
	              ."set all_bams = `find $opt->{dir_list}->[4] -name \\*all_contigs.first_pass.\\[0-9\\]\\*.bam`\nif (\"\$all_bams\" == '') exit\n";
	$opt->{cmd} .= sprintf("%s/submitSamCorrectVar.pl --max-sub 70 --fasta %s/all_contigs.first_pass.fa --bam \$all_bams --log %s/all_contigs.second_pass.raw.fa.samCorrectVariation.log --output %s/all_contigs.second_pass.raw.fa%s\n",
		$opt->{binpath}, $opt->{dir_list}->[4], $opt->{dir_list}->[4], $opt->{dir_list}->[4], $opt->{local} ? ' --local' : ''
	);
	$opt->{cmd} .= "cd-hit-est -i $opt->{dir_list}->[4]/all_contigs.second_pass.raw.fa -o $opt->{dir_list}->[4]/all_contigs.second_pass.fa -M 0 -d 0 -c 0.98 -T $opt->{env}->{n_cpu} > $opt->{dir_list}->[4]/all_contigs.second_pass.fa.cd-hit.log\n"
	              ."cat $opt->{dir_list}->[4]/all_contigs.second_pass.fa.clstr | ".q(perl -le '$/="\n>";map{s/^\d+.+>(\S+)\.\.\. \*$//m;$centroid=$1;print "$centroid\t$centroid";while(s/^\d+.+>(\S+)\.\.\..+//m){print "$1\t$centroid"}}<STDIN>')." > $opt->{dir_list}->[4]/all_contigs.second_pass.fa.cd-hit.history.log\n";
	my $correctRef_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step('Correct consensus - second pass', 'Correct insertions, deletions and substitutions on contigs. It follow the majority vote at each position of the alignment.')->get_or_create_metrics_filename('correctVariationLog');
	$opt->{cmd} .= "$opt->{binpath}/samCorrectVarMetrics2json.pl $opt->{dir_list}->[4]/all_contigs.second_pass.raw.fa.samCorrectVariation.log > $correctRef_metrics_file\n";
	my $clustering2_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step('cd-hit clustering - second pass')->get_or_create_metrics_filename('fastaAssembly');
	$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[4]/all_contigs.second_pass.fa > $opt->{dir_list}->[4]/cd-hit_second_pass_metrics.tsv\n"
	              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[4]/cd-hit_second_pass_metrics.tsv > $clustering2_metrics_file\n"
	              ."\\rm -f $opt->{dir_list}->[4]/cd-hit_second_pass_metrics.tsv\n";
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 
quit($opt) if ($opt->{quit} == 7);

# map reads back to contigs using BWA or STAR, run eXpress and filter contigs on length and coverage criterias
# write file 08.sh
my $aln_filter_analysis = $report->get_or_create_analysis('Alignment-Filtering', 'Align reads back to contigs and filter low coverage contigs.');
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
my $lowest_fpkm = $opt->{coverages}->[0];
$opt->{rmbt_ref}->{filtering} = $opt->{mapper} eq 'star' ? "$opt->{dir_list}->[5]/STAR_all_contigs.second_pass.fa" : "$opt->{dir_list}->[4]/all_contigs.second_pass.fa";
unless (-f "$opt->{dir_list}->[5]/transcripts_fpkm_$lowest_fpkm.fa" && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		$clean_msg .= clean_rmbt_directory($opt, $opt->{dir_list}->[5]);
		$clean_msg .= clean_directories($opt, 6);
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	} else {
	# runSamse only for indexing reference
		$opt->{cmd} = "if (! -e $opt->{dir_list}->[4]/all_contigs.second_pass.fa || -z $opt->{dir_list}->[4]/all_contigs.second_pass.fa) exit\n";
		$opt->{cmd} .= sprintf ("%s/runSamse.sh -r %s/all_contigs.second_pass.fa -f %s -o %s -m %s -t %d --force --index --sync%s\n",
			$opt->{binpath}, $opt->{dir_list}->[4], $opt->{alignR1}->[0], $opt->{dir_list}->[5],
			$opt->{mapper}, $opt->{env}->{n_cpu}, $opt->{local} ? ' --local' : ''
		);
	}
	# runSamse or runSampe to map back reads
	$opt->{cmd} .= sprintf("set reference = %s\n", $opt->{mapper} eq 'star' ? "`find $opt->{dir_list}->[5] -name STAR_all_contigs.second_pass.fa_\\*`" : "$opt->{dir_list}->[4]/all_contigs.second_pass.fa");
	my $rmbt_submit = 0;
	for (my $i = 0 ; $i < scalar(@{$opt->{alignR1}}) ; $i++) {
		next if (exists $opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}->{$opt->{alignR1}->[$i]} && $opt->{rmbt_status}->{$opt->{rmbt_ref}->{filtering}}->{$opt->{alignR1}->[$i]}->{status} == 1);
		my $fastq = $opt->{paired} ? basename($opt->{alignR1}->[$i]).'/'.basename($opt->{alignR2}->[$i]) : basename($opt->{alignR1}->[$i]);
		$opt->{cmd} .= sprintf("set bam = `%s/%s -r \$reference%s %s -o %s -m %s -t %d --flagstat --bam --sort_by_name%s%s`\n",
			$opt->{binpath}, $opt->{paired} ? 'runSampe.sh' : 'runSamse.sh', $opt->{mapper} eq 'star' ? ' -g' : '',
			$opt->{paired} ? "-1 $opt->{alignR1}->[$i] -2 $opt->{alignR2}->[$i]" : "-f $opt->{alignR1}->[$i]",
			$opt->{dir_list}->[5], $opt->{mapper}, $opt->{env}->{n_cpu}, $opt->{filter} ? ' -y' : '', $opt->{local} ? ' --local' : ''
		);
		my $flagstat_metrics_file = $report_db_folder."/".$aln_edit_analysis->get_or_create_step($opt->{mapper}." mapping")->get_or_create_metrics_filename('flagstatLog', $fastq);
		my $flagstat_metrics_basename = basename($flagstat_metrics_file, '.json');
		if ($opt->{local}) {
			$opt->{cmd} .= "$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file\n";
		} else {
			$opt->{cmd} .= "set bam_jid = `echo \$bam | sed -e 's/.*\\.\\([0-9]*\\)\\.bam/\\1/'`\n"
			              ."set hold_jid = `cat $opt->{dir_list}->[5]/err_log_\$bam_jid/runSam?e.\$bam_jid.jid | tr '\\n' ,`\n"
			              ."set metric = \"$opt->{binpath}/submitJob --name filtering_$flagstat_metrics_basename --stdout $dir/err_log --stderr $dir/err_log --options -hold_jid \$hold_jid --binary -- '$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file'\"\n"
			              ."echo \$metric\nset metric_out = `eval \$metric`\necho \$metric_out\n"
		}
		$rmbt_submit++;
	}
	unless ($opt->{local} || !$rmbt_submit) {
		$opt->{cmd} .= "set all_jid = `cat $opt->{dir_list}->[5]/err_log_*/runSam?e.*.jid | tr '\\n' ,`\n"
		              ."set pending = \"$opt->{binpath}/submitJob --name pending_filtering_rmbt --sync --stdout $dir/err_log --stderr $dir/err_log --options -hold_jid \$all_jid --binary -- echo filtering rmbt ended\"\n"
		              ."echo \$pending\nset pending_out = `eval \$pending`\necho \$pending_out\n"
	}
	$opt->{cmd} .= "set all_bams = `find $opt->{dir_list}->[5] -name \\*all_contigs.second_pass.\\[0-9\\]\\*.bam | tr '\\n' ,`\nif (\$all_bams == '') exit\n";
	$opt->{cmd} .= sprintf("express --no-update-check --no-bias-correct --logtostderr --output-dir %s%s %s/all_contigs.second_pass.fa \$all_bams\n",
		$opt->{dir_list}->[5], $opt->{stranded} ? ' --'.lc($opt->{strand}).'-stranded' : '', $opt->{dir_list}->[4]
	);
	$opt->{cmd} .= sprintf("%s/cov_length_filter.sh -f %s/all_contigs.second_pass.fa -o %s -x %s/results.xprs -b %s/all_contigs.raw.fa.transdecoder_cutter.orf_length.tsv -t %s -l %d -c %s\n",
		$opt->{binpath}, $opt->{dir_list}->[4], $opt->{dir_list}->[5], $opt->{dir_list}->[5], $opt->{dir_list}->[3], $opt->{type}, $opt->{length}, $opt->{optimize} ? $lowest_fpkm : $opt->{coverage}
	);
	my $flagstat_metrics_file = $report_db_folder."/".$aln_filter_analysis->get_or_create_step($opt->{mapper}." mapping")->get_or_create_metrics_filename('flagstatLog');
	$opt->{cmd} .= "$opt->{binpath}/alignmentMetrics2json.pl \$bam.flagstat > $flagstat_metrics_file\n";
	my $step_description = sprintf("Filter on coverage (fpkm: %s) and length (%s).", 
		$opt->{coverage}, $opt->{type} eq 'orf' ? "contigs >= $opt->{length} nucleotides" : "contigs with putative orf >= $opt->{length} nucleotides"
	);
	my $coverageFilter_step = $aln_filter_analysis->get_or_create_step('Coverage/length filter', $step_description);
	my $expressLog_metrics_file = $report_db_folder."/".$coverageFilter_step->get_or_create_metrics_filename('expressLog');
	$opt->{cmd} .= "$opt->{binpath}/expressMetrics2json.pl $opt->{dir_list}->[5]/results.xprs > $expressLog_metrics_file\n";
	foreach my $fpkm ( @{$opt->{coverages}} ){
		last if ($opt->{optimize} && $fpkm > $lowest_fpkm);
		my $expressFasta_metrics_file = $report_db_folder."/".$coverageFilter_step->get_or_create_metrics_filename('fastaAssembly', 'fpkm_'.$fpkm );
		$opt->{cmd} .= "$opt->{binpath}/fastaAssemblyMetrics.pl --input $opt->{dir_list}->[5]/transcripts_fpkm_$fpkm.fa > $opt->{dir_list}->[5]/fpkm-${fpkm}_metrics.tsv\n"
		              ."$opt->{binpath}/fastaMetrics2json.pl $opt->{dir_list}->[5]/fpkm-${fpkm}_metrics.tsv > $expressFasta_metrics_file\n"
		              ."\\rm -f $opt->{dir_list}->[5]/fpkm-${fpkm}_metrics.tsv\n";
	}
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 

# check each step and clean tree if all steps finished without errors
# write file 09.sh
my $score_analysis = undef;
if( !$opt->{no_rate} && $opt->{paired} ){
	$score_analysis = $report->get_or_create_analysis('Assembly scoring', 'Measures the quality of the assembly. A score is produced for the whole assembly, and for each contig.');
}
$opt->{step} = $opt->{steps}->[++$nstep];
$opt->{complete} = 1;
unless ((($opt->{no_rate} && -f "$dir/00-ASSEMBLY_COMPLETE") || -f "$dir/00-ASSEMBLY_RATING") && step_complete($opt)) {
	$opt->{complete} = 0;
	if ($opt->{restart}) {
		$opt->{restart} = 0;
		$resubmit = $opt->{step};
	}
	my $tab = "";
	if (-f "$dir/00-ASSEMBLY_COMPLETE") {
		$resubmit = 'transrate';
	} else {
		$opt->{cmd} = "$opt->{binpath}/check_assembly.pl $dir\n"
	             ."if (! -e $dir/00-ASSEMBLY_COMPLETE) then\n"
	             ."\techo 'Assembly not complete; exit'\n"
	             ."\texit\n"
	             ."else\n"
	             ."\t$opt->{binpath}/$opt->{dbg}_clean_dirs.sh -no $opt->{dir_list}->[0]\n\t$opt->{binpath}/tgicl_clean_dirs.sh $opt->{dir_list}->[2]\n";
		$tab = "\t";
	}
	unless ($opt->{no_rate}) {
		$opt->{cmd} .= sprintf("%stransrate --assembly=%s/transcripts_fpkm_%d.fa --left=%s --right=%s --threads=%d --output=%s/transrate %s> %s/00-ASSEMBLY_RATING\n",
			$tab, $opt->{dir_list}->[5], $lowest_fpkm, join(',', @{$opt->{alignR1}}), join(',', @{$opt->{alignR2}}), $opt->{env}->{n_cpu},
			$dir, $opt->{ref} ? "--reference=$opt->{ref} " : '', $dir
		);
		$opt->{cmd} .= "${tab}$opt->{binpath}/check_assembly.pl $dir $opt->{step}\n";
		# Report
		my $scroring_metrics_file = $report_db_folder."/".$score_analysis->get_or_create_step("Transrate score")->get_or_create_metrics_filename('assemblyScoring');
		$opt->{cmd} .= sprintf( "%s%s/transrateMetrics2json.pl %s/00-ASSEMBLY_RATING > %s\n", 
			$tab, $opt->{binpath}, $dir, $scroring_metrics_file
		);
		# Filter
		if ($opt->{optimize}) {
			$opt->{cmd} .= "${tab}if (! -e transrate/transcripts_fpkm_$lowest_fpkm/good.transcripts_fpkm_$lowest_fpkm.fa) exit\n"
			              ."${tab}ln -fs transrate/transcripts_fpkm_$lowest_fpkm/good.transcripts_fpkm_$lowest_fpkm.fa $dir/transcripts_fpkm_$lowest_fpkm.fa\n"
			              ."${tab}set tmp = `mktemp`\n"
			              ."${tab}grep '^>' $opt->{dir_list}->[5]/coding_transcripts_fpkm_$lowest_fpkm.fa | tr -d '>' > \$tmp\n" 
			              ."${tab}cat $dir/transcripts_fpkm_$lowest_fpkm.fa | $opt->{binpath}/fasta_extract.pl \$tmp > $dir/coding_transcripts_fpkm_$lowest_fpkm.fa\n"
			              ."${tab}\\rm -f \$tmp\n";
		}
	}
	$opt->{cmd} .= "endif\n" unless (-f "$dir/00-ASSEMBLY_COMPLETE");
	write_shell($opt);
}
print_msg($opt) if ($opt->{display}); 

# align contigs to reference using exonerate and blat
# write file 10.sh
$opt->{step} = $opt->{steps}->[++$nstep];
if ($opt->{ref}) {
	my $ref_analysis = $report->get_or_create_analysis('Reference alignment', 'Align reference proteins to contigs.');
	my $ref_jobs;
	foreach my $tool (qw(exonerate blat)) {
		next unless ($opt->{$tool});
		$opt->{complete} = 1;
		$opt->{tool} = $tool;
		my ($dir, $param) = $tool eq 'exonerate' ? ($opt->{dir_list}->[6], " --params '-m protein2dna'") : ($opt->{dir_list}->[7], " --params '-q=prot -t=dnax'");
		unless (-d $dir && step_complete($opt)) {
			$opt->{complete} = 0;
			my $target = sprintf("%s/transcripts_fpkm_%d.fa", $opt->{dir_list}->[5], $lowest_fpkm);
			my $output = sprintf("%s/%s.%s.%s.%s", $dir, basename($opt->{ref}), basename($target), $tool, $tool eq 'exonerate' ? 'tsv' : 'psl');
			my $best   = sprintf("%s/%s.%s.%s.best.tsv", $dir, basename($opt->{ref}), basename($target), $tool);
			$ref_jobs .= sprintf("\t%s/submit%s.pl --max-sub 70 --query %s --target %s --output %s --best %s%s%s\n",
				$opt->{binpath}, ucfirst($tool), $opt->{ref}, $target, $output, $best, is_prot_fasta($opt->{ref}) ? $param : '', $opt->{local} ? ' --local' : ''
			);
			$clean_msg .= clean_directory($dir, 0) if ($opt->{restart});
			my $ref_step = $ref_analysis->get_or_create_step( $tool, 'Align the reference sequences on contigs with '.$tool.'.' );
			my $ref_metrics_file = $report_db_folder."/".$ref_step->get_or_create_metrics_filename('alnProtLog');
			$ref_jobs .= "\t$opt->{binpath}/blat2json.pl $best > $ref_metrics_file\n";
		}
		print_msg($opt) if ($opt->{display}); 
	}
	if ($ref_jobs) {
		$opt->{complete} = 0;
		$opt->{cmd} = "if (! -e $dir/00-ASSEMBLY_COMPLETE) then\n\techo 'Assembly not complete; exit'\n\texit\nelse\n${ref_jobs}endif";
		write_shell($opt);
		if ($opt->{restart}) {
			$opt->{restart} = 0;
			$resubmit = $opt->{step};
		}
	}
}

# End
print $clean_msg if ($clean_msg);
print "Assembly complete: nothing to do\n" if ($opt->{restart} && ! $resubmit);
unlink(glob("$dir/00-ASSEMBLY_COMPLETE"), glob("$dir/00-ERROR_AT_*_STEP"), glob("$dir/*transcripts_fpkm_*.fa")) if ($resubmit && $resubmit ne 'reference' && $resubmit ne 'transrate');
quit($opt);

sub quit {
	my $opt = shift;
	# Create report
	create_report( $orig_report_folder, 'DRAP_assessment.html', $dir, $report_db_folder, $report );
	chmod(0775, glob("$opt->{outdir}/*.sh"));
	set_drap_config($opt->{outdir}, $opt);
	# remove renamed shell scripts if they are equal to new ones
	foreach my $script (glob("$opt->{outdir}/*.sh")) {
		next if ($script =~ /\d+:\d+:\d+\.sh/); # previously renamed script
		map { unlink() if (compare($script, $_) == 0) } glob($opt->{outdir}.'/'.basename($script, '.sh').".*-*-*_*:*:*.sh"); # remove previously renamed scripts identical to current one
	}
	exit 0;
}
