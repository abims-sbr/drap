#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 submitSamCorrectVar.pl

=head1 SYNOPSIS

 submitSamCorrectVar.pl \
   --fasta SEQ_FILE \
   --bam BAM_FILE[,...,BAM_FILE_n] \
   --output CORRECTED_FILE \
   --log LOG_FILE \
   --max-seq INT | max-sub INT \
   [--local]
   [--help]

=head1 DESCRIPTION

 Wrapper for samCorrectVariation.pl parallel execution.
 The execution is submitted locally or on HPC (see SchedulerFactory).
 
=head1 OPTIONS

=over 8

=item B<-h, --help>

Print help.

=item B<-i, --fasta>

The path to the sequence fasta file.

=item B<-i, --bam>

The path to the bam file(s).

=item B<-o, --output>

The path to the corrected file.

=item B<-o, --log>

The path to the log file.

=item B<--max-seq>

The max number of sequences per submission.

=item B<--max-sub>

The max number of submissions.

=item B<--local>

Run locally (no qsub). Without numeric option, use 1 CPU core
or env variable scheduler_local_cpu if defined.
Specify --local -1, to use all available CPU cores.
Specify --local n, to use n CPU cores.

=back

=head1 VERSION

 0.1.0

=head1 AUTHOR

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2017 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict ;
use warnings ;
use List::Util ;
use File::Spec ;
use File::Basename;
use Pod::Usage ;
use Getopt::Long ;
use FindBin ;
use lib ("$FindBin::Bin") ;
use Cwd 'getcwd' ;
use Cmd ;
use CmdSet ;
use SequenceFileReader ;



=head2 procedure gather_seq

 Usage        : gather( $in_files, $output )
 Function     : Concatenates the content of each sequence file.
 Args         : [array ref] The list of the files to concatenate.
                [str]       The path to the output file.

=cut
sub gather_seq {
	my ($in_files, $output) = @_ ;

	open( my $FH_output, ">", $output ) or die "Cannot open ".$output ;
	foreach my $current_file ( @{$in_files} ){
		open( my $FH_input, "<", $current_file ) or die "Cannot open ".$current_file ;
		while( my $line = <$FH_input> ){
			print $FH_output $line ;
		}
		close( $FH_input );
	}
	close( $FH_output );
}


=head2 procedure gather_log

 Usage        : gather( $in_files, $output )
 Function     : Concatenates the content of each log file.
 Args         : [array ref] The list of the files to concatenate.
                [str]       The path to the output file.

=cut
sub gather_log {
	my ($in_files, $output) = @_ ;
	open( my $FH_output, ">", $output ) or die "Cannot open ".$output ;
	foreach my $current_file ( @{$in_files} ){
		open( my $FH_input, "<", $current_file.".log" ) or die "Cannot open ".$current_file.".log" ;
		while( my $line = <$FH_input> ){
			print $FH_output $line;
		}
		close( $FH_input );
	}
	close( $FH_output );
}


=head2 function scatter

 Usage        : scatter( $in_files, $command, $cheduler )
 Function     : Submit the command on each in_files.
 Return       : [array] The list of ouput files.
 Args         : [array ref] The list of the input files.
                [str]       The command to execute. In this string the input 
                            path is replaced by ##input## and the output path
                            is replaced by ##output##.
                [AbstractScheduler] The scheduler used to submit commands.


=cut
sub scatter {
	my ($in_files, $command, $scheduler) = @_ ;
	my @output_files = () ;
	
	my $cmd_set = new CmdSet( undef, "parallel", "samCorrectVar", $scheduler, 'is' );
	foreach my $input_file ( @{$in_files} ){
		my $output_file = $input_file.".out" ;
		my $current_command = $command ;
		$current_command =~ s/##inputFasta##/$input_file.fa/g ;
		$current_command =~ s/##inputBam##/$input_file.bam/g ;
		$current_command =~ s/##output##/$output_file/g ;
		push( @output_files, $output_file );
		my $current_cmd_obj = new Cmd( $current_command );
		$cmd_set->add_cmd( $current_cmd_obj );
	}
	$cmd_set->submit();
	
	return @output_files
}


=head2 function splitter_by_submission

 Usage        : splitter_by_submission( $inputFasta, $inputBam, $max_submissions, 
                $tmp_path_prefix )
 Function     : Split the sequence file in N files.
 Return       : [array] The list of ouput files.
 Args         : [str] The path to the input file.
                [ref] Reference to the array containing path to the bam files.
                [int] The number of splits.
                [str] The temporary file folder path and file prefix.

=cut
sub splitter_by_submission {
	my ($inputFasta, $inputBam, $max_submissions, $tmp_path_prefix) = @_ ;
	
	my $nb_seq = 0 ;
	my $FH_input = SequenceFileReader->instantiate( $inputFasta );
	while( my $seq_record = $FH_input->next() ){
		$nb_seq++ ;
	}
	$FH_input->close() ;
	my $nb_seq_by_sub = $nb_seq % $max_submissions > 0 ? int($nb_seq/$max_submissions)+1 : $nb_seq/$max_submissions ;
	
	return splitter_by_sequences( $inputFasta, $inputBam, $nb_seq_by_sub, $tmp_path_prefix );
}


=head2 function splitter_by_sequences

 Usage        : splitter_by_sequences( $inputFasta, $inputBam, $max_sequences, 
                $tmp_path_prefix )
 Function     : Split the sequence file in files with N sequences.
 Return       : [array] The list of ouput files.
 Args         : [str] The path to the input file.
                [ref] Reference to the array containing path to the bam files.
                [int] The number sequences by split.
                [str] The temporary folder path.

=cut

sub splitter_by_sequences {
	my ($inputFasta, $inputBam, $max_sequences, $tmp_dir) = @_ ;
	my @output_files = ();
	
	my $FH_inputFasta = SequenceFileReader->instantiate( $inputFasta );
	my %inputFastaSeqs = $FH_inputFasta->getAll();
	open(my $FH_inputBam, 'samtools merge - '.join(' ', @$inputBam).' | samtools view - |') or die "Can't merge bam files: ".join(', ', @$inputBam)."\n$!";
	my $FH_outputFasta = undef ;
	my $FH_outputBam = undef ;
	my $written_nb_seq = $max_sequences ;
	my @bam_record = undef ;
	my $written_ref = '';
	while( my $bam_record = <$FH_inputBam> ){
		@bam_record = split("\t", $bam_record);
		next if ($bam_record[2] eq '*');
		if( $written_nb_seq == $max_sequences && $written_ref ne $bam_record[2]){
			if( defined($FH_outputBam) ){
				$FH_outputFasta->close();
				close($FH_outputBam);
			}
			my $output_path = $tmp_dir."/split_".(scalar(@output_files)) ;
			push( @output_files, $output_path );
			$FH_outputFasta = new FastaIO( "$output_path.fa", ">" );
			open($FH_outputBam, "|-", "samtools view -Sb -o $output_path.bam - ") or die "Can't write bam file: $output_path.bam\n$!";
			open(my $FH_bamHeader, "-|", "samtools view -H ".$inputBam->[0]) or die "Can't open bam file: ".$inputBam->[0].".bam\n$!";
			while ( my $bam_record = <$FH_bamHeader> ){
				print $FH_outputBam $bam_record;
			}
			close($FH_bamHeader);
			$written_nb_seq = 0 ;
		}
		@bam_record = split("\t", $bam_record);
		print $FH_outputBam $bam_record;
		unless ($written_ref eq $bam_record[2]) {
			$FH_outputFasta->write($inputFastaSeqs{$bam_record[2]});
			delete($inputFastaSeqs{$bam_record[2]});
			$written_nb_seq++ ;
		}
		$written_ref = $bam_record[2];
	}
	# write inputFastaSeq without bam record in the last fasta split
	if( %inputFastaSeqs ) {
		foreach my $inputFastaSeq ( keys %inputFastaSeqs ){
			$FH_outputFasta->write($inputFastaSeqs{$inputFastaSeq});
		}
	}
	$FH_outputFasta->close() ;
	close($FH_outputBam);

	return @output_files ;
}

MAIN:
{
	# Parameters
	my $fasta = undef ;
	my @bams = () ;
	my $output = undef ;
	my $log = undef ;
	my $max_submission = undef ;
	my $max_sequences = undef ;
	my $local = undef ;
	my $help = 0 ;
	GetOptions(
		'f|fasta=s'   => \$fasta,
		'b|bam=s{1,}' => \@bams,
		'o|output=s'  => \$output,
		'l|log=s'     => \$log,
		'max-seq=i'   => \$max_sequences,
		'max-sub=i'   => \$max_submission,
		'local:i'     => \$local,
		'h|help'      => \$help
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: '-f' is required.") if !defined($fasta) ;
	foreach my $file ($fasta, @bams) {
		pod2usage("$0: '".$file." doesn't exist.") if !(-e $file) ;
		pod2usage("$0: '".$file." is not readable.") if !(-r $file) ;
	}
	pod2usage("$0: '-o' is required.") if !defined($output) ;
	pod2usage("$0: '-l' is required.") if !defined($log) ;
	pod2usage("$0: '-max-seq' or '-max-sub' is required.") if !defined($max_submission) && !defined($max_sequences) ;
	pod2usage("$0: only one of '-max-seq' or '-max-sub' must be used.") if defined($max_submission) && defined($max_sequences) ;

	# Process
	@bams = split(/,/, join(',', @bams));
	my ($filename, $dir, $ext) = fileparse($output, qr/\.[^.]*/);
	my $tmp_folder = File::Spec->rel2abs($dir)."/tmp_".time."_".int(rand(10000));
	$fasta  = File::Spec->rel2abs($fasta);
	for (@bams) { $_ = File::Spec->rel2abs($_) };
	$output = File::Spec->rel2abs($output);
	$log    = File::Spec->rel2abs($log);
	my $current_folder = getcwd;
	mkdir( $tmp_folder );
	# Set scheduler
	my $scheduler;
	if ( defined($local) ){
		$scheduler = SchedulerFactory->instantiate($tmp_folder, 'local', $local || $ENV{scheduler_local_cpu} || 1);
	} else {
		$scheduler = SchedulerFactory->instantiate($tmp_folder);
	}
	# need to go into tmp_folder to clean tmp files forgotten by Bio::Tools::Run::StandAloneBlastPlus 
	chdir( $tmp_folder );
	my @in_files = undef ;
	if( defined($max_submission) ) {
		@in_files = splitter_by_submission( $fasta, \@bams, $max_submission, '.' );
	} else {
		@in_files = splitter_by_sequences( $fasta, \@bams, $max_sequences, '.' );
	}
	my @out_files = scatter( \@in_files, 'samtools mpileup -f ##inputFasta## ##inputBam## | samCorrectVariation.pl --log ##output##.log ##inputFasta## > ##output##', $scheduler );
	gather_seq( \@out_files, $output );
	gather_log( \@out_files, $log );
	
	# Remove temmporary files and folder
	for( my $idx = 0 ; $idx < scalar(@in_files) ; $idx++ ){
		unlink $in_files[$idx].".fa" ;
		unlink $in_files[$idx].".fa.fai" ;
		unlink $in_files[$idx].".bam" ;
		unlink $out_files[$idx] ;
		unlink $out_files[$idx].".log" ;
	}
	chdir( $current_folder );
	rmdir( $tmp_folder );
}