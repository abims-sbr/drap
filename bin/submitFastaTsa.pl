#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 submitFastaTsa.pl

=head1 SYNOPSIS

 submitFastaTsa.pl \
   --input SEQ_FILE
   --output FILTERED_FILE \
   --max-seq INT | max-sub INT \
   [--length INTEGER] \
   [--defline] \
   [--contaminant] \
   [--local]
   [--help]

=head1 DESCRIPTION

 Wrapper for fasta_tsa.pl parallel execution.
 The execution is submitted locally or on HPC (see SchedulerFactory).
 
=head1 OPTIONS

=over 8

=item B<-h, --help>

Print help.

=item B<-i, --input>

The path to the sequence file.

=item B<-o, --output>

The path to the filtered file.

=item B<-o, --length>

The contig length cutoff in nucleotide [200].

=item B<-o, --defline>

Reformat the FASTA definition line for TSA submission.

=item B<-o, --contaminant>

Screen all contaminant db. Not only vector/adaptor db.

=item B<--rrna>

Print the output of blastn to rrna contaminant db to this file.
This file is needed for the execution of the NCBI Foreign Contamination 
Screen workflow.

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

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

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
	print $FH_output "# Contigs modified by vecScreen filter:\n" ;
	foreach my $current_file ( @{$in_files} ){
		open( my $FH_input, "<", $current_file.".log" ) or die "Cannot open ".$current_file.".log" ;
		while( my $line = <$FH_input> ){
			print $FH_output $line if ($line =~ /olen=|too (short|noisy)/);
		}
		close( $FH_input );
	}
	close( $FH_output );
}


=head2 procedure gather_rrna

 Usage        : gather( $in_files, $output )
 Function     : Concatenates the content of each rrna blastn file.
 Args         : [array ref] The list of the files to concatenate.
                [str]       The path to the output file.

=cut
sub gather_rrna {
	my ($in_files, $output) = @_ ;
	open( my $FH_output, ">", $output ) or die "Cannot open ".$output ;
	foreach my $current_file ( @{$in_files} ){
		open( my $FH_input, "<", $current_file.".rrna" ) or die "Cannot open ".$current_file.".rrna" ;
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
	
	my $cmd_set = new CmdSet( undef, "parallel", "fastaTsa", $scheduler, 'is' );
	foreach my $input_file ( @{$in_files} ){
		my $output_file = $input_file.".out" ;
		my $current_command = $command ;
		$current_command =~ s/##input##/$input_file/g ;
		$current_command =~ s/##output##/$output_file/g ;
		push( @output_files, $output_file );
		my $current_cmd_obj = new Cmd( $current_command );
		$cmd_set->add_cmd( $current_cmd_obj );
	}
	$cmd_set->submit();
	
	return @output_files
}


=head2 function splitter_by_submission

 Usage        : splitter_by_submission( $input, $max_submissions, 
                $tmp_path_prefix )
 Function     : Split the sequence file in N files.
 Return       : [array] The list of ouput files.
 Args         : [str] The path to the input file.
                [int] The number of splits.
                [str] The temporary file folder path and file prefix.

=cut
sub splitter_by_submission {
	my ($input, $max_submissions, $tmp_path_prefix) = @_ ;
	
	my $nb_seq = 0 ;
	my $FH_input = SequenceFileReader->instantiate( $input );
	while( my $seq_record = $FH_input->next() ){
		$nb_seq++ ;
	}
	$FH_input->close() ;
	my $nb_seq_by_sub = $nb_seq % $max_submissions > 0 ? int($nb_seq/$max_submissions)+1 : $nb_seq/$max_submissions ;
	
	return splitter_by_sequences( $input, $nb_seq_by_sub, $tmp_path_prefix );
}


=head2 function splitter_by_sequences

 Usage        : splitter_by_sequences( $input, $max_sequences, 
                $tmp_path_prefix )
 Function     : Split the sequence file in files with N sequences.
 Return       : [array] The list of ouput files.
 Args         : [str] The path to the input file.
                [int] The number sequences by split.
                [str] The temporary folder path.

=cut
sub splitter_by_sequences {
	my ($input, $max_sequences, $tmp_dir) = @_ ;
	my @output_files = ();
	
	my $FH_input = SequenceFileReader->instantiate( $input );
	my $FH_output = undef ;
	my $writen_nb_seq = $max_sequences ;
	while( my $seq_record = $FH_input->next() ){
		if( $writen_nb_seq == $max_sequences ){
			if( defined($FH_output) ){
				$FH_output->close() ;
			}
			my $output_path = $tmp_dir."/split_".(scalar(@output_files)) ;
			push( @output_files, $output_path );
			$FH_output = undef ;
			if( ref $FH_input eq "FastaIO" ){
				$FH_output = new FastaIO( $output_path, ">" );
			} else {
				$FH_output = new FastqIO( $output_path, ">" );
			}
			$writen_nb_seq = 0 ;
		}
		$FH_output->write( $seq_record );
		$writen_nb_seq++ ;
	}
	if( defined($FH_output) ){
		$FH_output->close() ;
	}
	
	return @output_files ;
}


MAIN:
{
	# Parameters
	my $input = undef ;
	my $output = undef ;
	my $log = undef ;
	my $out_rrna = undef ;
	my $max_submission = undef ;
	my $max_sequences = undef ;
	my $length = 200;
	my $defline = undef;
	my $contaminant = undef;
	my $local = undef ;
	my $help = 0 ;
	GetOptions(
		'i|input=s'   => \$input,
		'o|output=s'  => \$output,
		'l|log=s'     => \$log,
		'rrna=s'      => \$out_rrna,
		'length=i'    => \$length,
		'defline'     => \$defline,
		'contaminant' => \$contaminant,
		'max-seq=i'   => \$max_sequences,
		'max-sub=i'   => \$max_submission,
		'local:i'     => \$local,
		'h|help'      => \$help
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: '-i' is required.") if !defined($input) ;
	pod2usage("$0: '".$input." doesn't exist.") if !(-e $input) ;
	pod2usage("$0: '".$input." is not readable.") if !(-r $input) ;
	pod2usage("$0: '-o' is required.") if !defined($output) ;
	pod2usage("$0: '-l' is required.") if !defined($log) ;
	pod2usage("$0: '-max-seq' or '-max-sub' is required.") if !defined($max_submission) && !defined($max_sequences) ;
	pod2usage("$0: only one of '-max-seq' or '-max-sub' must be used.") if defined($max_submission) && defined($max_sequences) ;

	# Checking
	foreach my $cmd (qw(vecscreen blastn)) {
		`which $cmd 2> /dev/null`;
		die "Command not found: $cmd\n" if ($?);
	}
	my $adaptors_for_screening_euks = $ENV{adaptors_for_screening_euks} || '/usr/local/bioinfo/src/drap/ncbi/adaptors_for_screening_euks.fa';
	my $contam_in_euks              = $ENV{contam_in_euks} || '/usr/local/bioinfo/src/drap/ncbi/contam_in_euks.fa';
	my $mito                        = $ENV{mito} || '/usr/local/bioinfo/src/drap/ncbi/mito.nt';
	my $rrna                        = $ENV{rrna} || '/usr/local/bioinfo/src/drap/ncbi/rrna';
	my @fasta = ($adaptors_for_screening_euks);
	push (@fasta, $contam_in_euks, $mito, $rrna) if defined($contaminant);
	foreach my $fasta (@fasta) {
		if (-f "$fasta.nhr" && -f "$fasta.nin" && -f "$fasta.nsq") {
			print "Found ".basename($fasta)." blast db at: $fasta\n";
		} else {
			die "No such blast db: $fasta\n";
		}
	}
	if (defined($contaminant) && dirname($contam_in_euks) ne dirname($mito) && dirname($contam_in_euks) ne dirname($rrna)) {
		die "Blast dbs ".basename($contam_in_euks).", ".basename($mito)." and ".basename($rrna). " must be in the same directory\n";
	}

	# Process
	my ($filename, $dir, $ext) = fileparse($output, qr/\.[^.]*/);
	my $tmp_folder = File::Spec->rel2abs($dir)."/tmp_".time."_".int(rand(10000)) ;
	mkdir( $tmp_folder );
	# Set scheduler
	my $scheduler;
	if ( defined($local) ){
		$scheduler = SchedulerFactory->instantiate($tmp_folder, 'local', $local || $ENV{scheduler_local_cpu} || 1);
	} else {
		$scheduler = SchedulerFactory->instantiate($tmp_folder);
	}
	my @in_files = undef ;
	if( defined($max_submission) ) {
		@in_files = splitter_by_submission( $input, $max_submission, $tmp_folder );
	} else {
		@in_files = splitter_by_sequences( $input, $max_sequences, $tmp_folder );
	}
	my @params = ('--minsize', $length, '--db', $adaptors_for_screening_euks);
	push(@params, '--defline') if defined($defline);
	push(@params, '--contaminant-dir ', dirname($contam_in_euks), ' --rrna ##output##.rrna') if defined($contaminant);
	my @out_files = scatter( \@in_files, 'fasta_tsa.pl '.join(' ',@params).' --log ##output##.log ##input## > ##output##', $scheduler );
	gather_seq( \@out_files, $output );
	gather_log( \@out_files, $log );
	gather_rrna( \@out_files, $out_rrna ) if defined($contaminant);
	
	# Remove temmporary files and folder
	for( my $idx = 0 ; $idx < scalar(@in_files) ; $idx++ ){
		unlink $in_files[$idx] ;
		unlink $out_files[$idx] ;
		unlink $out_files[$idx].".log" ;
		unlink $out_files[$idx].".rrna" ;
	}
	rmdir( $tmp_folder );
}