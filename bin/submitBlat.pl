#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 submitBlat.pl

=head1 SYNOPSIS

 submitBlat.pl \
   --query SEQ_FILE
   --target SEQ_FILE
   --output PSL_FILE \
   --max-seq INT | max-sub INT \
   [--params 'BLAT_PARAMS'] \
   [--best TSV_FILE] \
   [--mem MEMORY_GB] \
   [--vmem MEMORY_GB] \
   [--local]
   [--help]

=head1 DESCRIPTION

 Wrapper for blat parallel execution.
 The execution is submitted locally or on HPC (see SchedulerFactory).
 
=head1 OPTIONS

=over 8

=item B<-q, --query>

The path to the query fasta file.

=item B<-t, --target>

The path to the target fasta file.

=item B<-o, --output>

The path to the output psl file or the output directory.

=item B<-p, --params>

Additional Blat parameters between single quotes.
Ex: -p '-mask=lower -qMask=lower'

=item B<-b, --best>

The path to the "best hit for each query" output file.
blat_parser.pl must be in PATH.

=item B<--max-seq>

The max number of sequences per submission.

=item B<--max-sub>

The max number of submissions.

=item B<--mem>

The maximum memory used by each command on a sub-file (default: 4G).

=item B<--vmem>

The maximum virtual memory used by each command on a sub-file (default: 8G).

=item B<--local>

Run locally (no qsub). Without numeric option, use 1 CPU core
or env variable scheduler_local_cpu if defined.
Specify --local -1, to use all available CPU cores.
Specify --local n, to use n CPU cores.

=item B<-h, --help>

Print help.

=back

=head1 VERSION

 0.2.1

=head1 AUTHOR

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2016 INRA

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



=head2 procedure gather

 Usage        : gather( $in_files, $output )
 Function     : Concatenates the content of each files.
 Args         : [array ref] The list of the files to concatenate.
                [str]       The path to the output file.

=cut
sub gather {
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


=head2 function scatter

 Usage        : scatter( $in_files, $command, $scheduler, $mem, $vmem )
 Function     : Submit the command on each in_files.
 Return       : [array] The list of ouput files.
 Args         : [array ref] The list of the input files.
                [str]       The command to execute. In this string the input 
                            path is replaced by ##input## and the output path
                            is replaced by ##output##.
                [AbstractScheduler] The scheduler used to submit commands.
                [str]       The maximum memory used by each command on a 
                            sub-file.
                [str]       The maximum virtual memory used by each command on 
                            a sub-file. 

=cut
sub scatter {
	my ($in_files, $command, $scheduler, $mem, $vmem) = @_ ;
	my @output_files = () ;
	
	my $cmd_set = new CmdSet( undef, "parallel", "submitBlat", $scheduler, 'is' );
	foreach my $input_file ( @{$in_files} ){
		my $output_file = $input_file.".out" ;
		my $current_command = $command ;
		$current_command =~ s/##input##/$input_file/g ;
		$current_command =~ s/##output##/$output_file/g ;
		push( @output_files, $output_file );
		my $current_cmd_obj = new Cmd( $current_command, "", $mem, $vmem );
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
	my $query = undef ;
	my $target = undef ;
	my $output = undef ;
	my $params = "" ;
	my $best = undef ;
	my $max_submission = undef ;
	my $max_sequences = undef ;
	my $vmem = "8G" ;
	my $mem = "4G" ;
	my $local = undef ;
	my $help = 0 ;
	GetOptions(
		'q|query=s'  => \$query,
		't|target=s' => \$target,
		'o|output=s' => \$output,
		'p|params=s' => \$params,
		'b|best=s'   => \$best,
		'max-seq=i'  => \$max_sequences,
		'max-sub=i'  => \$max_submission,
		'vmem=s'     => \$vmem,
		'mem=s'      => \$mem,
		'local:i'    => \$local,
		'h|help'     => \$help
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: '-q' is required.") if !defined($query) ;
	pod2usage("$0: '".$query." doesn't exist.") if !(-e $query) ;
	pod2usage("$0: '".$query." is not readable.") if !(-r $query) ;
	pod2usage("$0: '-t' is required.") if !defined($target) ;
	pod2usage("$0: '".$target." doesn't exist.") if !(-e $target) ;
	pod2usage("$0: '".$target." is not readable.") if !(-r $target) ;
	pod2usage("$0: '-o' is required.") if !defined($output) ;
	pod2usage("$0: '-max-seq' or '-max-sub' is required.") if !defined($max_submission) && !defined($max_sequences) ;
	pod2usage("$0: only one of '-max-seq' or '-max-sub' must be used.") if defined($max_submission) && defined($max_sequences) ;

	# Checking
	if ($best) {
		`which blat_parser.pl 2> /dev/null`;
		die "Command not found: blat_parser.pl\n" if ($?);
	}
	my @params = ('-noHead');
	push(@params, split(/ /, $params));

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
		@in_files = splitter_by_submission( $query, $max_submission, $tmp_folder );
	} else {
		@in_files = splitter_by_sequences( $query, $max_sequences, $tmp_folder );
	}
	my @out_files = scatter( \@in_files, "blat $target ##input## ".join(' ',@params).' ##output##', $scheduler, $mem, $vmem );
	gather( \@out_files, $output );
	`cat $output | blat_parser.pl > $best` if ($best);

	# Remove temmporary files and folder
	for( my $idx = 0 ; $idx < scalar(@in_files) ; $idx++ ){
		unlink $in_files[$idx] ;
		unlink $out_files[$idx] ;
	}
	rmdir( $tmp_folder );
}
