#!/usr/bin/perl
# $Id$

=pod

=head1 NAME

 bwa

=head1 SYNOPSIS

 bwa.pl \
   --bam BAM_OUTPUT
   --reference REFERENCE_FASTA \
   --R1 R1_FILE [... R1_FILE_n] \
   [--R2 R2_FILE  [... R2_FILE_n]] \
   [--algorithm BWA_ALGORITHM] \
   [--mem MEMORY_GB] \
   [--vmem MEMORY_GB] \
   [--nb-threads NB_CPU] \
   [--local]
   [--help]

=head1 DESCRIPTION

 Launches bwa:
   1- Executes "bwa index" if the reference is not indexed
   2- for each R1 executes the alignment:
        aln then samse or aln R1, aln R2, then sampe
      or
        bwasw
      or
        mem
   3- Executes "samtools sort"
   4- Executes "samtools merge" if it exists several R1.
 The execution is submitted locally or on HPC (see SchedulerFactory).
 
=head1 OPTIONS

=over 8

=item B<-a, --algorithm>

The used bwa algorithm for alignment ("aln", "bwasw", "mem").

=item B<-b, --bam>

The path to the BAM produced.

=item B<-h, --help>

Print help.

=item B<--mem>

The maximum memory used by the alignment (default: 4G).

=item B<--R1>

The list of path to the R1 file(s).

=item B<--R2>

The list of path to the R2 file(s). The R2 files must be in the same order than
R1.

=item B<-r, --reference>

The path to the fasta reference.

=item B<-t, --nb-threads>

The maximum number of CPU cores used by the alignment (default: 4).

=item B<--vmem>

The maximum virtual memory used by the alignment (default: 8G).

=item B<--local>

Run locally (no qsub). Without numeric option, use 1 CPU core
or env variable scheduler_local_cpu if defined.
Specify --local -1, to use all available CPU cores.
Specify --local n, to use n CPU cores.

=back

=head1 VERSION

 2.0.0

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

MAIN:
{
	# Parameters
	my $reference = undef ;
	my @R1 = () ;
	my @R2 = () ;
	my $bam = undef ;
	my $algorithm = "aln" ;
	my $cpu = 1 ;
	my $vmem = "8G" ;
	my $mem = "4G" ;
	my $local = undef ;
	my $help = 0 ;
	GetOptions(
		'r|reference=s'  => \$reference,
		'R1=s{1,}'       => \@R1,
		'R2=s{1,}'       => \@R2,
		'b|bam=s'        => \$bam,
		'a|algorithm=s'  => \$algorithm,
		't|nb-threads=s' => \$cpu,
		'vmem=s'         => \$vmem,
		'mem=s'          => \$mem,
		'local:i'        => \$local,
		'h|help'         => \$help
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|VERSION|OPTIONS"
	) if ($help) ;
	pod2usage("$0: '-b' is required.") if !defined($bam) ;
	pod2usage("$0: '-r' is required.") if !defined($reference) ;
	pod2usage("$0: '".$reference." doesn't exist.") if !(-e $reference) ;
	pod2usage("$0: '".$reference." is not readable.") if !(-r $reference) ;
	pod2usage("$0: '-R1' is required.") if scalar(@R1) == 0 ;
	@R1 = split(/,/, join(',', @R1));
	foreach my $current_R1 (@R1){
		pod2usage("$0: '".$current_R1." doesn't exist.") if !(-e $current_R1) ;
		pod2usage("$0: '".$current_R1." is not readable.") if !(-r $current_R1) ;
	}
	@R2 = split(/,/, join(',', @R2));
	if( scalar(@R2) > 0 ) {
		foreach my $current_R2 (@R2){
			pod2usage("$0: '".$current_R2." doesn't exist.") if !(-e $current_R2) ;
			pod2usage("$0: '".$current_R2." is not readable.") if !(-r $current_R2) ;
		}
	}
	my ($bam_filename, $out_dir, $ext) = fileparse( File::Spec->rel2abs($bam), qr/\.[^.]*/ );
	# Set scheduler
	my $scheduler;
	if ( defined($local) ){
		$scheduler = SchedulerFactory->instantiate($out_dir, 'local', $local || $ENV{scheduler_local_cpu} || $cpu || 1);
	} else {
		$scheduler = SchedulerFactory->instantiate($out_dir);
	}

	# Index
	if( ! -e $reference.".bwt" ){
		my ($filename, $dir, $ext) = fileparse( $reference, qr/\.[^.]*/ );
		my $new_reference = $out_dir."/".$filename.$ext ;
		my $cmd = new CmdSet( undef, undef, "bwa.pl-0", $scheduler );
		$cmd->add_cmd( new Cmd("ln -sf ".$reference." ".$new_reference, "", $mem, $vmem) );
		$cmd->add_cmd( new Cmd("bwa index ".$new_reference, "", $mem, $vmem) );
		$cmd->submit();
		$reference = $new_reference ;
	}
	
	# Create BAMs
	my @unmerged_bams = () ;
	for( my $file_idx=0 ; $file_idx < scalar(@R1) ; $file_idx++ ){
		my $sam = $bam."_unsorted_".$file_idx.".sam" ;
		
		# Alignment
		if( $algorithm eq "aln" ){
			# aln
			my ($R1_filename, $R1_dir, $R1_ext) = fileparse( $R1[$file_idx], qr/\.[^.]*/ );
			my $sai_R1 = $out_dir."/".$bam_filename."_".$R1_filename.".sai" ;
			my $sai_R2 = undef ;
			my @cmd_list = ();
			push( @cmd_list, new Cmd("bwa aln -t ".$cpu." ".$reference." ".$R1[$file_idx]." > ".$sai_R1, $cpu, $mem, $vmem) );
			if( scalar(@R2) > $file_idx ){
				my ($R2_filename, $R2_dir, $R2_ext) = fileparse( $R2[$file_idx], qr/\.[^.]*/ );
				$sai_R2 = $out_dir."/".$bam_filename."_".$R2_filename.".sai" ;
				push( @cmd_list, new Cmd("bwa aln -t ".$cpu." ".$reference." ".$R2[$file_idx]." > ".$sai_R2, $cpu, $mem, $vmem) );
			}
			my $cmd_aln = new CmdSet( \@cmd_list, "parallel", "bwa.pl-1", $scheduler );
			$cmd_aln->submit();
			
			# sampe/samse
			my $second_cmd = new CmdSet( undef, undef, "bwa.pl-2", $scheduler );
			if( scalar(@R2) <= $file_idx ){
				$second_cmd->add_cmd( new Cmd("bwa samse ".$reference." ".$sai_R1." ".$R1[$file_idx]." > ".$sam, $cpu, $mem, $vmem) );
				$second_cmd->submit();
				unlink( $sai_R1 );
			} else {
				$second_cmd->add_cmd( new Cmd("bwa sampe ".$reference." ".$sai_R1." ".$sai_R2." ".$R1[$file_idx]." ".$R2[$file_idx]." > ".$sam, $cpu, $mem, $vmem) );
				$second_cmd->submit();
				unlink( $sai_R1 );
				unlink( $sai_R2 );
			}
		} else {
			my $cmd = new CmdSet( undef, undef, "bwa.pl-1", $scheduler );
			if( scalar(@R2) <= $file_idx ){
				$cmd->add_cmd( new Cmd("bwa ".$algorithm." -t ".$cpu." ".$reference." ".$R1[$file_idx]." > ".$sam, $cpu, $mem, $vmem) );
			} else {
				$cmd->add_cmd( new Cmd("bwa ".$algorithm." -t ".$cpu." ".$reference." ".$R1[$file_idx]." ".$R2[$file_idx]." > ".$sam, $cpu, $mem, $vmem) );
			}
			$cmd->submit();
		}
	
		# SAM to BAM
		my $unsorted_bam = $sam ;
		$unsorted_bam =~ s/\.sam$/\.bam/ ;
		push( @unmerged_bams, $bam."_".$file_idx ); 
		my $cmd = new CmdSet( undef, undef, "bwa.pl-3", $scheduler );
		$cmd->add_cmd( new Cmd("samtools view -bS -@ ".$cpu." ".$sam." > ".$unsorted_bam, $cpu, $mem, $vmem) );
		$cmd->add_cmd( new Cmd("samtools sort -@ ".$cpu." -O bam -o ".$bam." ".$unsorted_bam, $cpu, $mem, $vmem) );
		$cmd->add_cmd( new Cmd("mv ".$bam." ".$unmerged_bams[$file_idx], $cpu, $mem, $vmem) );
		$cmd->submit();
		unlink( $unsorted_bam );
		unlink( $sam );
	}

	# BAMs merge
	my $cmd = new CmdSet( undef, undef, "bwa.pl-4", $scheduler );
	if( scalar(@unmerged_bams) > 1 ){
		$cmd->add_cmd( new Cmd("samtools merge ".$bam." ".join(" ", @unmerged_bams)) );
		$cmd->submit();
		unlink( @unmerged_bams );
	} else {
		$cmd->add_cmd( new Cmd("mv ".$unmerged_bams[0]." ".$bam ) );
		$cmd->submit();
	}
}