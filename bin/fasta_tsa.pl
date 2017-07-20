#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use IPC::Run qw(run timeout);

my $opt = {};
$opt->{contaminant} = {
	'db' => ['contam_in_euks.fa', 'mito.nt', 'rrna'],
	'blastn_cmd_line' => {
		'contam_in_euks.fa' => [qw(blastn -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt 7 -db)],
		'mito.nt' => [qw(blastn -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 98.6 -soft_masking true -outfmt 7 -db)],
		'rrna' => [qw(blastn -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true -outfmt 7 -db)]
	},
	'valid_hit' => {
		'contam_in_euks.fa' => sub { return 1 if (($_[2]>=98.0 && $_[3]>=50)||($_[2]>=94.0 && $_[3]>=100)||($_[2]>=90.0 && $_[3]>=200)); return 0 },
		'mito.nt' => sub {return 1 if ($_[3]>=120); return 0 },
		'rrna' => sub {return 1 if ($_[3]>=100); return 0 }
	}
};

#------------------------------------------------------------
# options and parameters
#------------------------------------------------------------
# Options definitions (name,format,default,required)
my @getOptions = ( ['help',undef,undef,0],
	['man',undef,undef,0],
	['d|db','s',undef,0],
	['c|contaminant-dir','s',undef,0],
	['p|npercent','i',10,0],
	['g|ngap','i',10,0],
	['s|minsize','i',200,0],
	['defline',undef,0,0],
	['e|endtrim','i',20,0],
	['t|ntrim','i',10,0],
	['f|FCSreport','s',undef,0],
	['l|log','s',undef,0],
	['rrna','s',undef,0]
);
# defaults values
my $required = {};
map {$opt->{$_->[0]}=$_->[2];$required->{$_->[0]}=$_->[3];} @getOptions;
# build options list
my %getOptions = ();
map {$getOptions{defined($_->[1])?$_->[0].'='.$_->[1]:$_->[0]}=\$opt->{$_->[0]};} @getOptions;
# retrieve options
GetOptions(%getOptions) || pod2usage(2);
pod2usage(1) if ($opt->{'help'});
pod2usage(-exitstatus => 0, -verbose => 2) if ($opt->{'man'});

# check required
map {pod2usage("$_ is required") if ($required->{$_}&&!defined($opt->{$_}));} keys %$opt;

# restrict keys of %opt to long option names - key 'o|outdir' become 'outdir' and replace - by _
map { my $key = $_; s/.+\|//; s/-/_/g; $opt->{$_} = delete($opt->{$key}) } keys %$opt;

# read input parameters
if (defined($opt->{FCSreport})) {
	pod2usage("No such FCSreport $opt->{FCSreport}") unless (-e $opt->{FCSreport});
} else {
	pod2usage("d|db is required") if (!defined($opt->{db}));
	pod2usage("No such vecscreen db $opt->{db}") if (!-e $opt->{db});
	map {pod2usage("Can't find file $opt->{db}.$_\n") unless (-e "$opt->{db}.$_")} ('nhr', 'nin', 'nsq');
}
if (defined($opt->{contaminant_dir})) {
	map { 
		my $db = $_;
		map { pod2usage("Can't find file $opt->{contaminant_dir}/$db.$_\n") unless (-e "$opt->{contaminant_dir}/$db.$_") } ('nhr', 'nin', 'nsq');
		push(@{$opt->{contaminant}->{blastn_cmd_line}->{$db}}, "$opt->{contaminant_dir}/$db");
	} @{$opt->{contaminant}->{db}};
	if (defined($opt->{rrna})) { open(RRNA, ">$opt->{rrna}") or die "Can't open file $opt->{rrna}"; $opt->{rrnafh} = *RRNA; }
}
$#ARGV==0 || pod2usage("Missing reference file as parameter");

# log file
my $logfh= *STDERR;
if (defined($opt->{log})) { open(LOG, ">$opt->{log}") or die "Can't open file $opt->{log}"; $logfh= *LOG; }

# vecscreen command
my $version = `vecscreen -version 2> /dev/null`;
if ($version =~ /vecscreen: 2.0.0/) {
	$opt->{vecscreen_cmd_line} = \@{["vecscreen", "-db", $opt->{db}, "-outfmt", "1", "-text_output"]};
} else {
	my $help = `vecscreen --help 2> /dev/null`;
	if ($help =~ /-f\s+Output format:/ms) {
		$opt->{vecscreen_cmd_line} = \@{["vecscreen", "-d", $opt->{db}, "-f", "3"]};
	} else {
		die "Unknown vecscreen version\n";
	}
}

#------------------------------------------------------------
# read reference sequences and store it into a hash
#------------------------------------------------------------
my %REF = ();
my @REF = ();
open(REF,$ARGV[0]) || die "Can't open file $ARGV[0]";

my $ref = '';
while (my $line = <REF>) {
	chomp $line;
	if ($line =~ /^>(\S+)/) {
		if ($opt->{defline} && length($1) > 50) {
			die "The sequence identifier cannot exceed 50 characters: $1";
		}
		if (!exists($REF{$1})) {
			$ref = $1;
			$REF{$ref} = '';
			push(@REF, $ref);
		}
		else {
			die "Sequence $1 already seen";
		}
	}
	else {
		$line =~ s/\s//g;
		$REF{$ref} .= uc($line);
	}
}
close(REF);

if (defined($opt->{FCSreport})) {
	correct_reference($opt, \%REF, \@REF, $logfh);
} else {
	validate_reference($opt, \%REF, \@REF, $logfh);
}


#------------------------------------------------------------
# print reference sequences
#------------------------------------------------------------
my $moltype = $opt->{defline} ? ' [moltype=transcribed_RNA]' : '';
foreach $ref (@REF) {
	next unless (exists $REF{$ref});
	$REF{$ref}=~s/(.{60})/$1\n/g or warn "$ref\n";
	chomp $REF{$ref};
	print ">$ref$moltype\n$REF{$ref}\n";
}

sub correct_reference {
	my ($opt, $REF, $sortREF, $logfh) = @_;
	
	my @applied_corrections = ();
	my $corrections = get_corrections($opt->{FCSreport});
	foreach my $correction (split("\n", $corrections)) {
		my ($ref, $i, $d, $t) = split("\t", $correction);
		my $reflen = length($REF->{$ref});
		foreach my $interval (split(",", $i)) {
			my ($b, $e) = split(/\.\./, $interval);
			my $m = 1 + $e - $b;
			substr($REF->{$ref}, $b - 1, $m, $t x $m);
		}
		$applied_corrections[0] = "$i:$d";
		my $nvec += $REF->{$ref} =~ tr/V/N/; 
		my $ncontam += $REF->{$ref} =~ tr/M/N/;
		my $nmasked += $REF->{$ref} =~ tr/N/N/;

		# trim ends if gaps w/in --endtrim of ends
		my $ntrim = trim_ends($opt, $ref, $REF);

		# cut stretch of at least --ngap Ns
		my $ncut = cut_N_strech($opt, $ref, $REF);
		print_log($opt, $ref, $REF, $reflen, \@applied_corrections, $nvec, $ncontam, $nmasked, $ntrim, $ncut, $logfh);
	}
}

sub validate_reference {
	my ($opt, $REF, $sortREF, $logfh) = @_;

	my $rrnafh = defined($opt->{rrna}) ? $opt->{rrnafh} : undef;
	foreach my $ref (@$sortREF) {
		my $reflen = length($REF->{$ref}); 
		my $pass = 0;
		my @all_matches = ();
		my $nvec = my $ncontam = my $ntrim = my $ncut = my $nmasked = 0;
		my $matches = get_vecscreen_match($ref, $REF);
		$matches .= get_contaminant_match($ref, $REF) if (defined($opt->{contaminant_dir}));
		while (1) {
			$pass++;
			# mask vecsreen and contaminant matches
			if ($matches) {
				push(@all_matches, "Pass$pass");
				foreach my $match (split("\n", $matches)) {
					my ($b, $e, $d, $t) = split("\t", $match);
					my $m = 1 + $e - $b;
					# if region not already masked, keep match record and mask region
					if (substr($REF->{$ref}, $b - 1, $m) =~ /[ATGC]/) {
						push(@all_matches, "$b-$e:$d");
						substr($REF->{$ref}, $b - 1, $m, $t x $m);
					}
				}
			}
			$nvec += $REF->{$ref} =~ tr/V/N/; 
			$ncontam += $REF->{$ref} =~ tr/M/N/;
			$nmasked += $REF->{$ref} =~ tr/N/N/;

			# trim ends if gaps w/in --endtrim of ends
			$ntrim += trim_ends($opt, $ref, $REF);

			# cut stretch of at least --ngap Ns
			$ncut += cut_N_strech($opt, $ref, $REF);

			# check if reference still exists, check again vecscreen and contaminant matches and redo if needed
			last unless (length($REF->{$ref}) > 0);
			$matches = get_vecscreen_match($ref, $REF);
			$matches .= get_contaminant_match($ref, $REF, $rrnafh) if (defined($opt->{contaminant_dir}));
			last unless ($matches);
		}
		print_log($opt, $ref, $REF, $reflen, \@all_matches, $nvec, $ncontam, $nmasked, $ntrim, $ncut, $logfh);
	}
}

sub trim_ends {
	my ($opt, $ref, $REF) = @_;
	my $trimlen = $opt->{ntrim};
	my $endtrim = $opt->{endtrim};
	my $trim = ('N') x $trimlen;
	my $ntrim = 0;
	$ntrim += length($1) if ($REF->{$ref} =~ s/(N+)$//);
	for (my $trim3 = rindex($REF->{$ref}, $trim); $REF->{$ref} && $trim3 >= length($REF->{$ref}) - $endtrim; ) {
		$ntrim += length($REF->{$ref}) - $trim3;
		$REF->{$ref} = substr($REF->{$ref}, 0, $trim3);
		$ntrim += length($1) if ($REF->{$ref} =~ s/(N+)$//);
		$trim3 = rindex($REF->{$ref}, $trim);
	}
	$ntrim += length($1) if ($REF->{$ref} =~ s/^(N+)//);
	for (my $trim5 = index($REF->{$ref}, $trim); $REF->{$ref} && $trim5 >= 0 && $trim5 + $trimlen <= $endtrim; ) {
		$ntrim += $trim5 + $trimlen;
		$REF->{$ref} = substr($REF->{$ref}, $trim5 + $trimlen);
		$ntrim += length($1) if ($REF->{$ref} =~ s/^(N+)//);
		$trim5 = index($REF->{$ref}, $trim);
	}
	return $ntrim;
}


sub cut_N_strech {
	my ($opt, $ref, $REF) = @_;
	my $gaplen = $opt->{ngap};
	my $gap = ('N') x $gaplen;
	my $ncut = 0;
	for (my $gapb = index($REF->{$ref}, $gap); $gapb >= 0; ) {
		my $reflen = length($REF->{$ref});
		my $gape = $gapb + $gaplen;
		$gape++ while ($gape < $reflen && substr($REF->{$ref}, $gape, 1) eq "N"); 
		my $gapw = $gape - $gapb;
		my $keep = 6 + ($gapw % 3);
		my $cut = $gapw - $keep;
		$ncut += $cut;
		my $facut = substr($REF->{$ref}, 0, $gapb).substr("NNNNNNNNN", 0, $keep).substr($REF->{$ref}, $gape);
		$REF->{$ref} = $facut;
		$gapb = index($REF->{$ref}, $gap);
	}
	return $ncut;
}

sub print_log {
	my ($opt, $ref, $REF, $reflen, $all_matches, $nvec, $ncontam, $nmasked, $ntrim, $ncut, $logfh) = @_;
	my $npercent = $opt->{npercent}/100;
	my $new_reflen = length($REF->{$ref});
	my $new_nN = $REF->{$ref} =~ tr/N/N/;
	my $log = $new_reflen == $reflen ? $reflen : "$new_reflen; olen=$reflen; ncut=$ncut; nnn=$new_nN/$nmasked;";
	my @vlog;
	push(@vlog, "vecmask=$nvec") if ($nvec);
	push(@vlog, "contmask=$ncontam") if ($ncontam);
	push(@vlog, "[".join(';', @$all_matches)."]") if (@$all_matches);
	my $vlog = join(' ', @vlog);
	my $tlog = $ntrim ? "endtrim=$ntrim;" : '';
	if ($new_reflen < $opt->{minsize}) {
		print $logfh "#$ref\tlen=$log\t$tlog\t$vlog\tERR: too short:$new_reflen\n";
		delete($REF->{$ref});
	}
	elsif ($new_reflen*$npercent < $new_nN) {
		print $logfh "#$ref\tlen=$log\t$tlog\t$vlog\tERR: too noisy:${new_nN}n\n";
		delete($REF->{$ref});
	}
	else {
		print $logfh ">$ref\tlen=$log\t$tlog\t$vlog\n";
	}
}

sub get_corrections {
	my $report = shift;
	open(REPORT, $report) || die "Can't open file $report\n";
	my ($trim, $exclude, $desc);
	my $corrections;
	while (<REPORT>) {
		next unless ($desc || /^(Exclude|Trim):/);
		next if (/^Sequence name, length, /);
		next if (/^$/);
		if (/^Exclude:/) {
			$trim = 0;
			$exclude = 1;
			$desc = 1;
			next;
		}
		if (/^Trim:/) {
			$trim = 1;
			$exclude = 0;
			$desc = 1;
			next;
		}
		chomp;
		my @line = split("\t", $_);
		die "Unknown FCSreport format\n" unless (($#line == 2 && $exclude) || ($#line == 3 && $trim));
		$corrections .= join("\t", $line[0], $line[2], $line[3], $line[3] =~ /vector/ ? 'V' : 'M')."\n" if ($trim);
		$corrections .= join("\t", $line[0], '1..'.$line[1], $line[2], $line[2] =~ /vector/ ? 'V' : 'M')."\n" if ($exclude);
	}
	close REPORT;
	return $corrections||'';
}

sub get_contaminant_match {
	my ($ref, $REF, $rrnafh) = @_;
	my ($contaminant, $match, $suspect);
	my $contaminant_out = run_contaminant($ref, $REF, $rrnafh);
	foreach my $contaminant_db (@{$opt->{contaminant}->{db}}) {
		for (split("\n", $contaminant_out->{$contaminant_db})) {
			chomp;
			next if (/^#/);
			my @hit = split("\t", $_);
			$contaminant .= "$hit[6]\t$hit[7]\t$hit[1]\tM\n" if (&{$opt->{contaminant}->{valid_hit}->{$contaminant_db}}(@hit)); # match type is M
		}
	}
	return $contaminant||'';
}

sub get_vecscreen_match {
	my ($ref, $REF) = @_;
	my ($vecscreen, $match, $suspect);
	my @vecscreen_out = split("\n", run_vecscreen($ref, $REF));
	for (@vecscreen_out) {
		chomp;
		if (/match/) {
			$match = $_;
			$suspect = 0;
		}
		elsif (/^Suspect origin/) {
			$suspect = 1;
		}
		elsif (/^(\d+)\s+(\d+)/) {
			next if ($suspect);
			$vecscreen .= "$1\t$2\t$match\tV\n" unless ($match =~ /Weak/); # match type is V
		}
	}
	return $vecscreen||'';
}

sub run_contaminant {
	my ($ref, $REF, $rrnafh) = @_;
	my $input = ">$ref\n$REF->{$ref}\n";
	my ($return, $output);
	foreach my $contaminant_db (@{$opt->{contaminant}->{db}}) {
		$output = run_process($input, \@{$opt->{contaminant}->{blastn_cmd_line}->{$contaminant_db}});
		if (defined($rrnafh) && $contaminant_db eq 'rrna') { print $rrnafh $output };
		$return->{$contaminant_db} = $output;
	}
	return $return;
}

sub run_vecscreen {
	my ($ref, $REF) = @_;
	my $input = ">$ref\n$REF->{$ref}\n";
	return run_process($input, \@{$opt->{vecscreen_cmd_line}});
}

sub run_process {
	my ($input, $cmd) = @_;
	my ($output, $err);
	my $timeout = 30;
	my $try = 1;
	RUN: eval { run $cmd, \$input, \$output, \$err, timeout($timeout) };
	if ($@ =~ /timeout/) {
		$timeout *= 2;
		goto RUN unless ($timeout > 120);
		die("Timeout exceeded\nReference:\n$input\nCommand:\n".join(' ', @$cmd)."\n");
	} elsif ($err =~ /ERROR/) {
		$try++;
		goto RUN unless ($try > 3);
		die("Error: $err\nReference:\n$input\nCommand:\n".join(' ', @$cmd)."\n");
	} elsif ($@) {
		$try++;
		goto RUN unless ($try > 3);
		die("Unable to run command through IPC::RUN: $@\nReference:\n$input\nCommand:\n".join(' ', @$cmd)."\n");
	}
	return $output;
}

#============================================================

=head1 NAME

fasta_tsa.pl - clean reference sequences for NCBI TSA submission

=head1 SYNOPSIS

fasta_tsa.pl [options] --db adaptors_for_screening_euks.fa refseq.fa

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<-d, --db>

BLAST database name to detect vector contaminant using vecscreen.

Resource available at:
 ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa

=item B<-c, --contaminant-dir>

Directory containing other BLAST databases related to NCBI's Foreign Contamination Screens:
 -contam_in_euks.fa: contaminants database that contains vector sequences, bacterial insertion 
sequences, E. coli and phage genomes.
 -mito.nt: a database of mitochondrial genomes.

Resources available at:
 ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz
 ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/mito.nt.gz

=item B<-p, --npercent>

Set the n's rate required to discard a sequence [10].

=item B<-g, --ngap>

Set the n's length in a row required to compress a gap [10].

=item B<-e, --endtrim>

Trim sequences containing at least --ntrim Ns within the last --endtrim bases [20].

=item B<-t, --ntrim>

Trim sequences containing at least --ntrim Ns within the last --endtrim bases [10].

=item B<-s, --minsize>

Set the minimum size required to keep a sequence [200].

=item B<--defline>

Reformat the FASTA definition line for TSA submission.

=item B<--rrna>

Print the output of blastn to rrna contaminant db to this file.
This file is needed for the execution of the NCBI Foreign Contamination 
Screen workflow.

=item B<-l, --log>

Print log messages to this file.

=back

=head1 DESCRIPTION

Clean reference sequences for NCBI TSA submission. Mask Strong or Moderate vecscreen matches. Mask potential contamination matches if --db-contaminant option is provided. Compress gaps greater than 14 n's in a row. Trim ends containing at least 10 Ns within the last 20 bases. Discard sequences containing more than 10% n's or shorter than 200 bases. Add the defline option to format the FASTA definition line: check that unique identifiers do not exceed 50 characters and add the [moltype=transcribed_RNA] component.
Print on STDOUT the corrected reference sequences. If log file not provided, print log messages on STDERR. Inspired from http://arthropods.eugenes.org/genes2/evigene/scripts/rnaseq/asmrna2ncbitsa.pl.

=head1 DEPENDANCY

NCBI BLAST vecscreen from NCBI tools version 6.1 or later.
NCBI BLAST blastn from NCBI Blast+ package.

=head1 AUTHORS

Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=cut
