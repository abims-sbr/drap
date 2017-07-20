#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

#------------------------------------------------------------
# options and parameters
#------------------------------------------------------------
# Options definitions (name,format,default,required)
my @Options = ( ['help',undef,undef,0],
	['man',undef,undef,0],
	['m|mindepth','i',10,undef],
	['s|substitution-only',undef,undef,0],
	['i|indel-only',undef,undef,0],
	['l|log','s',undef,0]
);
# defaults values
my %Values=();
my %Required=();
map {$Values{$_->[0]}=$_->[2];$Required{$_->[0]}=$_->[3];} @Options;
# build options list
my %Options=();
map {$Options{defined($_->[1])?$_->[0].'='.$_->[1]:$_->[0]}=\$Values{$_->[0]};} @Options;
# retrieve options
GetOptions(%Options) || pod2usage(2);
pod2usage(1) if ($Values{'help'});
pod2usage(-exitstatus => 0, -verbose => 2) if ($Values{'man'});

# check required
map {pod2usage("$_ is required") if ($Required{$_}&&!defined($Values{$_}));} keys %Values;

# read input parameters
$#ARGV==0 || pod2usage("missing reference file as parameter");

# log file
my $logfh= *STDERR;
if($Values{'l|log'}) { open(LOG, ">$Values{'l|log'}") or die "Can't open file $Values{'l|log'}"; $logfh= *LOG; }

#------------------------------------------------------------
# read reference sequences and store it into a hash
#------------------------------------------------------------
my %REF=();
open(REF,$ARGV[0]) || die "Can't open file $ARGV[0]";

my $ref='';
while(my $line=<REF>) {
	chomp $line;
	if($line=~/^>(\S+)/){
		if(!exists($REF{$1})){
			$ref=$1;
			$REF{$ref}='';
		}
		else {
			die "Sequence $1 already seen";
		}
	}
	else {
		$line=~s/\s//g;
		$REF{$ref}.=uc($line);
	}
}
close(REF);

#------------------------------------------------------------
# read mpileup output
#------------------------------------------------------------
$ref='';
my $refseq='';
my %variation;

while (my $line=<STDIN>) {
	chomp $line;
	my @F=split/\s+/,$line;
	if ($F[0] ne $ref) { # next reference
		if ($ref ne '') { # correct reference
			correct_variation(\%variation, \%REF, $ref, $logfh) if (%variation);
		}
		$ref=$F[0];
		exists($REF{$ref})||die "unknown reference $ref";
		%variation = ();
	}
	next unless $F[3] >= $Values{'m|mindepth'};
	my (%currentIndel, %currentSubst, $type);
	# insertion pattern: \+[0-9]+[ACGTNacgtn]+ ; deletion pattern: -[0-9]+[ACGTNacgtn]+
	while ($F[4] =~ /([+-])(\d+)/) {
		$F[4] =~ s/([$1])$2([ACTGNactgn]{$2})//; # step needed even in mode substitution-only to make substitution collect regexp efficient
		$currentIndel{$1.uc($2)}++ unless $Values{'s|substitution-only'}; # collect indels
	}
	unless ($Values{'i|indel-only'}) {
		while ($F[4] =~ /([ACTGNactgn])/g) { $currentSubst{uc($1)}++ }; # collect substitutions
	}
	my $majIndel = shift(@{[sort { $currentIndel{$b} <=> $currentIndel{$a} } keys %currentIndel]});
	my $majSubst = shift(@{[sort { $currentSubst{$b} <=> $currentSubst{$a} } keys %currentSubst]});
	if ($majIndel && $currentIndel{$majIndel} > $F[3]/2) {
		$variation{$F[1]}->{indel} = $majIndel;
		$variation{$F[1]}->{depth} = $F[3];
		$variation{$F[1]}->{indel_evidence} = $currentIndel{$majIndel};
	}
	if ($majSubst && $currentSubst{$majSubst} > $F[3]/2) {
		$variation{$F[1]}->{subst} = $majSubst;
		$variation{$F[1]}->{depth} = $F[3];
		$variation{$F[1]}->{subst_evidence} = $currentSubst{$majSubst};
	}
}

if ($ref ne '') { # correct last reference
	correct_variation(\%variation, \%REF, $ref, $logfh) if (%variation);
}
close $logfh;

#------------------------------------------------------------
# print reference sequences
#------------------------------------------------------------
foreach $ref (sort keys %REF) {
	$REF{$ref}=~s/(.{60})/$1\n/g;
	chomp $REF{$ref};
	print ">$ref\n$REF{$ref}\n";
}

sub correct_variation {
	my ($variation, $REF, $ref, $logfh) = @_;
	my $action;
	my $index;
	foreach my $pos (sort { $b <=> $a } keys %$variation) {
		my $refseq=substr($REF->{$ref}, 0, $pos);
		if (exists $variation->{$pos}->{subst}) {
			$refseq=~s/(.)$/$variation->{$pos}->{subst}/e;
			$index=$pos;
			printf $logfh ("%s: replace %s by %s at position %d (%d/%d)\n",$ref, $1, $variation->{$pos}->{subst}, $pos, $variation->{$pos}->{subst_evidence}, $variation->{$pos}->{depth});
		}
		if (exists $variation->{$pos}->{indel}) {
			if ($variation->{$pos}->{indel} =~ /\+(.+)/) {
				$action='add';
				$refseq.=$1;
				$index=$pos;
			}
			elsif ($variation->{$pos}->{indel} =~ /-(.+)/) {
				$action='remove';
				$index=$pos+length($1);
			}
			printf $logfh ("%s: %s %s at position %d (%d/%d)\n",$ref, $action, $1, $pos+1, $variation->{$pos}->{indel_evidence}, $variation->{$pos}->{depth});
		}
		$refseq.=substr($REF->{$ref}, $index);
		$REF->{$ref}=$refseq;
	}
}


#============================================================

=head1 NAME

samCorrectVariation.pl - correct indels and/or substitutions in reference sequences from evidences seen in mpileup output

=head1 SYNOPSIS

samCorrectVariation.pl [options] refseq.fa < mpileup.out

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<-m, --mindepth>

Set the minimum depth required to engage in a correction (default 10).

=item B<-s, --substitution-only>

Correct only substitutions.

=item B<-i, --indel-only>

Correct only insertions or deletions.

=item B<-l, --log>

Print log messages to this file.

=back

=head1 DESCRIPTION

Collect insertions, deletions and substitutions at each position of the reference sequence.
Correct reference sequence to follow the majority vote at each position of the alignment if mindepth is reached.
Print on STDOUT the corrected reference sequences. If log file not provided, print log messages on STDERR.

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut
