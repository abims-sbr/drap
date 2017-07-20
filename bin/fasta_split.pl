#!/usr/bin/perl -w

=head1 Description

  
  Read a fasta file with multiple entries (parameter 1)
  Split this file into parts of maximum X (given as parameter 2)
     - Kb (parameter 3 = Kb)
     - entries (parameter 3 = entries)
  If fasta file is toto.tfa, parts will be named toto_0001.tfa ... toto_xxxx.tfa

=head1 AUTHORS

 Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

BEGIN
{
  ($prg)=($0=~/([^\/]+)$/);
  $dir=$0;
  $dir=~s/$prg$//;
  push @INC, $dir;
}

#------------------------------------------------------------

sub usage( $ )
  {
    printf STDERR "%s\n", $_[0];
    system("pod2text $0");
    exit(1);
  }

#------------------------------------------------------------

($#ARGV == 2) || &usage("bad parameter");

($ARGV[1] =~ /^\d+\.?\d*$/) || &usage("parameter is not a number: $ARGV[1]");

($ARGV[2] eq 'Kb' || $ARGV[2] eq 'entries') || &usage("bad size type $ARGV[2]");

if ($ARGV[2] eq 'Kb')
  {
    $max_length = $ARGV[1]*1000;
  }
else
  {
    $max_entries = $ARGV[1];
  }

$current_length = 0;
$nb_entries = 0;
if ($ARGV[0] =~ /\|\s*$/ || $ARGV[0] eq '-')
  {
	$basename = "split$$";
  }
else
  {
	$basename = $ARGV[0];
	$basename =~ s/\.\w+$//g;
	($basename) = ($basename =~ /([^\/]+)$/);
  }
$part_no = 1;

$inseq=0;

# try to open input file
open(FIN,$ARGV[0]) || open(FIN,"cat $ARGV[0]|") || die "Can't open $ARGV[0]";

# open output file
open(FOUT,">${basename}_0001.tfa") || die "can't create file ${basename}_0001.tfa";

while (<FIN>)
  {
    chomp;
    
    if (/^>/)
      {
	if ($inseq == 1)
	  {
	    if (($ARGV[2] eq 'Kb' && $current_length > $max_length) ||
		($ARGV[2] eq 'entries' && $nb_entries == $max_entries))
	      {
		close(FOUT);
		$part_no++;
		$current_length = 0;
		$nb_entries = 0;
		$part_name = sprintf("%s_%04d.tfa",$basename, $part_no);
		open(FOUT,">$part_name") || die "Can't create file $part_name";
	      }
	  }
	$inseq=1;
	$nb_entries++;
	    
	print FOUT "$_\n";
	next;
      }
    if ($inseq == 1)
    {
      print FOUT "$_\n";
      $current_length += length($_);
    }
  }

close(FOUT);
close(FIN);
