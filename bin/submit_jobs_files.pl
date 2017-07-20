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
use POSIX qw(ceil);
use File::Basename;
use AsmUtils;

my $dir = shift;
my $opt = get_drap_config($dir) or croak("Unable to read drap config file .drap_conf.json inside $dir");

(my $name = basename($dir)) =~ s/_(.)/uc($1)/eg;
my $n = 0;
my $job_id;
my $process_id = $$;
foreach my $step (@{$opt->{steps}}) {
	$n++;
	my $scripts = $opt->{scripts}->{$step};
	foreach my $i (0..$#$scripts) {
		my $script = $scripts->[$i];
		next unless (-e $script);
		printf("Submit %ssteps:\n", $opt->{restart} && -e $opt->{submit_log} ? 'again failed ' : '') unless ($job_id || $process_id != $$);
		if ($opt->{local}) {
			print "Execute $script\n";
			if ($step eq 'dbg' && $opt->{dbg} eq 'oases') {
				my $task = 1;
				open(SHELL, "$script") or croak "Can't open file $script";
				while (<SHELL>) {
					chomp;
					system(qq(csh -c "( csh -c '$_' > $dir/err_log/j$n-$name.o$process_id.$task ) >& $dir/err_log/j$n-$name.e$process_id.$task"));
					$task++;
				}
				close SHELL;
			} else {
				system(qq(csh -c "( csh $script > $dir/err_log/j$n-$name.o$process_id ) >& $dir/err_log/j$n-$name.e$process_id"));
			}
			$process_id++;
		} else {
 			my ($reservation, $array, $submit_cmd, $submit_msg);
			if ($step eq 'dbg') {
				$array = 1 if ($opt->{dbg} eq 'oases');
				my $index = $#$scripts > 0 ? sprintf("_%d",$i+1) : '';
				$reservation = $opt->{env}->{"$opt->{dbg}${index}_res"}||''; # get reservation from drap.cfg for specific step
				my $mem = $opt->{dbg} eq 'oases' ? $opt->{env}->{oases_ram} : ceil($opt->{dbg_mem}/$opt->{env}->{n_cpu});
				$reservation =~ s/\%MEM\%/$mem/;
				my $h_vmem = ceil($mem*3/2);
				$reservation =~ s/\%H_VMEM\%/$h_vmem/;
			} elsif ($step eq 'preprocess' && !$opt->{no_norm}) {
				$reservation = $opt->{env}->{normalize_res}||''; # get reservation from drap.cfg for specific step
				my $mem = ceil($opt->{norm_mem}/$opt->{env}->{n_cpu});
				$reservation =~ s/\%MEM\%/$mem/;
				my $h_vmem = ceil($mem*3/2);
				$reservation =~ s/\%H_VMEM\%/$h_vmem/;
			} elsif ($step =~ /postprocess/ && $opt->{no_rate}) {
				$reservation = '';
			} else {
				$reservation = $opt->{env}->{"${step}_res"}||''; # get reservation from drap.cfg for specific step
			}
			my $options = sprintf("%s%s%s", $reservation ? $reservation.' ' : '', $job_id ? "-hold_jid $job_id " : '', $opt->{debug} ? '-verify ' : '',);
			my $cmd = sprintf("%s/submitJob --print --complete %s--name j%d-%s --stdout %s/err_log --stderr %s/err_log --log %s %s-- %s",
				$opt->{binpath}, $array ? '--array ' : '', $n, $name, $dir, $dir, $opt->{submit_log}, $options ? "--options $options " : '', $script
			);
			($submit_cmd, $submit_msg, $job_id) = split("\n", process_cmd(1, $cmd));
			print "$submit_cmd\n$submit_msg\n";
		}
	}
}
exit 0;
