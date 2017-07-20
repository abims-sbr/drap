# $Id$

=pod
=head1 NAME

SchedulerFactory - Scheduler factory

=head1 DESCRIPTION

 This factory return the appropriate scheduler object.
 By default it returns a local schedule or the type of scheduler stored with 
 key type in section SCHEDULER in configuration file.

=head1 VERSION

 1.0.0

=head1 AUTHOR

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package SchedulerFactory ;
use strict ;
use warnings ;
use File::Basename ;
use ConfigFile ;
use SGEScheduler ;
use LocalScheduler ;


=head2 function instantiate

 Usage        : $scheduler = SchedulerFactory->instantiate($work_dir, $type, $scheduler_cpu)
 Function     : Returns the appropriate scheduler object.
 Returns      : [AbstractScheduler] An specialised class from the AbstractScheduler.
 Args         : [str] The work directory for scheduler.
                [str] The scheduler type ("local" or "sge").
                [str] The number of available cpu for local scheduler.
=cut
sub instantiate {
	my ($class, $working_dir, $scheduler_type, $scheduler_cpu) = @_ ;
	my $config_path = $main::config_file;
	
	# Select the scheduler type
	if( !defined($scheduler_type) ){
		$scheduler_type = "sge" ;
		if( $config_path && -e $config_path && -r $config_path ){ # If a config file exists
			my $config = new ConfigFile( $config_path );
			if( defined($config->{'SCHEDULER'}) && defined($config->{'SCHEDULER'}{'type'}) ){
				$scheduler_type = $config->{'SCHEDULER'}{'type'} ;
 			}
		}
	}
	# Select nb cpu for local scheduler
	if( $scheduler_type eq "local" && !defined($scheduler_cpu) ) {
		if( $config_path && -e $config_path && -r $config_path ){
			my $config = new ConfigFile( $config_path );
			if( defined($config->{'SCHEDULER'}) && defined($config->{'SCHEDULER'}{'local_cpu'}) ){
				$scheduler_cpu = $config->{'SCHEDULER'}{'local_cpu'} ;
			}
		}
	}
	
	# Create the scheduler
	my $scheduler = "" ;
	if( $scheduler_type eq "local" ){
		$scheduler = new LocalScheduler($working_dir, $scheduler_cpu) ;
	} elsif( $scheduler_type eq "sge" ){
		$scheduler = new SGEScheduler($working_dir) ;
	} else {
		die "Unknown scheduler '".$scheduler_type."'" ;
	}
	
	return $scheduler ;
}


1;