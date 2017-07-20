# $Id$

=pod
=head1 NAME

StepReport - Report workflow step.

=head1 DESCRIPTION

 The report store metrics on data.
 See description in package Report.

=head1 AUTHOR

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package StepReport ;
use strict ;
use JSON ;


our $_current_step_id = 0 ;
our $_current_metrics_id = 0 ;

=head2 function new

 Usage        : $step = new StepReport( $title[, $description] )
 Function     : Create, and return a new step in analysis.
 Returns      : [StepReport] The step.
 Args         : [str] The step title.
                [str] The step description.

=cut
sub new {
	my ($class, $title, $description) = @_ ;
 
	my %metrics = ();
	my $self = {
		'id'      => $_current_step_id,
		'title'   => $title,
		'description' => defined($description) ? $description : "",
		'metrics' => \%metrics
	};
 	$_current_step_id++ ;
	bless( $self, $class );
	
	return $self ;
}

=head2 function add_metrics

 Usage        : $metrics_id = $step->add_metrics($template[, $sample] )
 Function     : Create a new template of metrics in step.
 Returns      : [int] The metrics ID.
 Args         : [str] The template of metrics ('fastaMetrics', 'alnMetrics', 
                      ...)
                [str] The sample concerned by my metrics.
=cut
sub add_metrics {
	my ($self, $template, $sample) = @_ ;
	$sample = defined($sample) ? $sample : "default" ;
	$self->{metrics}{$template}{$sample} = $_current_metrics_id ;
	$_current_metrics_id++ ;
	return $self->{metrics}{$template}{$sample};
}

=head2 function get_or_create_metrics_filename

 Usage        : $filename = $step->get_or_create_metrics_filename( $template, [$sample] )
 Function     : Returns the metrics filename if the metrics exist. Otherwise 
                creates metrics and its filename.
 Returns      : [str] The metrics filename.
 Args         : [str] The metrics template ('fastaMetrics', 'alnMetrics', ...).
                [str] The sample concerned by metrics.

=cut
sub get_or_create_metrics_filename {
	my ($self, $template, $sample) = @_ ;
	return "metrics_".$self->get_or_create_metrics($template, $sample).".json" ;
}

sub get_metrics_id {
	my ($self, $template, $sample) = @_ ;
	$sample = defined($sample) ? $sample : "default" ;
	my $metrics_id = undef ;
	if( exists($self->{'metrics'}{$template}) && exists($self->{'metrics'}{$template}{$sample}) ){
		$metrics_id = $self->{'metrics'}{$template}{$sample} ;
	} else {
		die "The metrics '".$template."' for sample '".$sample."' doesn't exist in step '".$self->{'title'}."'." ;
	};
	
	return $metrics_id ;
}

=head2 function get_or_create_metrics

 Usage        : $metrics_id = $step->get_or_create_metrics( $template, [$sample] )
 Function     : Returns the metrics ID if it exist. Otherwise creates metrics 
                and its ID.
 Returns      : [int] The metrics ID.
 Args         : [str] The metrics template ('fastaMetrics', 'alnMetrics', ...).
                [str] The sample concerned by metrics.

=cut
sub get_or_create_metrics {
	my ($self, $template, $sample) = @_ ;
	
	my $metrics_id = undef ;
	eval {
		$metrics_id = $self->get_metrics_id( $template, $sample );
		1;
	} or do {
		$metrics_id = $self->add_metrics( $template, $sample );
	};

	return $metrics_id ;
}

=head2 function get_hash

 Usage        : %hash_representation = $step->get_hash()
 Function     : Return an hash representation of the step report.
 Returns      : [hash] The hash representation.
 Args         : None

=cut	
sub get_hash {
	my ($self) = @_ ;
	
	my %hash = (
		'id'      => $self->{id},
		'title'   => $self->{title},
		'description' => $self->{description},
		'metrics' => $self->{metrics}
	);
	
	return \%hash ;
}

=head2 procedure from_hash

 Usage        : $step->from_hash(%hash_representation)
 Function     : Set the step object from an hash representation of the step.
 Returns      : [str] The hash representation.
 Args         : None

=cut
sub from_hash {
	my ($self, $hash) = @_ ;
	
	$self->{'id'} = $$hash{'id'} ;
	$self->{'title'} = $$hash{'title'} ;
	$self->{'description'} = $$hash{'description'} ;
	foreach my $template ( keys %{$$hash{'metrics'}} ){
		foreach my $sample ( keys %{$$hash{'metrics'}{$template}} ){
			$self->{metrics}{$template}{$sample} = $$hash{'metrics'}{$template}{$sample} ;
			$_current_metrics_id++ ;
		}
	}
}

=head2 function get_json

 Usage        : $json_representation = $step->get_json()
 Function     : Return a json string representation of the step.
 Returns      : [str] The json representation.
 Args         : None

=cut
sub get_json {
	my ($self) = @_ ;
	return to_json( $self->get_hash(), { utf8 => 1, pretty => 1 } );
}


1;