# $Id$

=pod
=head1 NAME

AnalysisReport - Report workflow group of step.

=head1 DESCRIPTION

 The analysis is a group of steps in a workflow.
 See description in package Report.

=head1 AUTHORS

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package AnalysisReport ;
use strict ;
use JSON ;
use StepReport ;


our $_current_analysis_id = 0 ;

=head2 function new

 Usage        : $analysis = new AnalysisReport( $title[, $description] )
 Function     : Create, and return a new analysis.
 Returns      : [AnalysisReport] The analysis.
 Args         : [str] The analysis title.
                [str] The analysis description.

=cut
sub new {
	my ($class, $title, $description) = @_ ;

	my @steps = ();
	my $self = {
		'id'    => $_current_analysis_id,
		'title' => $title,
		'description' => $description,
		'steps' => \@steps
	};
	$_current_analysis_id++ ;
	bless( $self, $class );

	return $self ;
}

=head2 function add_step

 Usage        : $step = $analysis->add_step( $title[, $description] )
 Function     : Create, Add to analysis and return a new step.
 Returns      : [StepReport] The new analysis.
 Args         : [str] The step title.
                [str] The step description.

=cut	
sub add_step {
	my ($self, $title, $description) = @_ ;
	my $new_step = undef ;

	eval {
		my $step = $self->get_step_by_title( $title );
		die "The step '".$title."' in analysis '".$self->{'title'}."' already exists." ;
		1;
	} or do { # The analysis doen't exist
		$new_step = new StepReport( $title, $description );
		push( @{$self->{'steps'}}, $new_step );
	};

	return $new_step ;
}

=head2 function get_nb_steps

 Usage        : $nb = $analysis->get_nb_steps()
 Function     : Return the number of steps in the analysis.
 Returns      : [int] The number of steps.
 Args         : None

=cut
sub get_nb_steps {
	my ($self) = @_ ;
	return scalar( @{$self->{'steps'}} );		
}

=head2 function get_step_by_title

 Usage        : $step = $analysis->get_step_by_title( $title )
 Function     : Return the step with $title.
 Returns      : [StepReport] The step.
 Args         : [str] The step title.

=cut
sub get_step_by_title {
	my ($self, $title) = @_ ;
	my $find = undef ;
	my $idx = 0 ;
	my $nb_steps = $self->get_nb_steps() ;
	while( !defined($find) && $idx < $nb_steps ){
		if( ${$self->{'steps'}}[$idx]->{'title'} eq $title ){
			$find = ${$self->{'steps'}}[$idx] ;
		}
		$idx++ ;
	}
	if( !defined($find) ){
		die "the step '".$title."' doesn't exist in analysis '".$self->{'title'}."'." ;
	}
	
	return $find ;
}

=head2 function get_or_create_step

 Usage        : $step = $analysis->get_or_create_step( $title, [$description] )
 Function     : Returns the step if it exist. Otherwise creates and returns a new step.
 Returns      : [StepReport] The step.
 Args         : [str] The step title.
                [str] The step description.

=cut
sub get_or_create_step {
	my ($self, $title, $description) = @_ ;
	
	my $step = undef ;
	eval {
		$step = $self->get_step_by_title( $title );
		1;
	} or do {
		$step = $self->add_step( $title, $description );
	};

	return $step ;
}

=head2 function get_hash

 Usage        : %hash_representation = $analysis->get_hash()
 Function     : Return an hash representation of the analysis report.
 Returns      : [hash] The hash representation.
 Args         : None

=cut
sub get_hash {
	my ($self) = @_ ;

	my @steps = ();
	my %hash = (
		'id'    => $self->{id},
		'title' => $self->{title},
		'description' => $self->{description},
		'steps' => \@steps
	);
	foreach my $current_step ( @{$self->{steps}} ){
		push( @{$hash{'steps'}}, $current_step->get_hash() );
	}

	return \%hash ;
}

=head2 procedure from_hash

 Usage        : $analysis->from_hash(%hash_representation)
 Function     : Set the analysis object from an hash representation of the 
                analysis.
 Returns      : [str] The hash representation.
 Args         : None

=cut
sub from_hash {
	my ($self, $hash) = @_ ;
	
	$self->{id} = $$hash{'id'} ;
	$self->{title} = $$hash{'title'} ;
	$self->{description} = $$hash{'description'} ;
	foreach my $current_step ( @{$$hash{'steps'}} ){
		my $new_step = new StepReport();
		$new_step->from_hash( $current_step );
		push( @{$self->{steps}}, $new_step );
	}
}

=head2 function get_json

 Usage        : $json_representation = $analysis->get_json()
 Function     : Return a json string representation of the report.
 Returns      : [str] The json representation.
 Args         : None

=cut
sub get_json {
	my ($self) = @_ ;
	return to_json( $self->get_hash(), { utf8 => 1, pretty => 1 } );
}


1;