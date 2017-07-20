# $Id$

=pod
=head1 NAME

Report - Report workflow.

=head1 DESCRIPTION

 The report is used to store metrics during workflow execution.
 The report is composed of analysis and an analysis can be subdivide in 
 steps. The metrics on data along the workflow are stored by step.
 Example :
    Report: alignment
 	   - analysis: cleanning
 		   - step: filter on N
 		   - step: clean contaminant
 	   - analysis: alignment
 	       - step: alignment

=head1 AUTHOR

 Frédéric Escudié - Plateforme Bioinformatique de Toulouse
 (support.genopole@toulouse.inra.fr)

=head1 COPYRIGHT

 2015 INRA

=head1 LICENSE

 GNU GPLv3

=cut

package Report ;
use strict ;
use JSON ;
use AnalysisReport ;


=head2 function new

 Usage        : $report = new Report( [$cmd_line] )
 Function     : Creates, and returns a new report.
 Returns      : [Report] The report.
 Args         : [str] The workflow command line.

=cut
sub new {
	my ($class, $cmd_line) = @_ ;

	my @analysis = ();
	my $self = {
		'processed' => $cmd_line,
		'analysis'  => \@analysis
	};
	bless( $self, $class );

	return $self ;
}

=head2 function add_analysis

 Usage        : $analysis = $report->add_analysis( $title[, $description] )
 Function     : Creates, Adds to report and returns a new analysis.
 Returns      : [AnalysisReport] The new analysis.
 Args         : [str] The analysis title.
                [str] The analysis description.

=cut
sub add_analysis {
	my ($self, $title, $description) = @_ ;
	my $new_analysis = undef ;

	eval {
		my $analysis = $self->get_analysis_by_title( $title );
		die "The analysis '".$title."' already exists." ;
		1;
	} or do { # The analysis doen't exist
		$new_analysis = new AnalysisReport( $title, $description );
		push( @{$self->{analysis}}, $new_analysis );
	};

	return $new_analysis ;
}

=head2 function get_nb_analysis

 Usage        : $nb = $report->get_nb_analysis()
 Function     : Returns the number of analysis in the report.
 Returns      : [int] The number of analysis.
 Args         : None

=cut
sub get_nb_analysis {
	my ($self) = @_ ;
	return scalar( @{$self->{'analysis'}} );		
}

=head2 function get_analysis_by_title

 Usage        : $analysis = $report->get_analysis_by_title( $title )
 Function     : Returns the analysis with $title.
 Returns      : [AnalysisReport] The analysis.
 Args         : [str] The analysis title.

=cut
sub get_analysis_by_title {
	my ($self, $title) = @_ ;
	my $find = undef ;
	my $idx = 0 ;
	my $nb_analysis = $self->get_nb_analysis() ;
	while( !defined($find) && $idx < $nb_analysis ){
		if( ${$self->{'analysis'}}[$idx]->{'title'} eq $title ){
			$find = ${$self->{'analysis'}}[$idx] ;
		}
		$idx++ ;
	}
	if( !defined($find) ){
		die "The analysis '".$title."' doesn't exist." ;
	}

	return $find ;
}

=head2 function get_or_create_analysis

 Usage        : $analysis = $report->get_or_create_analysis( $title, [$description] )
 Function     : Returns the analysis if it exist. Otherwise creates and returns a new analysis.
 Returns      : [AnalysisReport] The analysis.
 Args         : [str] The analysis title.
                [str] The analysis description.

=cut
sub get_or_create_analysis {
	my ($self, $title, $description) = @_ ;
	
	my $analysis = undef ;
	eval {
		$analysis = $self->get_analysis_by_title( $title );
		1;
	} or do {
		$analysis = $self->add_analysis( $title, $description );
	};

	return $analysis ;
}

=head2 function get_hash

 Usage        : %hash_representation = $report->get_hash()
 Function     : Return an hash representation of the report.
 Returns      : [hash] The hash representation.
 Args         : None

=cut
sub get_hash {
	my ($self) = @_ ;

	my @analysis = ();
	my %hash = (
		'processed' => $self->{'processed'},
		'analysis'  => \@analysis
	);
	foreach my $current_analysis ( @{$self->{analysis}} ){
		push( @{$hash{'analysis'}}, $current_analysis->get_hash() );
	}

	return \%hash ;
}

=head2 procedure from_hash

 Usage        : $report->from_hash(%hash_representation)
 Function     : Set the report object from an hash representation of the report
                (see get_hash and load_json).
 Returns      : [str] The hash representation.
 Args         : None

=cut
sub from_hash {
	my ($self, $hash) = @_ ;
	
	$self->{'processed'} = $$hash{'processed'} ;
	foreach my $current_analysis ( @{$$hash{'analysis'}} ){
		my $new_analysis = new AnalysisReport();
		$new_analysis->from_hash( $current_analysis );
		push( @{$self->{analysis}}, $new_analysis );
	}
}

=head2 function get_json

 Usage        : $json_representation = $report->get_json()
 Function     : Return a json string representation of the report.
 Returns      : [str] The json representation.
 Args         : None

=cut
sub get_json {
	my ($self) = @_ ;
	return to_json( $self->get_hash(), { utf8 => 1, pretty => 1 } );
}

=head2 function load_json

 Usage        : $report = Report::load_json( $path )
 Function     : Return a report from a report json file.
 Returns      : [Report] The report.
 Args         : [str] The path to the json.

=cut
sub load_json {
	my ($path) = @_ ;

	# Load data
	local $/ = undef ;
	open( my $FH_json, '<', $path ) or die "Unable to open '".$path."'" ;
	my $json_text = <$FH_json> ;
	close( $FH_json );
	my $hash = from_json( $json_text, { utf8  => 1 } );

	# Set report
	my $report = new Report();
	$report->from_hash( $hash );
	
	return $report ;
}

1;
