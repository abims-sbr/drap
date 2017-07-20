/**
 * Copyright (C) 2015 INRA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * Returns the HTML table representation of remove chimera metrics.
 * @param step_data {Hash} By sample the remove chimera metrics.
 * @return {String} The HTML table.
 */
var chimera_table = function( step_data ){
	var categories = ["Sample", "Nb seq removed", "Nb seq trimmed", "Nb nt lost"] ;
	var data = new Array();
	for( var sample in step_data ){
		data.push( 
			new Array(
				sample == "default" ? "-" : sample,
				step_data[sample].nb_chimera_removed,
				step_data[sample].nb_chimera_trim,
				step_data[sample].nb_nt_lost					
			)
		);
	}
	
	return table( "Cleaning information", categories, data ) ;
}


/**
 * Returns the hash to init an Highchart's barplot from number of non-overlapping proteins by contig.
 * @example
 *     $('#test').highcharts( contigByNbProt_chart(fasta_data, prot_data, ordered_samples) )
 * @param prot_data {Hash} Number of non-overlapping proteins by contig.
 * @param fasta_data {Hash} By sample the length number of sequences in key 'nb_seq'.
 * @param ordered_samples {Array} The names of samples in display order.
 * @return {Hash} Parameter to use in Highchart's constructor.
 */
var contigByNbProt_chart = function( prot_data, fasta_data, ordered_samples ){
	var data = new Array() ;
	if( ordered_samples == null ){
		ordered_samples = Object.keys( prot_data );
	}
	for( var sample_idx = 0 ; sample_idx < ordered_samples.length ; sample_idx++ ){
		var sample = ordered_samples[sample_idx] ;
		var series = {
			'name': sample,
			'data': new Array(),
			'index': sample_idx
		};
		// Add sample data
		var categories = Object.keys(prot_data[sample].nb_contig_by_nb_prot).sort(function(a, b){return a-b});
		for( var idx = 0 ; idx < categories.length ; idx++ ){
			var curr_category = categories[idx] ;
			var contigs_in_category = prot_data[sample].nb_contig_by_nb_prot[curr_category];
			if( fasta_data == null ){
				series['data'].push( contigs_in_category );
			} else {
				var nb_contigs = fasta_data[sample].nb_seq ;
				series['data'].push( Math.round((contigs_in_category/nb_contigs)*10000)/100 );
			}
		}
		data.push( series );
	}
	var y_title = (fasta_data == null) ? "Nb contigs" : "% contigs" ;
	
	return barplot("Nb of non-overlapping proteins by contigs", "Nb protein(s)", y_title, categories, data, 30) ;
}


/**
 * Returns the HTML table representation of correctVariant metrics.
 * @param step_data {Hash} By sample the correctVariant metrics.
 * @return {String} The HTML table.
 */
var correctVariation_table = function( step_data ){
	var categories = ["Sample", "Nb consensus seq with bias", "Insertion case", "Deletion case", "Substitution case", "Reads divergent from consensus", "Consensus position different of all reads"] ;
	var data = new Array() ;
	for( var sample in step_data ){
		data.push(
			new Array(
				sample == "default" ? "-" : sample,
				step_data[sample].corrected_contigs,
				step_data[sample].insertion_case,
				step_data[sample].deletion_case,
				step_data[sample].substitution_case,
				step_data[sample].errors_in_reads,
				step_data[sample].contig_discrodant
			)
		);
	}
	
	return table( "Correction information", categories, data ) ;
}


/**
 * Returns the hash to init an Highchart's boxplot from FPKM distribution. 
 * @example
 *     var boxplot_param = express_chart( step_data, ordered_samples )
 *     $('#test').highcharts( boxplot_param[0], ordered_samples[1] );
 * @param step_data {Hash} FPKM distribution data by sample.
 * @param ordered_samples {Array} The names of samples in display order.
 * @return {Hash} Parameters to use in Highchart's constructor.
 */
var express_chart = function( step_data, ordered_samples ){
	var data = new Array();
	if( ordered_samples == null ){
		ordered_samples = Object.keys( step_data );
	}
	var assemblies = new Array();
	for( var sample_idx = 0 ; sample_idx < ordered_samples.length ; sample_idx++  ){
		var sample = ordered_samples[sample_idx] ;
		assemblies.push( sample );
		data.push(
			new Array(
				parseInt(step_data[sample].fpkm_distribution.min),
				parseInt(step_data[sample].fpkm_distribution.lower_quartile),
				parseInt(step_data[sample].fpkm_distribution.median),
				parseInt(step_data[sample].fpkm_distribution.upper_quartile),
				parseInt(step_data[sample].fpkm_distribution.max)
			)
		);
	}
	if( assemblies.length == 1 && assemblies[0] == "default" ){
		assemblies[0] = "All" ;
	}
	return boxplot("FPKM distribution", null, "FPKM", assemblies, [{'name': "runDrap", 'data':data }]);
}


/**
 * Returns the HTML table representation of the fasta metrics. 
 * @param step_data {Hash} By sample the fasta metrics.
 * @return {String} The HTML table.
 */
var fastaLength_table = function( step_data ){
	var categories = ["Sample", "Nb seq", "N50", "L50", "Lg sum", "Lg min", "Lg lower quartile", "Lg mean", "Lg median", "Lg upper quartile", "Lg max"] ;
	var data = new Array();
	for( var sample in step_data ){
		data.push(
			new Array(
				(sample == "default" ? "-" : sample),
				step_data[sample].nb_seq,
				step_data[sample].n50,
				step_data[sample].l50,
				step_data[sample].length_distribution.sum,
				step_data[sample].length_distribution.min,
				step_data[sample].length_distribution.lower_quartile,
				step_data[sample].length_distribution.mean,
				step_data[sample].length_distribution.median,
				step_data[sample].length_distribution.upper_quartile,
				step_data[sample].length_distribution.max
			)
		);
	}
	
	return table( "Sequences general information", categories, data ) ;
}


/**
 * Returns the hash to init an Highchart's boxplot from contigs length distribution. 
 * @example
 *     var boxplot_param = fastaLength_chart( step_data, ordered_samples )
 *     $('#test').highcharts( boxplot_param[0], ordered_samples[1] );
 * @param step_data {Hash} By sample the quartile length in key length_distribution.
 * @param ordered_samples {Array} The names of samples in display order.
 * @return {Hash} Parameter to use in Highchart's constructor.
 */
var fastaLength_chart = function( step_data, ordered_samples ){
	var data = new Array();
	if( ordered_samples == null ){
		ordered_samples = Object.keys( step_data );
	}
	var assemblies = new Array();
	for( var sample_idx = 0 ; sample_idx < ordered_samples.length ; sample_idx++  ){
		var sample = ordered_samples[sample_idx] ;
		var metrics = step_data[sample] ;
		assemblies.push(sample);
		data.push( new Array(
			parseInt(metrics.length_distribution.min),
			parseInt(metrics.length_distribution.lower_quartile),
			parseInt(metrics.length_distribution.median),
			parseInt(metrics.length_distribution.upper_quartile),
			parseInt(metrics.length_distribution.max)
		));
	}
	
	return boxplot( "Contigs length distribution", null, "Contigs length", assemblies, [{'name': "Contigs length", 'data': data}] );
}


/**
 * Returns the HTML table representation of mapping metrics.
 * @param step_data {Hash} By sample the mapping metrics.
 * @return {String} The HTML table.
 */
var mapping_table = function( step_data ){
	var categories = ["Sample", "Reads 1", "Reads 2", "Paired", "Mapped", "Properly paired", "Mate on other chr"] ;
	var data = new Array();
	for( var sample in step_data ){
		var mapped_prct = (step_data[sample].mapped/(step_data[sample].nb_read_1 + step_data[sample].nb_read_2))*100 ;
		var properly_paired_prct = 0 ;
		if( step_data[sample].paired != 0 ){
			properly_paired_prct = (step_data[sample].properly_paired/step_data[sample].paired)*100
		}
		data.push(
			new Array(
				sample == "default" ? "-" : sample,
				step_data[sample].nb_read_1,
				step_data[sample].nb_read_2,
				step_data[sample].paired,
				numberDisplay(mapped_prct) + "%",
				numberDisplay(properly_paired_prct) + "%",
				step_data[sample].mate_on_other_chr
			)
		);
	}
	
	return table( "Mapping information", categories, data ) ;
}


/**
 * Returns the hash to init an Highchart's barplot from orientation by contigs.
 * @example
 *     $('#test').highcharts( orientation_chart(orientation.data) );
 * @param step_data {Hash} By sample the number of contigs by % of R1 in forward orientation.
 * @param ordered_samples {Array} The names of samples in display order.
 * @return {Hash} Parameter to use in Highchart's constructor.
 */
var orientation_chart = function( step_data, ordered_samples ){
	var data = new Array() ;
	if( ordered_samples == null ){
		ordered_samples = Object.keys( step_data );
	}
	for( var sample_idx = 0 ; sample_idx < ordered_samples.length ; sample_idx++ ){
		var sample = ordered_samples[sample_idx] ;
		var series = {
				'name': sample,
				'data': new Array(),
				'index': sample_idx
		};
		var categories = Object.keys(step_data[sample].contigs_by_ratio).sort(function(a, b){return a-b});
		var nb_contigs = 0 ;
		for( var idx = 0 ; idx < categories.length ; idx++ ){
			var curr_category = categories[idx] ;
			nb_contigs += step_data[sample].contigs_by_ratio[curr_category]
		};
		for( var idx = 0 ; idx < categories.length ; idx++ ){
			var curr_category = categories[idx] ;
			series['data'].push([ parseInt(curr_category), step_data[sample].contigs_by_ratio[curr_category] ]);
		}
		data.push( series );
	}
	var chart = barplot("% of R1 in forward orientation by contig", "% of R1 in forward orientation", "Nb contigs", null, data, 30) ;
	chart['plotOptions']['column']['dataLabels']['enabled'] = false ;
	return chart ;
}


/**
 * Returns the HTML table representation of gene set completeness metrics.
 * @param step_data {Hash} The BUSCO metrics.
 * @return {String} The HTML table.
 */
var geneSet_table = function( step_data ){
	var categories = ["Sample", "lineage", "Complete", "Complete single-copy", "Complete duplicated", "Fragmented", "Missing", "Total searched", "BUSCO notation"] ;
	var data = new Array();
	for( var sample in step_data ){
		data.push(
			new Array(
					(sample == "default" ? "-" : sample),
					step_data[sample].lineage,
					step_data[sample].complete,
					step_data[sample].complete_single_copy,
					step_data[sample].complete_duplicated,
					step_data[sample].fragmented,
					step_data[sample].missing,
					step_data[sample].total,
					step_data[sample].notation
				)
		);
	}
	
	return table( "Recoved genes classification", categories, data ) ;
};

/**
 * Returns the HTML img tag of gene set completeness metrics.
 * @return {String} The HTML img tag.
 */

var geneSet_img = function( ){
	var src = '../busco/busco_figure.png';
	var width = '50%';
	return '<img src="' + src + '" class="img-responsive center-block" width="' + width + '">';
}

/**
 * Returns the HTML table representation of reference proteins alignment metrics.
 * @param step_data {Hash} The number of aligned proteins and the number of matched contigs.
 * @return {String} The HTML table.
 */
var proteinAln_table = function( step_data ){
	var categories = ["Sample", "Nb contig with protein(s)", "Nb protein(s) aligned on contig"] ;
	var data = new Array();
	for( var sample in step_data ){
		data.push(
			new Array(
					(sample == "default" ? "-" : sample),
					step_data[sample].nb_contig_with_prot,
					step_data[sample].nb_diff_proteins
				)
		);
	}
	
	return table( "Alignment count", categories, data ) ;
};


/**
 * Returns the hash to init an Highchart's boxplot from assembly contigs scores. 
 * @example
 *     var boxplot_param = scoring_chart( step_data, ordered_samples )
 *     $('#test').highcharts( boxplot_param[0], ordered_samples[1] );
 * @param step_data {Hash} By sample the scores distribution in key length_distribution.
 * @param ordered_samples {Array} The names of samples in display order.
 * @return {Hash} Parameter to use in Highchart's constructor.
 */
var scoring_chart = function( step_data, ordered_samples ){
	var data = new Array();
	if( ordered_samples == null ){
		ordered_samples = Object.keys( step_data );
	}
	var assemblies = new Array();
	for( var sample_idx = 0 ; sample_idx < ordered_samples.length ; sample_idx++  ){
		var sample = ordered_samples[sample_idx] ;
		var metrics = step_data[sample] ;
		assemblies.push(sample);
		data.push( new Array(
			metrics.contigs_scores.min,
			metrics.contigs_scores.lower_quartile,
			metrics.contigs_scores.median,
			metrics.contigs_scores.upper_quartile,
			metrics.contigs_scores.max
		));
	}

	return boxplot("Contigs scores distribution", null, "Contigs score", assemblies, [{'name': "Contigs score", 'data': data}]);
}


/**
 * Returns the HTML table representation of global scores.
 * @param step_data {Hash} By sample the global scores.
 * @return {String} The HTML table.
 */
var scoring_table = function( step_data ){
	var categories = ["Sample", "Score", "Optimal score", "Optimal cutoff"] ;
	var data = new Array();
	for( var sample in step_data ){
		data.push(
			new Array(
				sample == "default" ? "-" : sample,
				step_data[sample].score,
				step_data[sample].optimal_score,
				step_data[sample].optimal_cutoff
			)
		);
	}
	
	return table( "Global score", categories, data, 4 ) ;
}


/**
 * Returns the hash to init an Highchart's barplot from ORF count by contigs. 
 * @example
 *     $('#test').highcharts( transdecoder_chart(step.data.transdecoder) );
 * @param step_data {Hash} By sample the number of contigs by number of ORF in the contig.
 * @param ordered_samples {Array} The names of samples in display order.
 * @return {Hash} Parameter to use in Highchart's constructor.
 */
var transdecoder_chart = function( step_data, ordered_samples ){
	var data = new Array() ;
	if( ordered_samples == null ){
		ordered_samples = Object.keys( step_data );
	}
	for( var sample_idx = 0 ; sample_idx < ordered_samples.length ; sample_idx++ ){
		var sample = ordered_samples[sample_idx] ;
		var series = {
				'name': sample,
				'data': new Array(),
				'index': sample_idx
		};
		var categories = Object.keys(step_data[sample].nb_ORF_distribution).sort(function(a, b){return a-b});
		var nb_contigs = 0 ;
		for( var idx = 0 ; idx < categories.length ; idx++ ){
			var curr_category = categories[idx] ;
			nb_contigs += step_data[sample].nb_ORF_distribution[curr_category]
		};
		for( var idx = 0 ; idx < categories.length ; idx++ ){
			var curr_category = categories[idx] ;
			series['data'].push([ curr_category, Math.round((step_data[sample].nb_ORF_distribution[curr_category]/nb_contigs)*10000)/100 ]);
		}
		data.push( series );
	}
	return barplot("Nb ORFs by contigs", "Nb ORF(s)", "% contigs", null, data, 30) ;
}