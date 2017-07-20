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
 * Returns the HTML table representation of the data. 
 * @param pTitle {String} The title of the table.
 * @param pCategories {Array} The title of each column.
 * @param pData {Array} 2D matrix with row and column data.
 * @param pDecimal_precision {Int} The number of kept decimals (default: null).
 * @return {String} The HTML table representation.
 */
var table = function( pTitle, pCategories, pData, pDecimal_precision  ) {
	// Header
	var table_header = '    <tr>\n<th colspan="' + pCategories.length + '">' + pTitle + '</th>    </tr>\n' ;
	var table_header_line = "" ;
	for(var idx = 0 ; idx < pCategories.length ; idx++){
		table_header_line += "      <th>" + pCategories[idx] + "</th>\n" ;
	}
	table_header += "    <tr>\n" + table_header_line + "    </tr>\n" ;
	table_header = "  <thead>\n" + table_header + "  </thead>\n" ;
	
	// Body
	var table_body = '' ;
	for(var data_idx = 0 ; data_idx < pData.length ; data_idx++){
		var table_body_row = "" ;
		for(var category_idx = 0 ; category_idx < pCategories.length ; category_idx++){
			if( $.isNumeric(pData[data_idx][category_idx]) ) {
				table_body_row += "      <td>" + numberDisplay(pData[data_idx][category_idx], pDecimal_precision) + "</td>\n" ;
			} else {
				table_body_row += "      <td>" + pData[data_idx][category_idx] + "</td>\n" ;
			}
		}
		table_body += "    <tr>\n" + table_body_row + "    </tr>\n" ;
	}
	table_body = "  <tbody>\n" + table_body + "  </tbody>\n" ;

	return '<table class="table table-striped">\n' + table_header + table_body + "</table>\n" ;
}

/**
 * Returns hash use to init HightChart object (without 'type'). 
 * @param pTitle {String} The title of the chart.
 * @param pXTitle {String} The xAxis title.
 * @param pYTitle {String} The yAxis title.
 * @param pXCategories {Array} The title of each category (x scale labels).
 * @param pData {Array} The HightChart series.
 * @return {Hash} The hash.
 * @warning This method use HightChart xAxis.categories.
 */
var chartplot = function( pTitle, pXTitle, pYTitle, pXCategories, pData ) {
	var chart = {
	        title: {
	            text: pTitle
	        },
	        xAxis: {},
	        yAxis: {
	            title: {
	                text: pYTitle
	            }
	        },
	        series: pData,
	        credits: {
	        	enabled: false
	        }
	}
 	if( pXCategories != null ){
 		chart['xAxis']['categories'] = pXCategories ;
 	}
	if( pXTitle != null ){
		chart['xAxis']['title'] = { text: pXTitle } ;
	}
	if( pData.length <= 1 ) {
		chart['legend'] = {'enabled': false};
	} else {
		chart['legend'] = {'enabled': true};
	}
	return chart ;
} 

/**
 * Returns hash use to init HightChart line. 
 * @param pTitle {String} The title of the chart.
 * @param pXTitle {String} The xAxis title.
 * @param pYTitle {String} The yAxis title.
 * @param pXCategories {Array} x scale labels.
 * @param pData {Array} The HightChart series.
 * @return {Hash} Parameters to use in Highchart's constructor.
 */
var lineplot = function(pTitle, pXTitle, pYTitle, pXCategories, pData) {
	var chart = chartplot( pTitle, null, pYTitle, pXCategories, pData );
	chart['chart'] = { type: 'line' };
	return chart ;
}

/**
 * Returns hash use to init HightChart boxplot. 
 * @param pTitle {String} The title of the chart.
 * @param pXTitle {String} The xAxis title.
 * @param pYTitle {String} The yAxis title.
 * @param pXCategories {Array} x scale labels.
 * @param pData {Array} The HightChart series.
 * @return {Array} The first element is an hash parameters and the second element is a zoom function. Use in the same order in Highchart's constructor.
 */
var boxplot = function(pTitle, pXTitle, pYTitle, pXCategories, pData) {
	var chart = chartplot( pTitle, null, pYTitle, pXCategories, pData );
	chart['chart'] = { type: 'boxplot', zoomType: 'y' };
	chart['tooltip'] = { followPointer: true };
	chart['plotOptions'] = { boxplot: {
        fillColor: '#F0F0E0',
        lineWidth: 1,
        medianColor: '#0C5DA5',
        medianWidth: 3,
        stemColor: '#A63400',
        stemDashStyle: 'dot',
        stemWidth: 1,
        whiskerColor: '#3D9200',
        whiskerLength: '20%',
        whiskerWidth: 4
    } };
	chart['yAxis']['min'] = 0 ;
	// Find the smallest min and lower quartile and the tallest max and upper quartile
	var min = null ;
	var max = null ; 
	var min_lower_quartile = null ;
	var max_upperquartile = null ;
	for( var series_idx=0 ; series_idx < pData[0]['data'].length ; series_idx++ ){
		if( min != null ){
			min = Math.min( min, pData[0]['data'][series_idx][0] );
			min_lower_quartile = Math.min( min_lower_quartile, pData[0]['data'][series_idx][1] );
			max_lower_quartile = Math.max( max_lower_quartile, pData[0]['data'][series_idx][3] );
			max = Math.max( max, pData[0]['data'][series_idx][4] );
		} else {
			min = pData[0]['data'][series_idx][0] ;
			min_lower_quartile = pData[0]['data'][series_idx][1] ;
			max_lower_quartile = pData[0]['data'][series_idx][3] ;
			max = pData[0]['data'][series_idx][4] ;
		}
	}
	// Initial zoom level
	var zoom_fct = function(chart){} ;
	if( max > 1000 && ((max - min) > 3*(max_lower_quartile - min_lower_quartile)) ){
		var zoom_ymax = max_lower_quartile + (max_lower_quartile - min_lower_quartile) ;
		zoom_fct = function(chart){
        	this.yAxis[0].setExtremes( chart['yAxis']['min'], zoom_ymax );
        	this.showResetZoom() ;
        }
	}
	
	return( [chart, zoom_fct] );
};

/**
 * Returns hash use to init HightChart column. 
 * @param pTitle {String} The title of the chart.
 * @param pXTitle {String} The xAxis title.
 * @param pYTitle {String} The yAxis title.
 * @param pXCategories {Array} x scale labels.
 * @param pData {Array} The HightChart series.
 * @return {Hash} Parameters to use in Highchart's constructor.
 */
var barplot = function( pTitle, pXTitle, pYTitle, pXCategories, pData ) {
	var chart = chartplot( pTitle, pXTitle, pYTitle, pXCategories, pData );
	chart['chart'] = { type: 'column' };
	chart['plotOptions'] = {
		column: {
			dataLabels: {
				enabled: true
			}
		}
	}
	return chart ;
}

/**
 * Returns hash use to init HightChart pie. 
 * @param pTitle {String} The title of the chart.
 * @param pData {Array} The HightChart series.
 * @param pUnity {String} The unity of the count to display in hover.
 * @return {Hash} Parameters to use in Highchart's constructor.
 */
var pieplot = function( pTitle, pData, pUnity ) {
	var series = [{
        type: 'pie',
        name: pUnity,
        data: pData
	}]
	var plot = chartplot( pTitle, null, null, null, series );
	plot['tooltip'] = {
        pointFormat: '{series.name}: <b>{point.y:,.0f}</b>'
    };
	plot['plotOptions'] = {
        pie: {
            allowPointSelect: true,
            cursor: 'pointer',
            dataLabels: {
                enabled: true,
                format: '<b>{point.name}</b>: {point.percentage:.1f}%',
                style: {
                    color: (Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black'
                }
            }
        }
    };
	return plot ;
};