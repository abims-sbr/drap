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
 * Returns the string representation of the number. 
 * @param pValue {Float} The number to process.
 * @param pDecimal_precision {Int} The number of kept decimals (default: 2).
 * @return {String} The string representation (example: 12856892.11111 => 12,856,892.11).
 */
var numberDisplay = function( pValue, pDecimal_precision ){
	var precision = pDecimal_precision ;
	if( pDecimal_precision == null ){
		precision = 2 ;
	}
	var new_val = "" ;
	if( ("" + pValue + "").indexOf(".") != -1 ){
		new_val = pValue.toFixed(precision).replace(/^(\d)(?=(\d{3})+\b)\.?/g, '$1,');
	} else {
		new_val = pValue.toFixed().replace(/(\d)(?=(\d{3})+\b)/g, '$1,');
	}
    return new_val ;
}