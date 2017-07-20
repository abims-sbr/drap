#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

# Usage: cov_length_filter.sh -f FASTA -x FILE [OPTIONS]
# 
# -f FASTA       MANDATORY contig fasta file
# -x FILE        MANDATORY eXpress result file
# -b FILE        ORF length tsv file (MANDATORY if type is orf)
# -o FASTA       Output directory
# -t STRING      Type must be contig or orf [contig]
# -l INTEGER     Length cutoff in nucleotides for the specified type [200]
# -c STRING      Coverage cutoff range [1,2,3,4,5]
# --log FILE     Log file
#
# END HELP

foreach a ($argv)
	if (! $?option) then
		set option = $a
	else
		switch ($option)
			case "-f" :
				set fasta = `readlink -f $a`
				breaksw
			case "-x" :
				set express = `readlink -f $a`
				breaksw
			case "-b" :
				set orfLength = `readlink -f $a`
				breaksw
			case "-o" :
				set outputDir = $a
				breaksw
 			case "-t" :
				set type = $a
				breaksw
			case "-l" :
				set len = $a
				breaksw
			case "-c" :
				set cov = $a
				breaksw
			case "--log" :
				set logfile = $a
				breaksw
		endsw
		unset option
	endif
end	

if !($?fasta) then
	echo "PLEASE GIVE CONTIG FASTA FILE"
	goto HELP
else if ($fasta == "") then
	echo "NO SUCH FASTA FILE $fasta"
	goto HELP
else if !(-f $fasta) then
	echo "NO SUCH FASTA FILE $fasta"
	goto HELP
endif

if !($?outputDir) then
	set outputDir = `pwd`
else if ($outputDir == "") then
	set outputDir = `pwd`
else
	set outputDir = `readlink -f $outputDir`
endif

if !($?logfile) then
	set logfile = $outputDir/`basename $0 .sh`.log
else if ($logfile == "") then
	set logfile = $outputDir/`basename $0 .sh`.log
else
	mkdir -p `dirname $logfile`
	set logfile = `readlink -f $logfile`
endif

if !($?express) then
	echo "PLEASE GIVE EXPRESS RESULT FILE"
	goto HELP
else if !(-f $express) then
	echo "NO SUCH EXPRESS RESULT FILE $express"
	goto HELP
endif

if !($?type) then
	set type = 'contig'
else if !($type == "contig" || $type == "orf") then
	echo "TYPE MUST BE contig OR orf"
	goto HELP
endif

if ($type == "orf") then
	if !($?orfLength) then
		echo "PLEASE GIVE ORF LENGTH TSV FILE"
		goto HELP
	else if !(-f $orfLength) then
		echo "NO SUCH ORF LENGTH TSV FILE $orfLength"
		goto HELP
	endif
endif

if !($?len) then
	set len = 200
endif

if !($?cov) then
	set cov = '1,2,3,4,5'
endif

if !(-d $outputDir) then
	mkdir -p $outputDir
endif

set script = `readlink -f $0`
set scriptpath = `dirname $script`
set cfgpath = `dirname $scriptpath`/cfg

if ($type == "orf") then
	$scriptpath/get_longest_orf.pl -f $fasta -stats | awk -v l=$len '$4>=l' | cut -f 1 | sort > $outputDir/all_contigs_length_threshold.lst
	cat $orfLength | awk -v l=$len '$2>=l' | cut -f 1 | sort > $outputDir/coding_contigs_length_threshold.lst
else
	cat $fasta | $scriptpath/fasta_length.pl | awk -v l=$len '$2>=l' | cut -f 1 | sort > $outputDir/all_contigs_length_threshold.lst
endif
set coverage
set filtered
if ($type == "orf") set coding
foreach c (`echo $cov | tr ',' ' '`)
	set tmp = `mktemp`
	sed 1d $express | awk -v c=$c '$11>=c' | cut -f2 | sort > $tmp
	set coverage = "${coverage},"`cat $tmp | wc -l`
	comm -12 $tmp $outputDir/all_contigs_length_threshold.lst > $outputDir/all_contigs_cov_length_thresholds.lst
	cat $fasta | $scriptpath/fasta_extract.pl $outputDir/all_contigs_cov_length_thresholds.lst > $outputDir/transcripts_fpkm_$c.fa
	set filtered = "${filtered},"`cat $outputDir/transcripts_fpkm_$c.fa | grep -c '^>'`
	if ($type == "orf") then
		comm -12 $tmp $outputDir/coding_contigs_length_threshold.lst > $outputDir/coding_contigs_cov_length_thresholds.lst
		cat $fasta | $scriptpath/fasta_extract.pl $outputDir/coding_contigs_cov_length_thresholds.lst > $outputDir/coding_transcripts_fpkm_$c.fa
		set coding = "${coding},"`cat $outputDir/coding_transcripts_fpkm_$c.fa | grep -c '^>'`
	endif
	rm -f $tmp
end

printf "length:        %s contigs pass filter\n" `cat $outputDir/all_contigs_length_threshold.lst | wc -l` > $logfile
if ($type == "orf") printf "coding length: %s coding contigs pass filter\n" `cat $outputDir/coding_contigs_length_threshold.lst | wc -l` >> $logfile
printf "coverage:      %s contigs pass filter\n" `echo $coverage | sed -e 's/,//'` >> $logfile
printf "input:         %s contigs\n" `cat $fasta | grep -c '^>'` >> $logfile
printf "output:        %s contigs\n" `echo $filtered | sed -e 's/,//'` >> $logfile
if ($type == "orf") printf "               %s coding contigs\n" `echo $coding | sed -e 's/,//'` >> $logfile
rm -f $outputDir/*.lst

exit 0

HELP:
set n=`grep -n '^# END HELP' $0 | cut -d : -f 1`
set n=`echo "$n 1 - pq"|dc`
head -$n $0 | tail -n +3
exit 1
