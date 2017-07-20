#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

set fastq = $1
set clean = $2
set script = `readlink -f $0`
set scriptpath = `dirname $script`

if (`readlink -f $fastq | xargs file | grep -c compressed` == 0) then
	set cmd = 'cat'
else
	set cmd = 'zcat'
endif

if ($clean == 1 && `$scriptpath/test_illumina_filter.sh $fastq` != 1) then
	$cmd $fastq | fastq_illumina_filter -vN -o /dev/null | tail -1 | cut -d' ' -f2 | tr -d ','
else
	set n = `$cmd $fastq | wc -l`
	echo "$n 4 / pq" | dc
endif
