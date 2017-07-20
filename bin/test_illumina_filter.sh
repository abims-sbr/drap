#!/bin/bash
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

fastq=$1

if [ `readlink -f $fastq | xargs file | grep -c compressed` == 0 ]; then
        cmd='cat'
else
        cmd='zcat'
fi

$cmd $fastq | head -4 | fastq_illumina_filter -N >& /dev/null
echo $?
