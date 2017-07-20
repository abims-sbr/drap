#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

set fastq = $1
set dir = $2

set basename = `basename $fastq`
if (`readlink -f $fastq | xargs file | grep -c compressed` == 0) then
	ln -s $1 $dir/$basename.unzip
else
	zcat $1 > $dir/$basename.unzip
endif
