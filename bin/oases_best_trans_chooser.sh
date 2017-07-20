#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

# To execute inside RawData directory

if ($1 == "") then
        echo "Give sample directory"
        exit
endif

set script = `readlink -f $0`
set scriptpath = `dirname $script`

cd $1
foreach dir (k*)
	cd $dir
	echo "Clean locus in $dir"
	$scriptpath/Oasesv2.0.4BestTransChooser.py -l 0.8
	cd ..
end
