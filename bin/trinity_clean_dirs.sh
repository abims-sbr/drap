#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

if ("$1" == "") then
	echo 'trinity_clean_dirs.sh [-no] <directory>'
	exit
endif

if ("$1" == "-no") then
	set no_Seq = 1
	set new_argv = `echo $argv | cut -d' ' -f2-`
	set argv = "$new_argv"
endif

foreach dir ($argv)
	echo "Clean $dir"
	find $dir -mindepth 1 -maxdepth 1 -not -name Trinity.\* -exec rm -rf {} \;
	if !($?no_Seq) then
		echo "Clean Sequences file"
		find $dir/.. -maxdepth 1 -type f -name Sequences -exec truncate -s 0 {} \;
	endif
end
