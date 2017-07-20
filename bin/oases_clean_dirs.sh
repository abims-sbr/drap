#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

if ("$1" == "") then
	echo 'oases_clean_dirs.sh [-no] <directory>'
	exit
endif

if ("$1" == "-no") then
	set no_Seq = 1
	set new_argv = `echo $argv | cut -d' ' -f2-`
	set argv = "$new_argv"
endif

foreach dir ($argv)
	echo "Clean $dir"
	echo "Clean graph files"
	find $dir -type f \( -name Graph2 -or -name LastGraph -or -name PreGraph -or -name Roadmaps \) -exec truncate -s 0 {} \;
	if !($?no_Seq) then
		echo "Clean Sequences file"
		find $dir/.. -maxdepth 1 -type f -name Sequences -exec truncate -s 0 {} \;
	endif
end
