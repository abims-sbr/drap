#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

# tgicl_clean_dirs.sh $dir/c-asm

foreach dir ($argv)
	echo " directory $dir cleaned"
	find $dir -maxdepth 1 -type d \( -name err_log -or -name asm_\* \) -exec rm -rf {} \;
	find $dir -type f -not \( -name tgicl.cfg -or -name \*.log -or -name runAssembly.sh -or -name all_dbg.fa -or -name all_contigs.raw.fa -or -name -or -name all_dbg.fa.clstr \) -exec rm -rf {} \;
end
