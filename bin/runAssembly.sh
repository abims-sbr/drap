#!/bin/csh
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

find . -name all_tgicl_singlets.lst -exec rm -f {} \;
tgicl all_dbg.fa -c TGICL_CPU -l 60 -p 96 -s 100000
if (! -d asm_1) then
	set clustering_start = `cat tgicl_all_dbg.fa.log | grep -c '^>>> --- clustering \[all_dbg\.fa\] started at'`
	set clustering_stop = `cat tgicl_all_dbg.fa.log | grep -c '^<<< --- clustering \[all_dbg\.fa\] finished at'`
	set nb_empty_Z_file = `cat err_tgicl_all_dbg.fa.log | grep -c 'gzip: cluster_[0-9]*/\*\.Z: No such file or directory'`
	set global_Z_file_size = `gzip -l all_dbg.fa_cl_tabhits_001.Z | sed -e 's/ \+/\t/g;1d' | cut -f3`
	if ($clustering_start == 1 && $clustering_stop == 1 && -z all_dbg.fa_cl_clusters && $global_Z_file_size == 0 && $nb_empty_Z_file == `echo "TGICL_CPU 1 -pq"|dc`) then
		echo
		echo "######################"
		echo "No contigs to assemble, skip step"
		echo "######################"
		echo
		ln -s all_dbg.fa all_contigs.raw.fa
		grep '^>' all_dbg.fa | tr -d '>' | awk '{print $0"\t"$0}' > all_contigs.raw.fa.tgicl.history.log
		exit
	endif
else
	if (! -e all_tgicl_singlets.lst) then
		cdbyank -l all_dbg.fa.cidx | sort > seqnames_all
		cat asm_*/ACE | grep '^AF '| cut -d ' ' -f 2 | sort -u > seqnames_in_asms
		comm -23 seqnames_all seqnames_in_asms > all_tgicl_singlets.lst
	endif
	if (-e masked.lst) then
		sort -o masked.lst masked.lst
		sort -o all_tgicl_singlets.lst_unclean all_tgicl_singlets.lst
		comm -23 all_tgicl_singlets.lst_unclean masked.lst > all_tgicl_singlets.lst
	endif
	cat asm_*/contigs > all_tgicl_contigs.fa
	BINPATH/fasta_extract.pl < all_dbg.fa all_tgicl_singlets.lst > all_tgicl_singlets.fa
	cat all_tgicl_contigs.fa all_tgicl_singlets.fa > all_contigs.raw.fa
	cat */ACE | perl -lne '$c=$1 if /^CO (\w+) .+/;print "$1\t$c" if /^AF (\S+) .+/' > all_contigs.raw.fa.tgicl.history.log
	cat all_tgicl_singlets.lst | awk '{print $0"\t"$0}' >> all_contigs.raw.fa.tgicl.history.log
endif
