#!/bin/csh

# Usage: runSamse -r FASTA -f FASTQ [options]
#
# -r FASTA       REFERENCE fasta file or directory
# -f FASTQ       FASTQ reads file(s) - separated by coma if more than one
# -o DIRECTORY   OUTPUT directory [reference directory]
# -m STRING      Map with bwa or star [bwa]
# -g             REFERENCE is STAR genomeDir [unset]
#                Need to specify OUTPUT directory
# -s             SAMTOOLS sort, index, flagstat, idxstats with default parameters [unset]
# -t INTEGER     Nb threads for bwa aln or STAR [4]
# -n INTEGER     Nb reads to extract from fastq(s) [unset]
# -y             Filter Casava-filtered sequences [unset]
# -p 'PARAMS'    STAR mapping params BETWEEN SINGLE QUOTES
#                Default --alignIntronMin 10 --alignIntronMax 25000 --outFilterMultimapNmax 10000
# -q STRING      Run align jobs with queue other than default queue [workq]
# --flagstat     Generate samtools flagstat file
# --sort_by_name Sort bam by read name (no flagstat no idxstats) [unset]
# --sort_by_pos  Sort bam by reference position (no flagstat no idxstats) [unset]
# --index        Index the REFERENCE fasta file and exit
# --force        Index the REFERENCE fasta file even if it has been already done
# --bam          Return only bam result file name to stdout [unset]
# --sync         Wait for last submitted job to end [false]
# --local        Run commands sequentially (no qsub) [false]
#
# Output files written into reference directory
# Index the REFERENCE fasta file if it has not been done
# CAREFUL: do not launch more than one runSamse with same non indexed REFERENCE fasta file to prevent concurrent indexing processes
#
# COPYRIGHT 2015 INRA
# AUTHORS: Cedric Cabau - INRA Toulouse - sigenae-support@listes.inra.fr
# LICENSE: GNU GPLv3

# END HELP

set bwa_cmd = 'bwa'
set star_cmd = 'STAR'
set script = `readlink -f $0`
set scriptpath = `dirname $script`

# default resources if env variables are unset
if (! $?pe_smp) set pe_smp = 'thread'
if (! $?star_idx_res) set star_idx_res = "-R y -pe $pe_smp 4 -l mem_free=16G,h_vmem=64G"
if (! $?star_map_res) set star_map_res = "-l mem_free=4G,h_vmem=16G"
if (! $?star_genome_generate_ram) set star_genome_generate_ram = 64
if (! $?star_genome_generate_cpu) set star_genome_generate_cpu = 4
if (! $?star_genome_sa_index) set star_genome_sa_index = 14
if (! $?bwa_map_res) set bwa_map_res = "-l mem_free=8G,h_vmem=16G"
if (! $?samtools_idx_res) set samtools_idx_res = "-R y -l mem_free=16G,h_vmem=32G"

foreach a ($argv:q)
	if (! $?option) then
		if ("$a" == "--sync") then
			set sync = 'true'
		else if ("$a" == "--local") then
			set local = 'true'
		else if ("$a" == "--force") then
			set force = 'true'
		else if ("$a" == "--index") then
			set index = 'true'
		else if ("$a" == "--flagstat") then
			set flagstat = 'true'
		else if ("$a" == "--bam") then
			set bam = 'true'
		else if ("$a" == "--sort_by_name") then
			set name = 'true'
		else if ("$a" == "--sort_by_pos") then
			set pos = 'true'
		else if ("$a" == "-g") then
			set ref_is_genome_dir = 'true'
		else if ("$a" == "-s") then
			set stats = 'true'
		else if ("$a" == "-y") then
			set filter = 'true'
		else
			set option = $a
		endif
	else
		switch ($option)
			case "-r" :
				set reference = $a
				breaksw
			case "-f" :
				set fastq = $a
				breaksw
			case "-o" :
				set output = $a
				breaksw
			case "-m" :
				set mapper = $a
				breaksw
			case "-t" :
				set threads = $a
				breaksw
			case "-n" :
				set nb = $a
				breaksw
			case "-p" :
				set star_params = "$a"
				breaksw
			case "-q" :
				set queue = $a
				breaksw
			case "--bwa_cmd" :
				set bwa_cmd = $a
				breaksw
		endsw
		unset option
	endif
end

if !($?bam) set bam = 'false'
if !($?output) then
	if ($?reference) then
		set output = `dirname $reference`
	else
		set msg = 'PLEASE GIVE REFERENCE FASTA FILE OR OUPUT DIRECTORY'
		goto HELP
	endif
else
	mkdir -p $output
endif
set job_id = $$
echo "runSamse $argv:q" >> $output/runSamse.$job_id.cmd

if !($?fastq) then
	set msg = 'PLEASE GIVE READ FASTQ FILE'
	goto HELP
else
	foreach file ( `echo $fastq | tr ',' '\n'` )
		if !(-e $file) then
			set msg = "NO SUCH FASTQ FILE $file"
			goto HELP
		endif
	end
endif

if !($?reference) then
	set msg = 'PLEASE GIVE REFERENCE FASTA FILE'
	goto HELP
else if !(-e $reference) then
	set msg = "NO SUCH FILE or DIRECTORY $reference"
	goto HELP
endif

if !($?mapper) then
	set mapper = 'bwa'
else
	set mapper = `echo $mapper | awk '{print tolower($1)}'`
endif
if !($mapper == 'bwa' || $mapper == 'star') then
	set msg = "MAPPER MUST BE bwa OR star"
	goto HELP
endif

if !($?ref_is_genome_dir) then
	set ref_is_genome_dir = 'false'
endif
if !($ref_is_genome_dir == 'false') then
	if !(-d $reference) then
		set msg = "-g OPTION only AVAILABLE if REFERENCE is a DIRECTORY"
		goto HELP
	else if !($mapper == 'star') then
		set msg = "-g OPTION only AVAILABLE with star"
		goto HELP
	else if !($?output) then
		set msg = "-o OPTION is MANDATORY with -g OPTION"
		goto HELP
	endif
endif

if !($?stats)    set stats    = 'false'
if !($?name)     set name     = 'false'
if !($?pos)      set pos      = 'false'
if !($?threads)  set threads  = 4
if !($?filter)   set filter   = 'false'
if ($?default_queue)    set queue    = $default_queue
if !($?queue)    set queue    = 'workq'
if !($?sync)     set sync     = 'false'
if !($?local)    set local    = 'false'
if !($?force)    set force    = 'false'
if !($?index)    set index    = 'false'
if !($?flagstat) set flagstat = 'false'
if !($?star_params) set star_params = '--alignIntronMin 10 --alignIntronMax 25000 --outFilterMultimapNmax 10000 --outSAMunmapped Within'
if ($filter != 'false' && $mapper != 'bwa') then
	set msg = "-y OPTION only AVAILABLE with bwa"
	goto HELP
endif
if ($flagstat == 'true' && $stats == 'true') set flagstat = 'false'
if (($name == 'true' || $pos  == 'true') && $stats == 'true') set stats = 'false'
if ($mapper == 'star') then
	set ref_basename = `basename $reference`
	set star_genome_dir = "$output/STAR_{$ref_basename}_$job_id"
	if ($ref_is_genome_dir == 'true') set star_genome_dir = $reference
endif

mkdir -p $output/err_log_$job_id

if ($?nb) then
	set input_dir = $output/runSamse.$job_id.fastq
	mkdir -p $input_dir
	foreach file ( `echo $fastq | tr ',' '\n'` )
		if !(`echo $file | sed 's/\.gz$//;s/.*\.//'` == 'fastq' || `echo $file | sed 's/\.gz$//;s/.*\.//'` == 'fq') then
			echo 'OPTION -n ONLY AVAILABLE WITH FASTQ FILES'
			goto HELP
		endif
		if (`echo $file:e` == 'gz') then
			set cmd = 'zcat'
		else
			set cmd = 'cat'
		endif
		set head = `expr $nb \* 4`
		set fname = `basename $file .gz`
		set illumina = ''
		if ($filter == 'true') set illumina = '| fastq_illumina_filter -N'
		$cmd $file $illumina | head -$head > $input_dir/$fname
		if ($?new_file) then
			set new_file = `echo "${new_file},$input_dir/$fname"`
		else
			set new_file = $input_dir/$fname
		endif
	end
	set fastq = $new_file
endif

set need_index = 1;
if ($mapper == 'bwa' && -e "$reference.amb") set need_index = 0
if ($mapper == 'star' && $ref_is_genome_dir == 'true') set need_index = 0
if ($force == 'true') set need_index = 1
if ($need_index == 1) then
	set index_opt = ''
	if ($mapper == 'bwa') then
		set index_cmd = "$bwa_cmd index $reference"
	else
		mkdir -p $star_genome_dir
		set index_cmd = "$star_cmd --runMode genomeGenerate --genomeSAindexNbases $star_genome_sa_index --limitGenomeGenerateRAM ${star_genome_generate_ram}000000000 --genomeDir $star_genome_dir --genomeFastaFiles $reference --outFileNamePrefix $output/ --runThreadN $star_genome_generate_cpu"
		set index_opt = "--options $star_idx_res "
	endif
	if ($local == 'true') then
		if ($bam == 'false') echo $index_cmd
		echo $index_cmd >> $output/runSamse.$job_id.log
		( $index_cmd > $output/err_log_$job_id/index.o$job_id ) >& $output/err_log_$job_id/index.e$job_id
	else
		if ($index == 'true' && $sync == 'true') then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		set index_cmd = "$scriptpath/submitJob --name index $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSamse.$job_id.log $index_opt--binary -- $index_cmd"
		echo $index_cmd >> $output/runSamse.$job_id.log
		set index_id = `$index_cmd`
		if ($bam == 'false') tail -2 $output/runSamse.$job_id.log
		echo $index_id > $output/err_log_$job_id/runSamse.$job_id.jid
	endif
endif

if ($index == 'true') exit 0

if ($?index_id) then
	set hold = "-hold_jid $index_id"
else
	set hold = ''
endif
if ($threads > 1) then
	set pe = "-R y -pe $pe_smp $threads"
else
	set pe = ''
endif
set ext = `echo $reference:e`
set ref = `basename $reference .$ext`
set n = `echo $fastq | grep -c ','`

foreach file ( `echo $fastq | tr ',' '\n'` )
	set outname = `echo $file | sed -e 's/\.gz$//;s/\.fastq$//;s/\.fq$//' | xargs basename`
	if ($bam == 'true') echo "$output/${outname}_to_$ref.$job_id.bam"
	if ($mapper == 'bwa') then
		set illumina = ''
		if ($filter == 'true') set illumina = '-Y'
		set samse_cmd = "$bwa_cmd aln -t $threads $illumina $reference $file | $bwa_cmd samse $reference - $file"
		set map_res = "$bwa_map_res"
		set threads = 1
	else
		set read_cmd = ''
		if (`echo $file:e` == 'gz') set read_cmd = "--readFilesCommand zcat"
		set samse_cmd = "$star_cmd --genomeDir $star_genome_dir --readFilesIn $file --runThreadN $threads $star_params --outFileNamePrefix $output/${outname}_to_$ref.$job_id. $read_cmd --outStd SAM"
		set map_res = "$star_map_res"
	endif
	if ($stats == 'true' || $pos == 'true') then
		echo "$samse_cmd | samtools view -bS -@ $threads - | samtools sort -@ $threads -m 4G -O bam -T $output/tmp.$job_id -o $output/${outname}_to_$ref.$job_id.bam -" >> $output/err_log_$job_id/samse.$job_id.sh
		if ($stats == 'true') then
			echo "samtools index $output/${outname}_to_$ref.$job_id.bam;samtools flagstat $output/${outname}_to_$ref.$job_id.bam > $output/${outname}_to_$ref.$job_id.bam.flagstat;samtools idxstats $output/${outname}_to_$ref.$job_id.bam > $output/${outname}_to_$ref.$job_id.bam.idxstats" >> $output/err_log_$job_id/stats.$job_id.sh
		endif
	else if ($name == 'true') then
		echo "$samse_cmd | samtools view -bS -@ $threads - | samtools sort -@ $threads -m 4G -n -O bam -T $output/tmp.$job_id -o $output/${outname}_to_$ref.$job_id.bam -" >> $output/err_log_$job_id/samse.$job_id.sh
	else
		echo "$samse_cmd | samtools view -bS -@ $threads - > $output/${outname}_to_$ref.$job_id.bam" >> $output/err_log_$job_id/samse.$job_id.sh
	endif
	if ($flagstat == 'true') then
		echo "samtools flagstat $output/${outname}_to_$ref.$job_id.bam > $output/${outname}_to_$ref.$job_id.bam.flagstat" >> $output/err_log_$job_id/flagstat.$job_id.sh
	endif
end

if (($stats == 'false') && !($?nb) && ($sync == 'true')) then
	set syncp = '--sync'
else
	set syncp = ''
endif
set samse_cmd = "$output/err_log_$job_id/samse.$job_id.sh"
if ($local == 'true') then
	if ($bam == 'false') echo $samse_cmd
	echo $samse_cmd >> $output/runSamse.$job_id.log
	set task = 1
	foreach line ("`cat $samse_cmd`")
		( csh -c "$line" > $output/err_log_$job_id/samse.o$job_id.$task ) >& $output/err_log_$job_id/samse.e$job_id.$task
		@ task++
	end
else
	set samse_cmd = "$scriptpath/submitJob --array --name samse $syncp --queue $queue --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSamse.$job_id.log --options $pe $hold $map_res -- $samse_cmd"
	echo $samse_cmd >> $output/runSamse.$job_id.log
	set samse_id = `$samse_cmd`
	if ($bam == 'false') tail -2 $output/runSamse.$job_id.log
	echo $samse_id >> $output/err_log_$job_id/runSamse.$job_id.jid
endif
if ($stats == 'true') then
	set stats_cmd = "$output/err_log_$job_id/stats.$job_id.sh"
	if ($local == 'true') then
		if ($bam == 'false') echo $stats_cmd
		echo $stats_cmd >> $output/runSamse.$job_id.log
		set task = 1
		foreach line ("`cat $stats_cmd`")
			( csh -c "$line" > $output/err_log_$job_id/stats.o$job_id.$task ) >& $output/err_log_$job_id/stats.e$job_id.$task
			@ task++
		end
	else
		if (!($?nb) && ($sync == 'true')) then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		if ($n > 0) then
			set hold = "-hold_jid_ad $samse_id"
		else
		set hold = "-hold_jid $samse_id"
		endif
		set stats_cmd = "$scriptpath/submitJob --array --name stats $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSamse.$job_id.log --options $hold $samtools_idx_res -- $stats_cmd"
		echo $stats_cmd >> $output/runSamse.$job_id.log
		set stats_id = `$stats_cmd`
		if ($bam == 'false') tail -2 $output/runSamse.$job_id.log
		echo $stats_id >> $output/err_log_$job_id/runSamse.$job_id.jid
	endif
else if ($flagstat == 'true') then
	set flagstat_cmd = "$output/err_log_$job_id/flagstat.$job_id.sh"
	if ($local == 'true') then
		if ($bam == 'false') echo $flagstat_cmd
		echo $flagstat_cmd >> $output/runSamse.$job_id.log
		set task = 1
		foreach line ("`cat $flagstat_cmd`")
			( csh -c "$line" > $output/err_log_$job_id/flagstat.o$job_id.$task ) >& $output/err_log_$job_id/flagstat.e$job_id.$task
			@ task++
		end
	else
		if (!($?nb) && ($sync == 'true')) then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		if ($n > 0) then
			set hold = "-hold_jid_ad $samse_id"
		else
			set hold = "-hold_jid $samse_id"
		endif
		set flagstat_cmd = "$scriptpath/submitJob --array --name flagstat $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSamse.$job_id.log --options $hold -- $flagstat_cmd"
		echo $flagstat_cmd >> $output/runSamse.$job_id.log
		set flagstat_id = `$flagstat_cmd`
		if ($bam == 'false') tail -2 $output/runSamse.$job_id.log
		echo $flagstat_id >> $output/err_log_$job_id/runSamse.$job_id.jid
	endif
endif

endif

if ($?nb) then
	set clean_cmd = "rm -rf $input_dir"
	if ($local == 'true') then
		if ($bam == 'false') echo $clean_cmd
		echo $clean_cmd >> $output/runSamse.$job_id.log
		( $clean_cmd > $output/err_log_$job_id/clean.o$job_id ) >& $output/err_log_$job_id/clean.e$job_id
	else
		if ($sync == 'true') then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		set clean_cmd = "$scriptpath/submitJob --name clean $syncp --stdout $output/err_log_$job_id/rmSamse.$job_id.log --options -j y -hold_jid $samse_id --binary -- $clean_cmd"
		echo $clean_cmd >> $output/runSamse.$job_id.log
		set clean_id = `$clean_cmd`
		if ($bam == 'false') tail -2 $output/runSamse.$job_id.log
	endif
endif

exit 0

HELP:
set n=`grep -n '^# END HELP' $0 | cut -d : -f 1`
set n=`echo "$n 1 - pq"|dc`
if ($bam == 'false') then
	echo "ERROR: $msg" | sed -e 's/#/\n/g'
	head -$n $0 | tail -n +3
endif
echo "ERROR: $msg" | sed -e 's/#/\n/g' > $output/runSamse.$job_id.log
exit 1
