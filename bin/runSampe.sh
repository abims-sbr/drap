#!/bin/csh

# Usage: runSampe -r FASTA -1 FASTQ -2 FASTQ [options]
#
# -r FASTA       REFERENCE fasta file or directory
# -1 FASTQ       FASTQ read file - R1
# -2 FASTQ       FASTQ mate file - R2
# -o DIRECTORY   OUTPUT directory [reference directory]
# -m STRING      Map with bwa or star [bwa]
# -g             REFERENCE is STAR genomeDir [unset]
#                Need to specify OUTPUT directory
# -s             SAMTOOLS sort, index, flagstat, idxstats with default parameters [unset]
# -t INTEGER     Nb threads for bwa aln or STAR [4]
# -y             Filter Casava-filtered sequences (only with bwa) [unset]
# -p 'PARAMS'    STAR mapping params BETWEEN SINGLE QUOTES
#                Default --alignIntronMin 10 --alignIntronMax 25000 --outFilterMultimapNmax 10000
# -q STRING      Run jobs with queue other than default queue [workq]
# --flagstat     Generate samtools flagstat file [unset]
# --sort_by_name Sort bam by read name (no flagstat no idxstats) [unset]
# --sort_by_pos  Sort bam by reference position (no flagstat no idxstats) [unset]
# --index        Index the REFERENCE fasta file and exit [unset]
# --force        Index the REFERENCE fasta file even if it has been already done [unset]
# --bam          Return only bam result file name to stdout [unset]
# --sync         Wait for last submitted job to end [unset]
# --local        Run commands sequentially (no qsub) [unset]
#
# Output files written into reference directory
# Index the REFERENCE fasta file if it has not been done
# CAREFUL:
# - do not launch more than one runSampe with same non indexed REFERENCE fasta file to prevent concurrent indexing processes
# - names of paired fastq files MUST CONTAIN following strings:
#	 _R1_/_R2_ or _R1./_R2. or .R1./.R2. or _1_/_2_ or _1./_2. or .1./.2.
# - or names of paired fastq files MUST START with following strings:
#	 R1_/R2_ or R1./R2. or 1_/2_ or 1./2.
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
			case "-1" :
				set R1 = $a
				breaksw
			case "-2" :
				set R2 = $a
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
echo "runSampe $argv:q" >> $output/runSampe.$job_id.cmd

if !($?R1) then
	set msg = 'PLEASE GIVE R1 FASTQ FILE'
	goto HELP
else if !(-f $R1) then
	set msg = "NO SUCH FILE $R1"
	goto HELP
endif

if !($?R2) then
	set msg = 'PLEASE GIVE R2 FASTQ FILE'
	goto HELP
else if !(-f $R2) then
	set msg = "NO SUCH FILE $R2"
	goto HELP
endif

set baseR1 = `basename $R1`
set baseR2 = `basename $R2`
if (`echo $baseR1 | awk '$0~/[_.]R?1\.|_R?1_|^R?1[._]/'` == "") then
	set msg = "NAME OF PAIRED FASTQ files MUST CONTAIN following strings:#_R1_/_R2_ or _R1./_R2. or .R1./.R2. or _1_/_2_ or _1./_2. or .1./.2.#OR MUST START with following strings:#R1_/R2_ or R1./R2. or 1_/2_ or 1./2."
	goto HELP
endif

set basePair = `perl -le '@R1=split(//,shift);@R2=split(//,shift);map{unless($R1[$_]eq$R2[$_]){$d++;$i=$_}}0..$#R1;if($#R1==$#R2&&$R1[$i]==1&&$R2[$i]==2&&$d==1){$R1[$i]=12;print join("",@R1)}' $baseR1 $baseR2`
if ($basePair == '') then
	set msg = "NAMES OF PAIRED FASTQ files MUST ONLY DIFFER BY 1 or 2 AT THE SAME POSITION"
	goto HELP
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

if !($?ref_is_genome_dir) set ref_is_genome_dir = 'false'
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
if !($?queue)    set queue    = 'workq';
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
		echo $index_cmd >> $output/runSampe.$job_id.log
		( $index_cmd > $output/err_log_$job_id/index.o$job_id ) >& $output/err_log_$job_id/index.e$job_id
	else
		if ($index == 'true' && $sync == 'true') then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		set index_cmd = "$scriptpath/submitJob --name index $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSampe.$job_id.log $index_opt--binary -- $index_cmd"
		echo $index_cmd >> $output/runSampe.$job_id.log
		set index_id = `$index_cmd`
		if ($bam == 'false') tail -2 $output/runSampe.$job_id.log
		echo $index_id > $output/err_log_$job_id/runSampe.$job_id.jid
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
set aln_opt
if !("$hold" == '' && "$pe" == '') set aln_opt = "--options $hold $pe "
set ext = `echo $reference:e`
set ref = `basename $reference .$ext`

set R1out = `echo $baseR1 | sed -e 's/\.gz$//;s/\.fastq$//;s/\.fq$//'`
set R2out = `echo $baseR2 | sed -e 's/\.gz$//;s/\.fastq$//;s/\.fq$//'`
set paiRout = `echo $basePair | sed -e 's/\.gz$//;s/\.fastq$//;s/\.fq$//'`
if ($bam == 'true') echo "$output/${paiRout}_to_$ref.$job_id.bam"

if ($mapper == 'bwa') then
	set illumina = ''
	if ($filter == 'true') set illumina = '-Y'
	echo "$bwa_cmd aln -t $threads $illumina -f $output/${R1out}_to_$ref.$job_id.sai $reference $R1" > $output/err_log_$job_id/aln.$job_id.sh
	echo "$bwa_cmd aln -t $threads $illumina -f $output/${R2out}_to_$ref.$job_id.sai $reference $R2" >> $output/err_log_$job_id/aln.$job_id.sh
	set aln_cmd = "$output/err_log_$job_id/aln.$job_id.sh"
	if ($local == 'true') then
		if ($bam == 'false') echo $aln_cmd
		echo $aln_cmd >> $output/runSampe.$job_id.log
		set task = 1
		foreach line ("`cat $aln_cmd`")
			( csh -c "$line" > $output/err_log_$job_id/aln.o$job_id.$task ) >& $output/err_log_$job_id/aln.e$job_id.$task
			@ task++
		end
	else
		set aln_cmd = "$scriptpath/submitJob --array --name aln --queue $queue --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSampe.$job_id.log $aln_opt-- $aln_cmd"
		echo $aln_cmd >> $output/runSampe.$job_id.log
		set aln_id = `$aln_cmd`
		if ($bam == 'false') tail -2 $output/runSampe.$job_id.log
		echo $aln_id >> $output/err_log_$job_id/runSampe.$job_id.jid
		set hold = "-hold_jid $aln_id"
		set pe = ''
		set map_res = "$bwa_map_res"
		set threads = 1
	endif
	set sampe_cmd = "$bwa_cmd sampe $reference $output/${R1out}_to_$ref.$job_id.sai $output/${R2out}_to_$ref.$job_id.sai $R1 $R2"
else
	set read_cmd = ''
	if (`echo $R1:e` == 'gz') set read_cmd = "--readFilesCommand zcat"
	set sampe_cmd = "$star_cmd --genomeDir $star_genome_dir --readFilesIn $R1 $R2 --runThreadN $threads $star_params --outFileNamePrefix $output/${paiRout}_to_$ref.$job_id. $read_cmd --outStd SAM"
	set map_res = "$star_map_res"
endif
set syncp = ''
if ($stats == 'false' && $name == 'false' && $pos == 'false') then
	if ($sync == 'true') set syncp = '--sync'
	set sort = "> $output/${paiRout}_to_$ref.$job_id.bam"
else
	set by_name = ''
	if ($name == 'true') set by_name = '-n '
	set sort = "| samtools sort -@ $threads -m 4G -O bam -T $output/tmp.$job_id ${by_name}-o $output/${paiRout}_to_$ref.$job_id.bam -"
endif
echo "$sampe_cmd | samtools view -bS -@ $threads - $sort" > $output/err_log_$job_id/sampe.$job_id.sh
set sampe_cmd = "$output/err_log_$job_id/sampe.$job_id.sh"
if ($local == 'true') then
	if ($bam == 'false') echo $sampe_cmd
	echo $sampe_cmd >> $output/runSampe.$job_id.log
	( csh $sampe_cmd > $output/err_log_$job_id/sampe.o$job_id ) >& $output/err_log_$job_id/sampe.e$job_id
else
	set sampe_cmd = "$scriptpath/submitJob --name sampe --queue $queue $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSampe.$job_id.log --options $pe $hold $map_res -- $sampe_cmd"
	echo $sampe_cmd >> $output/runSampe.$job_id.log
	set sampe_id = `$sampe_cmd`
	if ($bam == 'false') tail -2 $output/runSampe.$job_id.log
	echo $sampe_id >> $output/err_log_$job_id/runSampe.$job_id.jid
endif
if ($stats == 'true') then
	echo "samtools index $output/${paiRout}_to_$ref.$job_id.bam\nsamtools flagstat $output/${paiRout}_to_$ref.$job_id.bam > $output/${paiRout}_to_$ref.$job_id.bam.flagstat\nsamtools idxstats $output/${paiRout}_to_$ref.$job_id.bam > $output/${paiRout}_to_$ref.$job_id.bam.idxstats" > $output/err_log_$job_id/stats.$job_id.sh
	set stats_cmd = "$output/err_log_$job_id/stats.$job_id.sh"
	if ($local == 'true') then
		if ($bam == 'false') echo $stats_cmd
		echo $stats_cmd >> $output/runSampe.$job_id.log
		( csh $stats_cmd > $output/err_log_$job_id/stats.o$job_id ) >& $output/err_log_$job_id/stats.e$job_id
	else
		if ($sync == 'true') then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		set stats_cmd = "$scriptpath/submitJob --name stats $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSampe.$job_id.log --options -hold_jid $sampe_id $samtools_idx_res -- $stats_cmd"
		echo $stats_cmd >> $output/runSampe.$job_id.log
		set stats_id = `$stats_cmd`
		if ($bam == 'false') tail -2 $output/runSampe.$job_id.log
		echo $stats_id >> $output/err_log_$job_id/runSampe.$job_id.jid
	endif
endif
if ($flagstat == 'true') then
	echo "samtools flagstat $output/${paiRout}_to_$ref.$job_id.bam > $output/${paiRout}_to_$ref.$job_id.bam.flagstat" > $output/err_log_$job_id/flagstat.$job_id.sh
	set flagstat_cmd = "$output/err_log_$job_id/flagstat.$job_id.sh"
	if ($local == 'true') then
		if ($bam == 'false') echo $flagstat_cmd
		echo $flagstat_cmd >> $output/runSampe.$job_id.log
		( csh $flagstat_cmd > $output/err_log_$job_id/flagstat.o$job_id ) >& $output/err_log_$job_id/flagstat.e$job_id
	else
		if ($sync == 'true') then
			set syncp = '--sync'
		else
			set syncp = ''
		endif
		set flagstat_cmd = "$scriptpath/submitJob --name flagstat $syncp --stdout $output/err_log_$job_id --stderr $output/err_log_$job_id --log $output/runSampe.$job_id.log --options -hold_jid $sampe_id -- $flagstat_cmd"
		echo $flagstat_cmd >> $output/runSampe.$job_id.log
		set flagstat_id = `$flagstat_cmd`
		if ($bam == 'false') tail -2 $output/runSampe.$job_id.log
		echo $flagstat_id >> $output/err_log_$job_id/runSampe.$job_id.jid
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
echo "ERROR: $msg" | sed -e 's/#/\n/g' > $output/runSampe.$job_id.log
exit 1
