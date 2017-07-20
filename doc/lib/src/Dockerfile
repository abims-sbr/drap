FROM ubuntu
MAINTAINER sigenae "sigenae-support@listes.inra.fr"
ENV REFRESHED_AT 2017-04-20

## -- install packages -- ##
RUN apt-get update && apt-get install -y \
	bedtools \
	berkeley-express \
	bioperl \
	bowtie2 \
	bwa \
	bzip2 \
	cd-hit \
	csh \
	curl \
	dc \
	default-jre \
	emboss \
	exonerate \
	g++ \
	gcc \
	git-hub \
	khmer \
	libc6-i386 \
	libcompress-raw-zlib-perl \
	libgetopt-long-descriptive-perl \
	libipc-run-perl  \
	libjson-perl \
	libncurses5-dev \
	libpod2-base-perl \
	libterm-extendedcolor-perl \
	cpanminus \
	make \
	nano \
	ncbi-blast+ \
	ncbi-tools-bin \
	parallel \
	python \
	python-biopython \
	python-cutadapt \
	python-scipy \
	r-base \
	r-cran-ggplot2 \
	rna-star \
	rssh \
	ruby-dev \
	subversion \
	transdecoder \
	wget \
	zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

## -- Last version of List::Util Perl module -- ##
RUN cpanm List::Util

## -- create links -- ##
RUN cd /usr/bin \
	&& ln -s $PWD/usr/bin/env /bin/env \
	&& ln -s $PWD/python /bin/python \
	&& ln -s $PWD/cdhit /usr/local/bin/cd-hit \
	&& ln -s $PWD/cdhit-est /usr/local/bin/cd-hit-est \
	&& ln -s $PWD/berkeley-express /usr/local/bin/express

## -- set ENV -- ##
ENV PATH /usr/lib/khmer/bin:$PATH
ENV USER root

## -- install velvet/oases with long kmers -- ##
RUN cd /usr/local/src \
	&& git clone --recursive https://github.com/dzerbino/velvet.git \
	&& cd velvet \
	&& make 'MAXKMERLENGTH=70' 'LONGSEQUENCES=1' \
	&& ln -s $PWD/velveth /usr/local/bin/velveth \
	&& ln -s $PWD/velvetg /usr/local/bin/velvetg 
RUN cd /usr/local/src \
	&& git clone --recursive https://github.com/dzerbino/oases  \
	&& cd oases && make 'MAXKMERLENGTH=70' 'LONGSEQUENCES=1' \
 	&& ln -s $PWD/oases /usr/local/bin/oases

## -- install trimGalore -- ##
RUN cd /usr/local/src \
	&& git clone --recursive https://github.com/FelixKrueger/TrimGalore.git \
	&& cd TrimGalore \
	&& ln -s $PWD/trim_galore /usr/local/bin/trim_galore

## -- install blat -- ##
RUN cd /usr/local/bin \
	&& wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
	&& chmod +x blat

## -- install fastq_illumina_filter -- ##
RUN cd /usr/local/bin \
	&& wget http://cancan.cshl.edu/labmembers/gordon/fastq_illumina_filter/release/0.1/fastq_illumina_filter-Linux-x86_64 \
	&& chmod +x fastq_illumina_filter-Linux-x86_64 \
	&& ln -s fastq_illumina_filter-Linux-x86_64 fastq_illumina_filter

## -- install seqclean -- ##
## unable to dowload archives from occams.dfci.harvard.edu with curl; use wget ##
RUN cd /usr/local/src \
	&& wget ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/seqclean/seqclean_x86_64.tar.gz \
	&& tar -xzf seqclean_x86_64.tar.gz && rm -f seqclean_x86_64.tar.gz \
	&& ln -s $PWD/seqclean/seqclean /usr/local/bin/seqclean
ENV PATH /usr/local/src/seqclean/bin:$PATH

## -- install tgicl -- ##
RUN cd /usr/local/src \
	&& wget ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/tgicl/tgicl_linux.tar.gz  \
	&& tar -xzf tgicl_linux.tar.gz && rm -f tgicl_linux.tar.gz \
	&& wget ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/tgicl/gclib.tar.gz \
	&& tar -xzf gclib.tar.gz && rm -f gclib.tar.gz \
	&& wget ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/tgicl/zmsort.tar.gz \
	&& tar -xzf zmsort.tar.gz && rm -f zmsort.tar.gz \
	&& curl -SL http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar | tar -x \
	&& cd zmsort && make && cp zmsort ../tgicl_linux/bin && cd .. \
	&& cp CAP3/cap3 tgicl_linux/bin \
	&& ln -s $PWD/tgicl_linux/tgicl /usr/local/bin/tgicl
ENV PATH /usr/local/src/tgicl_linux/bin:$PATH
ENV PERL5LIB /usr/local/src/tgicl_linux:$PERL5LIB

## -- install trinity -- ##
RUN cd /usr/local/src \
	&& curl -SL https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.tar.gz | tar -xz \
	&& cd trinityrnaseq-Trinity-v2.4.0 \
	&& make && make plugins \
	&& ln -s $PWD/Trinity /usr/local/bin/Trinity \
	&& ln -s $PWD/util/insilico_read_normalization.pl /usr/local/bin/insilico_read_normalization.pl \
	&& rm -rf $PWD/sample_data
ENV PATH /usr/local/src/trinityrnaseq-Trinity-v2.4.0/trinity-plugins/BIN:$PATH

## -- install transRate -- ##
RUN cd /usr/local/src \
	&& curl -SL https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz | tar -xz \
	&& cd transrate-1.0.3-linux-x86_64 \
	&& ./transrate --install-deps ref \
	&& rm -f bin/librt.so.1
ENV PATH /usr/local/src/transrate-1.0.3-linux-x86_64:/usr/local/src/transrate-1.0.3-linux-x86_64/bin:/usr/local/src/transrate-1.0.3-linux-x86_64/lib/app/ruby/2.2.0/bin:$PATH

## -- install busco -- ##
RUN cd /usr/local/src \
	&& git clone --recursive https://gitlab.com/ezlab/busco.git \
	&& cd busco \
	&& curl -SL http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz | tar -xz \
	&& python setup.py install

RUN cd /usr/local/src/busco \
	&& mv config/config.ini.default config/config.ini \
	&& sed -i -e 's|/home/osboxes.*|/usr/bin/|' config/config.ini \
	&& chmod +x scripts/run_BUSCO.py scripts/generate_plot.py \
	&& ln -s $PWD/scripts/run_BUSCO.py /usr/local/bin/run_BUSCO.py \
	&& ln -s $PWD/scripts/generate_plot.py /usr/local/bin/generate_plot.py

## -- loading the adapter and contamination databases for TSA compliance -- ##
RUN mkdir /usr/local/src/banks \
	&& cd /usr/local/src/banks \
	&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa \
	&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz && gzip -d contam_in_euks.fa.gz \
	&& wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/mito.nt.gz && gzip -d mito.nt.gz \
	&& makeblastdb -in adaptors_for_screening_euks.fa -dbtype nucl \
	&& makeblastdb -in contam_in_euks.fa -dbtype nucl \
	&& makeblastdb -in mito.nt -dbtype nucl

## -- install DRAP -- ##
RUN cd /usr/local/src \
	&& svn checkout svn://scm.mulcyber.toulouse.inra.fr/svnroot/drap/tags/drap-v1.91 \
	&& cd drap-v1.91 \
	&& ln -s $PWD/runCheck /usr/local/bin/runCheck \
	&& ln -s $PWD/runDrap /usr/local/bin/runDrap \
	&& ln -s $PWD/runMeta /usr/local/bin/runMeta \
	&& ln -s $PWD/runAssessment /usr/local/bin/runAssessment
ENV PERL5LIB /usr/local/src/drap-v1.91/bin:$PERL5LIB
ENV PATH /usr/local/src/drap-v1.91/bin:$PATH

## -- configure DRAP -- ##
RUN cd /usr/local/src/drap-v1.91 \
	&& sed -i -e 's|/path/to/adaptors_for_screening_euks.fa|/usr/local/src/banks/adaptors_for_screening_euks.fa|' cfg/drap.cfg \
	&& sed -i -e 's|/path/to/contam_in_euks.fa|/usr/local/src/banks/contam_in_euks.fa|' cfg/drap.cfg \
	&& sed -i -e 's|/path/to/mito.nt|/usr/local/src/banks/mito.nt|' cfg/drap.cfg \
	&& sed -i -e 's|/path/to/busco_directory|/usr/local/src/busco|' cfg/drap.cfg \
	&& sed -i -e 's|type = sge|type = local|' cfg/drap.cfg \
	# correcting tgicl conflict version
	&& sed -i -e 's/tgicl -F/tgicl/' bin/runAssembly.sh
