# DRAP : De novo RNA-seq Assembly Pipeline

**This is a fork** of [DRAP v1.91: http://www.sigenae.org/drap](http://www.sigenae.org/drap) adapted to the cluster of the Staion Biologique de Roscoff.

Modifications are quick and dirty and not intensively tested.
At the moment only `runMeta` works correctly with `sge` set in the `cfg/drap.cfg` file . Everything else must be used `local` option set (but can be submitted as a job with a qsub)

License: [GNU GPLv3](http://www.gnu.org/licenses/gpl-3.0.en.html)

## Description

Short read RNASeq de novo assembly is a well established method to study transcription of organisms lacking a reference genome sequence. Available software packages such as Trinity and Oases have proven to be able to build high quality contigs from short reads. But there is still room for improvement on different points such as:

- compactness: they often produce different contigs which are included in one another or overlapping one another,
- chimerism: the contigs contain different kinds on chimera such as duplicated open reading frames,
- substitution, insertion, deletion errors: the consensus sequences build by the assembler contain errors which can be partly corrected using the read alignments.


DRAP includes three modules:

- runDrap chains an Oases or Trinity assembly of reads from a given sample with several compaction and correction steps. It produces several assembly files with different FPKM threshold for total contigs or contigs comprising an open reading frame. A report file presents the resulting assembly and alignment metrics.
- runMeta gathers all the samples assemblies and fusions the results in a unique representative contig set. It also removes the redundancy between sets and produces a general reports including assembly and alignment metrics.
- runAssessment processes different contigs sets build from the same read sets to generate assembly and alignment metrics which are collected in report. It helps to choose the best assembly.

## Docker install

Go to the original install page: [http://www.sigenae.org/drap/install.html](http://www.sigenae.org/drap/install.html)

## Local install

### Dependencies:

- Programming languages interpreters and modules
    - bash
    - csh
    - perl 5.*.* with non standard modules:
        - Bio::Search::Hit::GenericHit
        - Bio::Search::Tiling::MapTiling
        - Bio::SearchIO::Writer::HSPTableWriter
        - Bio::SeqIO
        - Bio::Tools::Run::StandAloneBlast
        - IPC::Run
        - JSON
        - List::Util
        - Term::ANSIColor
    - python 2.7.* with non standard modules:
        - Bio
        - NumPy
        - SciPy

- External softwares
    - [bedtools >= 2.22.1](http://bedtools.readthedocs.org/en/latest/content/installation.html)
    - [blat >= 35](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
    - [bowtie >= 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
    - [busco >= 3.0](http://busco.ezlab.org/)
    - [bwa >= 0.7.15](http://sourceforge.net/projects/bio-bwa/files/)
    - [cd-hit >= 4.6](https://github.com/weizhongli/cdhit)
    - [cutadapt >= 1.8.3](https://github.com/marcelm/cutadapt)
    - [dc >= 1.3.95 (in bc >=1.06)](https://www.gnu.org/software/bc/)
    - [exonerate >= 2.2.0](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
    - [express >= 1.5.1](https://github.com/adarob/eXpress)
    - [fastq_illumina_filter >= 0.1](http://cancan.cshl.edu/labmembers/gordon/fastq_illumina_filter/#download)
    - [getorf >= EMBOSS: 6.4.0.0](http://emboss.sourceforge.net/download)
    - [khmer >= 2.0](https://github.com/dib-lab/khmer)
    - [NCBI Blast+ >= 2.2.29](ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.2.29/)
    - [NCBI Tools >= 12.0.0](ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/CURRENT/)
    - [oases >= 0.2.06](https://github.com/dzerbino/oases/tree/master)
    - [parallel >= parallel-20141022](http://ftp.gnu.org/gnu/parallel)
    - [rsync >= 3.0.6](http://www.htslib.org/download/)
    - [samtools >= 1.3](ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/seqclean)
    - [seqclean](ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/seqclean/)
    - [STAR >= 2.4.0i](https://github.com/alexdobin/STAR/releases)
    - [tgicl < 2](ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/tgicl/) -> yes, it is **not** [tgicl >= 2.1](https://sourceforge.net/projects/tgicl)
    - [TransDecoder >= 2.0.1](http://transdecoder.github.io/)
    - [TransRate >= 1.0.1](http://hibberdlab.com/transrate/installation.html)
    - [trim_galore >= 0.4.0](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
    - [Trinity >= 2.4.0](https://github.com/trinityrnaseq/trinityrnaseq/releases)


Details about how those softwares are used can be see in [doc/third_party_tools.html](./doc/third_party_tools.html)

### Configure & Test

see the [doc/install.html](./doc/install.html) and [doc/quick_start.html](./doc/quick_start.html) pages.

## Hacks

### Run runMeta without having used runDrap before.

You must use a command similar to the following:

```
#!/bin/bash

DRAP_PATH="/usr/local/genome2/drap"
WORKING_DIR="$(pwd)"
OUT_FOLDER="$WORKING_DIR"

$DRAP_PATH/runMeta \
 --cfg-file $WORKING_DIR/cfg/drap.cfg \
 --drap-dirs $OUT_FOLDER/trinity_splA,$OUT_FOLDER/trinity_splB \
 --ref $DRAP_PATH/test/data/Danio_rerio.pep.fasta \
 --outdir $OUT_FOLDER/meta_trinity \
```

where:
`--drap-dirs` are the folders obtained from `runDrap`.
Each of those folders must contains at least le following contents in order to successfully run `runMeta`:

 - `.drap_conf.json` (used in the steps `06-meta_index.sh`, `07-meta_rmbt.sh` and `09-meta_postprocess.sh` of `runMeta`): a json file containing at least the following elements:
 ```
 {
   "alignR1" : [
      "/path/to/sampleA_R1.fastq.gz"
   ],
   "alignR2" : [
      "/path/to/sampleA_R2.fastq.gz"
   ],
   "coverages" : [
      "1",
      "3",
      "5",
      "10"
   ],
   "paired" : 1,
   "strand" : null
  }
```

 - `transcripts_fpkm_X.fa` (used in `01-meta_merge.sh`): which is the file of transcripts to be mapped.
 The `X` in `transcripts_fpkm_X.fa` name must be the minimal value in the list associated to the "coverages" key in the file `.drap_conf.json`. The value of `X` is the *coverage cutoff* used by `express` when transcripts (from `transcripts_fpkm.fa`) are filtered.

**Notes**:

 - The `alignR2` key  can be ommited when the `paired` key is set to 0.
 - If we set, for example, `"coverages" : ["2"]` and `transcripts_fpkm_2.fa`, the file `.drap_conf.json` produced by `runMeta` will put back `"coverages" : ["1", "3", "5", "10"]`.
