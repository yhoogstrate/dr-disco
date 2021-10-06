# Dr. Disco

### Detecting exon-to-exon and genomic breakpoints in RNA-Seq data

[![Build Status](https://travis-ci.org/yhoogstrate/dr-disco.svg?branch=master)](https://travis-ci.org/yhoogstrate/dr-disco)
 
 - Free software: GNU General Public License v3 (GPLv3)

## Introduction

Fusion genes are often driver events in several types of cancer, generally detected with DNA-sequencing (DNA-seq) and more recently with RNA-sequencing (RNA-seq).
In classical RNA-seq experiments, mRNA is extracted by targeting their poly-A tails, or reverse transcribed using oligo-dT primers. 
Multiple studies have shown that fusion transcripts can be successfully detected in such data, but the location of the DNA break remains often undetected because the majority of the breaks are in introns and poly-A RNA-seq data lacks intronic sequencing reads.
By reverse transcribing rRNA-depleted RNA using random hexamer priming, the RNA library will also include pre-mRNA, covering introns in subsequent alignments.

Here we propose a computational method for detecting intron-to-intron breaks, corresponding to DNA breaks, on top of exon-to-exon splice junctions of fusion genes, using RNA-seq data only.
To detect genomic breakpoints, intronic and exonic data are kept separated in a graph and further analysed.
The analyses of the graph has built-in solutions for sequencing fragments rather than entire transcripts, for sequencing only the fragments' ends and for multiple exon-to-exon boundaries on top of genomic breaks.
Unlike most RNA-seq fusion gene detection software, it is not restricted to exons or gene annotations, allowing detection of novel splice junctions, fusions to non-gene regions and fusions of non-polyadenylated transcripts.

## Installation

### Local contained installation (recommended)

However, using virtual environments gives more control over dependencies and local installations may be more convenient. This can be achieved as follows:

```
git clone https://github.com/yhoogstrate/dr-disco.git
cd dr-disco
virtualenv -p python3 .venv
source .venv/bin/activate
python setup.py install
```

### Local installation

Dr. Disco makes use of python3. A typical system wide installation could be achieved as follows:

```
git clone https://github.com/yhoogstrate/dr-disco.git
cd dr-disco
pip install -r requirements.txt ;
python3 setup.py install --user
```

### System installation

```
git clone https://github.com/yhoogstrate/dr-disco.git
cd dr-disco
sudo pip install -r requirements.txt ; 
sudo python3 setup.py install --user
```

### Pre-compiled

Dr. Disco is also available at BioConda but this does not automatically ship with the blacklist files. Also, because of the changes in the built-system of bioconda, this version may be out of sync with the git repo.

[![Bio-Conda installer](https://cdn.rawgit.com/yhoogstrate/dr-disco/master/share/bioconda-badge.svg)](share/bioconda-badge.svg)

### Usage

The method is designed to analyse a single sample at the time, but can easily be run using linux parallel to scale up.

A `dr-disco` pipeline typically consists of the following steps:

1. *Run STAR* in fusion mode to obtain the Chimeric sam file - By the time of writing we used STAR 2.4.2
2. *Fix* the Chimeric sam file and prepare it for analysis
3. *Detect* fusion events and merge junctions corresponding to similar fusion transcripts
4. *Classify* fusion genes based on statistical parameters (and reference genome specific blacklist)
5. *Integrate* integrates results from the same genomic event and annotates gene names

#### Check mate pair orentation first
STAR does not seem to be able to appropriately determine Chimeric reads if the orientation of the mates is not default. Hence, before running STAR, make sure whether the strand orientation is `FR` (_R1: forward, _R2: _reverse). If this is not the case, you can use **fastx_reverse_complement** (http://hannonlab.cshl.edu/fastx_toolkit/) to convert to the appropriate strands.

 - For *forward*, *reverse* (default) you can proceed with the original _R1 and _R2
 - For *forward*, *forward* you need to use the original _R1 and create and use the reverse complement of _R2.
 - For *reverse*, *forward* you need to create and use the reverse complement of _R1 and o_R2.
 - For *reverse*, *reverse* you need to create and use the reverse complement of _R1 and the original of _R2 (but I am not aware of such data...).

#### Illumina NextSeq data cleaning

The illumina NextSeq is claimed to be responsible for measuring poly-G sequences (https://www.biostars.org/p/294612/) which may result in vast amounts of discordant reads aligned to poly-G regions in the genome.
It is recommended to clean such datasets and erase poly-G suffixes, for instance with the tool fastp.

### Step 1: STAR

Running RNA-STAR with fusion settings produces a file: ``<...>.Chimeric.out.sam``. This file contains discordant reads (split and spanning). It is recommended to run STAR with corresponding settings:

```
STAR --genomeDir ${star_index_dir} \  
     --readFilesIn ${left_fq_filenames} ${right_fq_filenames} \  
     --outFileNamePrefix {$prefix} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMstrandField intronMotif \
     --outFilterIntronMotifs None \
     --alignIntronMax 200000 \
     --alignMatesGapMax 200000 \
     --alignSJDBoverhangMin 10 \
     --alignEndsType Local \
     --chimSegmentMin 12 \
     --chimJunctionOverhangMin 12 \
     --sjdbGTFfile {$gene_model_gtf} \
     --sjdbOverhang 100 \
     --quantMode GeneCounts \
     --twopass1readsN -1 \
     --twopassMode Basic
```

Due to the `chimSegmentMin` and `chimJunctionOverhangMin` settings, STAR shall produce the additional output file(s). This file may need to be converted to BAM format, e.g. by running `samtools view -bhS Chimeric.out.sam > Chimeric.out.bam`. This generated `ChimericSorted.bam` can then be used as input file for *Dr. Disco*. Note that samtools has changed its commandline interface over time and that the stated command might be slighly different for different versions of samtools.

### Step 2: dr-disco fix

The first step of Dr. Disco is fixing the BAM file. Altohugh fixing implies it is broken, the alignment provided by STAR is incompatible with the IGV split view and misses the 'SA:' sam flag. In `dr-disco fix` this is solved and allows a user to view discordant reads in IGV in more detail and by using the spit-view. For the split view it may be convenient to color the reads by subtype (split or spanning, and how the direction of the break is), as the mate and strand are related and discordant mates are usually shifted a bit with respect to the breakpoint. In the fixing steps such annotations are added as read groups. Also, this step ensures proper indexing and sorting. If you have as input file `<...>.Chimeric.out.bam`, you can generate the fixed bam file with *Dr. Disco* as follows:

```
Usage: dr-disco fix [OPTIONS] INPUT_ALIGNMENT_FILE OUTPUT_ALIGNMENT_FILE

Options:
  -t, --temp-dir PATH  Path in which temp files are stored (default: /tmp)
  --help               Show this message and exit.
```

Here the INPUT_ALIGNMENT_FILE will be `<...>.Chimeric.out.bam` and OUTPUT_ALIGNMENT_FILE will be `<...>.Chimeric.fixed.bam`.
Hence, you can generate the fixed bam-file with:

```
dr-disco fix '<...>.Chimeric.out.bam' 'sample.fixed.bam'
```

### Step 3: dr-disco detect

To estimate breakpoints using *Dr. Disco*, you can proceed with:

```
Usage: dr-disco detect [OPTIONS] BAM_INPUT_FILE OUTPUT_FILE

Options:
  -m, --min-e-score INTEGER  Minimal score to initiate pulling sub-graphs
                             (larger numbers boost performance but result in
                             suboptimal results) [default=8]
```

Here the `-m` argument controls the level of merging sub-graphs. If datasets become very large there shall be many subgraphs. However, because this is such a time consuming process, it can be desired to skip some of them. Rule of thumb: for every 10.000.000 reads, increase this value with 1. So for 10.000.000-20.000.000 mate pairs choose 1. So for a dataset of 75.000.000 mate pairs, proceed with:

`dr-disco detect -m 7 'sample.fixed.bam' dr-disco.sample-name.all.txt`


### Step 4: dr-disco classify

The results from `dr-disco detect` contain many false positives. These may be derived from many different types of issues, like homolpolymeric sequences in the reference genome, rRNA, multimaps and optimal chimeric alignments just by chance. Most of these issues are recognisable by patterns imprinted in some variables in the output table. Then `dr-disco classify` filters the results and you can proceed with:

```
Usage: dr-disco classify [OPTIONS] TABLE_INPUT_FILE TABLE_OUTPUT_FILE

Options:
  --only-valid                 Only return results marked as 'valid'
  --blacklist-regions TEXT     Blacklist these regions (BED file)
  --blacklist-junctions TEXT   Blacklist these region-to-region junctions
                               (custom format, see files in ./share/)
  --min-chim-overhang INTEGER  Minimum alignment length on each side of the
                               junction. May need to be set to smaller values
                               for read lengths smaller than 75bp. Larger
                               values are more stringent. [default=50]
  --ffpe                       Lowers the threshold for the relative amount of
                               mismatches, as often found in FFPE material.
                               Note that enabling this option will
                               consequently result in more false positives.
  --help                       Show this message and exit.
```

After analysing our results after classification there were some positives that were systemetically recurrent, but also visual inspection suggested them to be true. These are often new exons or new genes (and marked as discordant probably due to non canonical splice junctions) or mistakes in the reference genome. Hence, as they are in terms of data exactly identical to fusion genes, these can not be filtered out. Therefore we added the option to blacklist certain regions or junctions. There is an important difference between them as a blacklist region prohibits any fusion to be within a certain region while a blacklist junction prohibits a fusion where the one break is in one region and the other break in the other. We have added quite some of these regions to the blacklist. In particular the blacklist for hg38 is rather complete. The blacklist files for hg19 and hh38 can be found here: [./share/](https://github.com/yhoogstrate/dr-disco/tree/master/share).

It may be inconvenient to include the results that are classified as false. By using `--only-valid` the results are restricted to only those classified positive. Hence, a typical way of using `dr-disco classify` for a hg38 reference genome would be:

```
dr-disco classify \
    --only-valid \
    --blacklist-regions "$dr_disco_dir/share/blacklist-regions.hg38.bed" \
    --blacklist-junctions "$dr_disco_dir/blacklist-junctions.hg38.txt" \
    dr-disco.sample-name.all.txt \
    dr-disco.sample-name.filtered.txt
```

### Step 5: dr-disco integrate

The results of Dr. Disco classify can be integrated into fusion genes by running `dr-disco integrate`. The corresponding GTF file should contain gene annotations that correspond to the used reference genome. The reference FASTA is used to determine splice junction motifs and their edit distance. The file should have a *.fai index file. I am not sure if the pyfaidx library will auto-generate one if it does not exist.

```
Usage: dr-disco integrate [OPTIONS] TABLE_INPUT_FILE TABLE_OUTPUT_FILE

Options:
  --gtf TEXT    Use gene annotation for estimating fusion genes and improve
                classification of exonic (GTF file)
  --fasta TEXT  Use reference sequences to estimate edit distances to splice
                junction motifs (FASTA file)
  --help        Show this message and exit.
```

Hence, we could proceed with:

```
dr-disco integrate \
    --gtf ENSEMBL 'GRCh38_latest_genomic.gff' \
    dr-disco.sample-name.filtered.txt \
    dr-disco.sample-name.integrated.txt
```
