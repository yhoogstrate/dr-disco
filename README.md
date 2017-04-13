Dr. Disco
=========

Detection exon-to-exon and genomic breakpoints in RNA-Seq data
--------------------------------------------------------------

[![Build Status](https://travis-ci.org/yhoogstrate/dr-disco.svg?branch=master)](https://travis-ci.org/yhoogstrate/dr-disco)
 
 - Free software: GNU General Public License v3 (GPLv3)

Introduction
------------

Fusion genes are often driver events in several types of cancer, generally detected with DNA-sequencing (DNA-seq) and more recently with RNA-sequencing (RNA-seq).
In classical RNA-seq experiments, mRNA is extracted by targeting their poly-A tails, or reverse transcribed using oligo-dT primers. 
Multiple studies have shown that fusion transcripts can be successfully detected in such data, but the location of the DNA break remains often undetected because the majority of the breaks are in introns and data lacks intronic sequencing reads.
By reverse transcribing rRNA-depleted RNA using random hexamer priming, the RNA library will also include pre-mRNA, covering introns in subsequent alignments.

Here we introduce a computational method for detecting intron-to-intron breaks, corresponding to DNA breaks, on top of exon-to-exon splice junctions of fusion genes, using RNA-seq data only.
To detect genomic breakpoints, intronic and exonic data are kept separated in a graph and further analysed, with software package  *Dr. Disco*.
The analyses of the graph has built-in solutions for sequencing fragments rather than entire transcripts, for sequencing only the fragments' ends and for multiple exon-to-exon boundaries on top of genomic breaks.
Unlike most RNA-seq fusion gene detection software, it is not restricted to exons or gene annotations, allowing detection of novel splice junctions, fusions to non-gene regions and fusions of non-polyadenylated transcripts.

In this study its relevance is demonstrated by generating a map of TMPRSS2-ERG breakpoints in a cohort of 50 prostate cancer samples, by using only RNA-seq data.
TMPRSS2-ERG was predicted in 30 of the 50 patients, of which 28 were found to have intron-to-intron boundaries.
A map shows in both genes a region enriched with breakpoints, suggesting preference for genomic recombinations in these sub regions.
In one intron in ERG, most breakpoints are found only in the first part, a transcription factor binding site rich region.
On top of the advantages RNA-seq has, like detecting splicing patterns and expression profiles, random primed RNA-seq data allows to detect genomic breakpoints.
Apart from detecting genomic breaks of fusion genes, the software also predicted multiple deletions within a single intron of TMPRSS2.
Because such deletions do not affect the mRNA, they will not be detected in poly-A RNA-seq data.
To support downstream analysis and allow integration with workflow management systems, the software makes use of community standard file formats and is therefore also compatible with genome browsers to allow data visualization.
The software is available in Bio-Conda and the Galaxy platform under a free and open software license.

Installation
------------
```
git clone https://github.com/yhoogstrate/dr-disco.git
cd dr-disco
pip install -r requirements.txt ; # Use pip2 in case you have a py3 system
python2 setup.py install --user
```

[![Bio-Conda installer](https://cdn.rawgit.com/yhoogstrate/dr-disco/master/share/bioconda-badge.svg)](share/bioconda-badge.svg)

Usage
=====

Running RNA-STAR with fusion settings produces a ``....Chimeric.out.sam`` file, containing all discordant reads. This alignment does not properly link mates together because of some incorrect SAM-tags. This tool, `dr-disco fix-chimeric`, is able to solve these issues and allows you to view the discordant reads in IGV in much more detail. Because the tool makes the **disco**rdant alignments **healthy** again, we've named it **Dr. Disco**.

Usage: STAR
-----------

It is recommended to run STAR with corresponding settings:

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
     --twopass1readsN 18446744073709551615 \
     --twopassMode Basic
```

Due to the `chimSegmentMin` and `chimJunctionOverhangMin` settings, STAR shall produce additional output files including `Chimeric.out.sam`. This file needs to be converted to BAM format, e.g. by running `samtools view -bhS Chimeric.out.sam > Chimeric.out.bam`. This generated `ChimericSorted.bam` is the input file for *Dr. Disco*. Note that samtools has changed its commandline interface over time and that the stated command might be slighly different for different versions of samtools. Regardless, the only thin important is to generate a BAM file of the discordant reads sam file (sorting and indices will be taken care of by dr-disco itself).

Usage: dr-disco fix
-------------------
The first step of Dr. Disco is fixing the BAM file. Fixing? Is it broken then? It is not in particular broken, but in order to view the file with IGV and get reads properly linked, certain links have to put in place. Also, it is convenient to color the reads by subtype, as the mate and strand are related and discordant mates are usually shifted a bit with respect to the breakpoint. In the fixing steps such annotations are added as read groups. Also, this step ensures proper indexing and sorting. If you have as input file `Chimeric.out.bam`, you can generate the fixed bam file with *Dr. Disco* as follows:

```
Usage: dr-disco fix [OPTIONS] OUTPUT_BAM_FILE INPUT_BAM/SAM_FILE

Options:
  -t, --temp-dir PATH  Pathin which temporary files will be stored (default:
                       /tmp)
  --help               Show this message and exit.
```

Hence, `dr-disco fix Chimeric.out.fixed.bam Chimeric.out.bam` will do the conversion for you.

Usage: dr-disco detect
----------------------
To estimate break points using *Dr. Disco*, you can proceed with:

```
Usage: dr-disco detect [OPTIONS] OUTPUT_FILE BAM_INPUT_FILE

Options:
  -m, --min-e-score INTEGER  Minimal score to initiate pulling sub-graphs
                             (larger numbers boost performance but result in
                             suboptimal results) [default=8]
  --help                     Show this message and exit.
```

Here the `-m` argument controls merging of sub-graphs. If datasets become very large there shall be many subgraphs. However, because this is such a time consuming process, it can be desired to skip some of these.
Rule of thumb: for every 10.000.000 reads, increase this value with 1. So for 10.000.000-20.000.000 mate pairs choose 1. So for a dataset of 75.000.000 mate pairs, proceed with:

`dr-disco detect -m 7 dr-disco.sample-name.out.txt Chimeric.out.fixed.bam`


Currently the results contain many false positives. We're working on understanding the mechanisms behind the noise and we are trying to implement additional background filters.

