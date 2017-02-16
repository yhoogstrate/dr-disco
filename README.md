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



Usage: dr-disco fix
-------------------
If you have as input file `....Chimeric.out.sam`, you should run Dr. Disco as follows:

```
samtools view -bS '....Chimeric.out.sam' > '....Chimeric.out.bam'
```
```
Usage: dr-disco fix [OPTIONS] OUTPUT_BAM_FILE INPUT_BAM/SAM_FILE

Options:
  -t, --temp-dir PATH  Pathin which temporary files will be stored (default:
                       /tmp)
  --help               Show this message and exit.
```

Usage: dr-disco intronic
------------------------
To estimate break points at the intronic level, you can proceed with:

```
export CIRCOS_DIR=~/circos-vxx-xx/
```

Followed by:
```
Usage: dr-disco intronic [OPTIONS] OUTPUT_FILE FUSION_CANDIDATES_INPUT_FILE
                         BAM_INPUT_FILE

Options:
  --help  Show this message and exit.
```
