Dr. Disco
=========

Running RNA-STAR with fusion settings will produce a ``....Chimeric.out.sam`` file, containing all discordant reads. This alignment does not properly link mates together by having some incorrectly used tags. This tool, fix-chimeric-out-bam.py, is able to solve these issues allowing you to view the discordant reads in IGV in much more detail. Because the tool makes the *disco*rdant alignment health again, we've named it **Dr. Disco**.

Installation
------------
None yet, please make sure pysam 0.9, and click are avaiable to python, e.g. using `pip install ...`. Also samtools >= 1.3 is required.

Usage: dr-disco fix
-------------------
If you have as input file `....Chimeric.out.sam`, you should run Dr. Disco as follows:

```
samtools view -bS '....Chimeric.out.sam' > '....Chimeric.out.bam' ;
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
Usage: dr-disco intronic [OPTIONS] OUTPUT_FILE FUSION_CANDIDATES_INPUT_FILE
                         BAM_INPUT_FILE

Options:
  --help  Show this message and exit.
```
