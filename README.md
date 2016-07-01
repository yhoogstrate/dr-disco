Dr. Disco
=========

Running RNA-STAR with fusion settings will produce a ``....Chimeric.out.sam`` file, containing all discordant reads. This alignment does not properly link mates together by having some incorrectly used tags. This tool, fix-chimeric-out-bam.py, is able to solve these issues allowing you to view the discordant reads in IGV in much more detail. Because the tool makes the *disco*rdant alignment health again, we've named it **Dr. Disco**.

Installation
------------
None yet, please make sure pysam 0.9, and click are avaiable to python, e.g. using `pip install ...`. Also samtools >= 1.3 is required.

Usage
-----
If you have as input file `....Chimeric.out.sam`, you should run Dr. Disco as follows:

```
samtools view -bS '....Chimeric.out.sam' > '....Chimeric.out.bam' ;
./fix-chimeric-out-bam.py '....Chimeric.out.bam' '....Chimeric.fixed.bam'
```

If you have bam files in advance, you can proceed with:

```
./fix-chimeric-out-bam.py '....Chimeric.out.bam' '....Chimeric.fixed.bam'
```
