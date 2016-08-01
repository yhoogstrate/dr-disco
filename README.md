Dr. Disco
=========

Running RNA-STAR with fusion settings produces a ``....Chimeric.out.sam`` file, containing all discordant reads. This alignment does not properly link mates together because of some incorrect SAM-tags. This tool, `dr-disco fix-chimeric`, is able to solve these issues and allows you to view the discordant reads in IGV in much more detail. Because the tool makes the **disco**rdant alignments **healthy** again, we've named it **Dr. Disco**.

Installation
------------
```
git clone https://github.com/yhoogstrate/dr-disco.git
cd dr-disco
pip install -r requirements.txt ; # Use pip2 in case you have a py3 system
python2 setup.py install
```

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
Usage: dr-disco intronic [OPTIONS] OUTPUT_FILE FUSION_CANDIDATES_INPUT_FILE
                         BAM_INPUT_FILE

Options:
  --help  Show this message and exit.
```
