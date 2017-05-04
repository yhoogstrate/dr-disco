#!/bin/bash

source .venv/bin/activate
python setup.py install

samtools view -bS tests/detect-intronic/test_26.sam > /tmp/test.bam
dr-disco fix /tmp/test.fixed.bam /tmp/test.bam
dr-disco detect test.out /tmp/test.fixed.bam
cat test.out
