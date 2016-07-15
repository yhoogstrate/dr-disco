#!/bin/bash

##for D in `ls samples/*.fixed.bam`
##do
##    echo "$D";
##    samtools view "$D" "chr21:39,722,880-40,070,424" | grep "chr21,42" > /tmp/a.txt
##    samtools view "$D" "chr21:42,834,478-42,882,085" | grep -P "chr21,39|chr21,40" >> /tmp/a.txt
##    cut -f 1 /tmp/a.txt | sort | uniq | wc -l
##
##done

D="samples/7046-004-019_discordant.fixed.bam"
samtools view "$D" "chr21:39,722,880-40,070,424" | grep "chr21,42" > /tmp/a.txt
samtools view "$D" "chr21:42,834,478-42,882,085" | grep -P "chr21,39|chr21,40" >> /tmp/a.txt
cat /tmp/a.txt

