#!/bin/sh

sort -k 1.4,1 -s share/blacklist-regions.hg38.bed > tmp/a.bed
mv tmp/a.bed share/blacklist-regions.hg38.bed

sort -k 1.4,1 -s share/blacklist-junctions.hg38.txt > tmp/b.txt
mv tmp/b.txt share/blacklist-junctions.hg38.txt

