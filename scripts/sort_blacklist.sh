#!/bin/sh

sort -k 1.4,1 -s share/blacklist-regions.hg38.bed > tmp/a.bed
mv tmp/a.bed share/blacklist-regions.hg38.bed

