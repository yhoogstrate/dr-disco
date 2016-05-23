#!/usr/bin/env python

import pysam
#numpy



## Sep 01: fix bam files:
## - SA tag at the end to link the multiple supplementary alignments IF present - by using samtools -n 
## - Fix the mate of the other: 
##  [     A   >                [  B ]   [  C >
##  A -> B
##  B -> C
##  C <- B
## - Update NH tags per mate
## - Mark as deletion if possible (N is for splicing, so do D/P)
## - Set read group to sample name


def find_breakpoints(r1,r2, bam_file_concordant, bam_file_discordant):
	sam_file_discordant = pysam.AlignmentFile(bam_file_discordant, "rb")
reads = []

erg_bp =['chr21', 39859063, 39859285]
erg_bp =['chr21', 39859069, 39859271]

for read in sam_file_discordant.fetch(erg_bp[0],erg_bp[1],erg_bp[2]):
	reads.append(read)
	
	readname = "D00476:156:C6VVJANXX:7:2115:5471:96808"
	for read in reads:
		if read.qname == readname:
			print read.qname
			print "\t",read.get_blocks()
			
			if read.is_paired:
				print "\t*",read.next_reference_start
				print "\t-",sam_file_discordant.mate(read)
				print "\tb",read.get_blocks()
				#print "\tn",read.next_reference_id
				#print "\tx",read.get_aligned_pairs()
	
	sam_file_discordant.close()


#hg19
r1_tmprss = "chr21:42,834,478-42,882,085"
r2_erg = "chr21:39,737,183-40,035,618"
offset = 5000

find_breakpoints(r1_tmprss,r2_erg,"samples/7046-004-041_concordant.bam","samples/7046-004-041_discordant.bam")
