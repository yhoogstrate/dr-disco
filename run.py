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

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

def exp_vector(name,data,fh):
	fh.write(name+"=c(")
	start = True
	for chunk in chunker(data,100):
		if start:
			fh.write(",".join([str(x) for x in chunk]))
			start = False
		else:
			fh.write(",\n\t"+",".join([str(x) for x in chunk]))
	fh.write(")\n")


def exp_table(table,column_names,column_data):
	fh = open(table,"w")
	x = len(column_data)
	y = len(column_data[0])
	fh.write("\t".join(column_names)+"\n")
	for i in range(y):
		#for j in range(x):
		column = [str(column_data[j][i]) for j in range(x)]
		fh.write("\t".join(column)+"\n")
	fh.close()

def get_vector_type_1(region,bam):
	"""
	Vector for disco alignments
	"""
	
	n = region[2] - region[1] + 1
	matrix = [[0,0,0]] * n
	vector_s = [0] * n
	vector_e = [0] * n
	
	# splice junction targets
	vector_j = [0] * n
	
	previous_count = 0
	for i in range(region[1],region[2]+1):
		for r in bam.fetch(region[0],i,i+1):
			start = r.reference_start - region[1]
			end = start
			for cigar in r.cigartuples:
				if cigar[0] == 3:
					if end < n and end > 0:
						vector_j[end] += 1
				
				if cigar[0] in [0,2,3,7,8]:
					end += cigar[1]
				
				if cigar[0] == 3:
					if end < n:
						vector_j[end] += 1
				"""
				cigar_m = r.get_cigar_stats()[0][0]
				cigar_i = r.get_cigar_stats()[0][1]#insertion has no weight
				cigar_d = r.get_cigar_stats()[0][2]
				cigar_n = r.get_cigar_stats()[0][3]
				cigar_s = r.get_cigar_stats()[0][4]
				cigar_h = r.get_cigar_stats()[0][5]
				#cigar_p = r.get_cigar_stats()[0][6]
				cigar_S = r.get_cigar_stats()[0][7]
				cigar_x = r.get_cigar_stats()[0][8]
				"""
			
			#end = start + cigar_m + cigar_d + cigar_n + cigar_S + cigar_x
			# - cigar_s - cigar_h
			
			if start >= 0:# spliced alignments may start before the start point
				vector_s[start] += 1
			if end < n:
				vector_e[end] += 1
			
			#import sys
			#sys.exit(1)
	
	exp_table("dump.txt",["s","e","j"],[vector_s,vector_e,vector_j])
	#fh = open("dump.R","w")
	#exp_vector("s",vector_s,fh)
	#exp_vector("e",vector_e,fh)
	#exp_vector("j",vector_j,fh)
	#fh.close()
	
	#plot(1:length(s),s,pch=20,cex=0.2,col="blue")
	##points(1:length(e),s,pch=20,cex=0.2,col="green")
	#points(1:length(j),s,pch=20,cex=0.2,col="red")

#def find_breakpoints(r1,r2, bam_file_concordant, bam_file_discordant):
	#sam_file_discordant = pysam.AlignmentFile(bam_file_discordant, "rb")
#reads = []


#tmprss2=['chr21',42833093,42891172]

#erg_bp =['chr21', 39859063, 39859285]
#erg_bp =['chr21', 39859069, 39859271]

#for read in sam_file_discordant.fetch(erg_bp[0],erg_bp[1],erg_bp[2]):
	#reads.append(read)
	
	#readname = "D00476:156:C6VVJANXX:7:2115:5471:96808"
	#for read in reads:
		#if read.qname == readname:
			#print read.qname
			#print "\t",read.get_blocks()
			
			#if read.is_paired:
				#print "\t*",read.next_reference_start
				#print "\t-",sam_file_discordant.mate(read)
				#print "\tb",read.get_blocks()
				##print "\tn",read.next_reference_id
				##print "\tx",read.get_aligned_pairs()
	
	#sam_file_discordant.close()


#hg19
#r1_tmprss = "chr21:42,834,478-42,882,085"
#r2_erg = "chr21:39,737,183-40,035,618"
#offset = 5000

#r1 = ['chr21',42851839,42851841]
#r1 = ['chr21',42851839,42852839]
r1 = ['chr21',42834187,42882196]
r1 = ['chr21',39737183,40035618]#erg
bam_file_discordant = "samples/7046-004-043_discordant.bam"
sam_file_discordant = pysam.AlignmentFile(bam_file_discordant, "rb")

get_vector_type_1(r1,sam_file_discordant)

#find_breakpoints(r1_tmprss,r2_erg,)
