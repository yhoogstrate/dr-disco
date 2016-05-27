#!/usr/bin/env python

import pysam,math



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

def valid_insert_size(read,max_size):
	if read.get_tag('RG') == 'discordant_mates':
		#print read.get_tag('SA')
		if read.reference_id == read.next_reference_id:
			dist = abs(read.reference_start - read.next_reference_start)
			#print read
			#print "template length",read.template_length
			#print "abs dist caclulated",
			#print
			
			if dist < max_size:
				return False
	return True

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
			
			passing = True
			
			# List of filters:
			#if r.get_tag('RG') == "silent_mate":
				#passing = False
				#break
			
			if not valid_insert_size(r,126):
				passing = False
			
			#if r.get_tag('RG') == "discordant_reads":
				## check for abs(ins size) > 126
				# passing = False
			
			if passing:
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
	
	exp_table("dump.txt",["s","e","j"],[vector_s,vector_e,vector_j])


def get_entropy(reads,nreferences):
	n = len(reads)
	if n <= 5:
		return None
	else:
		frequency_table = [0] * nreferences
		
		for read in reads:
			#print read
			#print read.reference_id
			frequency_table[read.next_reference_id] += 1
		prob = [float(x)/n for x in frequency_table]
		
		#print frequency_table
		#import sys
		#sys.exit()
		
		#frequency_table = [5] * nreferences
		#n = sum(frequency_table)
		#prob = [float(x)/n for x in frequency_table]
		
		entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob if p > 0 ])
		#entropy = max(0,entropy)# stupid -0
		
		# For plotting it's easier to normalize it to [0 <-> 1]
		entropy = entropy / (math.log(nreferences) / math.log(2.0))
		
		return entropy

def get_vector_type_2(region,bam,window_size):
	"""
	Vector for windowed entropy
	"""
	
	fh = open("entropy.txt","w")
	
	n = region[2] - region[1] + 1
	matrix = [[0,0,0]] * n
	vector_s = [0] * n
	vector_e = [0] * n
	
	# splice junction targets
	vector_j = [0] * n
	
	previous_count = 0
	for i in range(region[1],region[2],window_size):#+1? wi
		j = min(region[2],i+window_size)
		
		#window = i,j
		reads = []
		for r in bam.fetch(region[0],i,j):
			covered_bases = r.get_overlap(i,j)
			if covered_bases > 0:
				reads.append(r)
		
		entropy = get_entropy(reads,bam.nreferences)
		
		if entropy != None:
			fh.write(str(i-region[1])+"\t"+str(j-region[1])+"\t"+str(entropy)+"\n")
		else:
			fh.write(str(i-region[1])+"\t"+str(j-region[1])+"\t-1\n")
	
	fh.close()




#r1 = ['chr21',42851839,42851841]
#r1 = ['chr21',42851839,42852839]
r1 = ['chr21',42834187,42882196]
r1 = ['chr21',39736266,39928045]#erg
#r1 = ['chr21',39813226,39862468]#erg-testing
bam_file_discordant = "samples/7046-004-043_discordant.fixed.bam"
sam_file_discordant = pysam.AlignmentFile(bam_file_discordant, "rb")

get_vector_type_1(r1,sam_file_discordant)
get_vector_type_2(r1,sam_file_discordant,800)

#find_breakpoints(r1_tmprss,r2_erg,)
