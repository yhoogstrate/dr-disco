#!/usr/bin/env python

import pysam,math

windowsize = 800

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


def exp_table(fh,column_names,column_data):
	#fh = open(table,"w")
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

def get_vector_type_1(fh,fh_windowed,region,bam,window_size):
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
	#for i in range(region[1],region[2]+1):
	for r in bam.fetch(region[0],region[1],region[2]+1):
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
						vector_j[end+1] += 1
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
				vector_s[start+1] += 1
			if end < n:
				vector_e[end+1] += 1
	
	exp_table(fh,["s","e","j"],[vector_s,vector_e,vector_j])
	
	#total = [ vector_s[i]+vector_e[i] for i in range(len(vector_s))]
	#export_window(fh_windowed,total,window_size)

def export_window(fh_windowed,vec,size):
	n = len(vec)
	j = -1
	
	for i in range(n):
		if i % size == 0:
			if j >= 0:
				fh_windowed.write(str(j*size)+ "\t" + str((j+1)*size-1) + "\t" + str(c)+"\n")
			c = 0
			j += 1
		
		c += vec[i]
	
	if j >= 0:
		fh_windowed.write(str(j*size)+ "\t" + str((j+1)*size-1) + "\t" + str(c)+"\n")


#make_window([1,1, 0,0, 5,6, 7,8, 9,0, 1],3)

#import sys
#sys.exit()


def get_entropy(reads,nreferences):
	n = len(reads)
	if n <= 5:
		return None
	else:
		frequency_table = [0] * nreferences
		
		for read in reads:
			frequency_table[read.next_reference_id] += 1
		prob = [float(x)/n for x in frequency_table]
		
		#frequency_table = [5] * nreferences
		#n = sum(frequency_table)
		#prob = [float(x)/n for x in frequency_table]
		
		entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob if p > 0 ])
		#entropy = max(0,entropy)# stupid -0
		
		# For plotting it's easier to normalize it to [0 <-> 1]
		entropy = entropy / (math.log(nreferences) / math.log(2.0))
		
		return entropy

def get_entropy_mates(reads,references):
	#n = 
	nreferences = len(references)
	if len(reads) <= 5:
		return None
	else:
		frequency_table = [0] * nreferences
		
		for read in reads:
			chrs = [references.index(x.split(",",2)[0]) for x in read.get_tag('SA').split(";")]
			for _chr in chrs:
				frequency_table[_chr] += 1
			
		prob = [float(x)/sum(frequency_table) for x in frequency_table]
		
		entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob if p > 0 ])
		#entropy = max(0,entropy)# stupid -0
		
		# For plotting it's easier to normalize it to [0 <-> 1]
		entropy = entropy / (math.log(nreferences) / math.log(2.0))
		
		return entropy

def get_vector_type_2(fh,region,bam,window_size):
	"""
	Vector for sliding window entropy
	"""
	
	n = region[2] - region[1] + 1
	
	# splice junction targets
	vector_j = [0] * n
	offset = int(round(0.5 * window_size))
	
	previous_count = 0
	for k in range(region[1]+offset,region[2]-offset+1):#+1? w
		i = k - 200
		j = min(region[2],i+window_size)
		
		#window = i,j
		reads = []
		for r in bam.fetch(region[0],i,j):
			covered_bases = r.get_overlap(i,j)
			if covered_bases > 0:
				reads.append(r)
		
		entropy = get_entropy(reads,bam.nreferences)
		#entropy_m = get_entropy_mates(reads,bam.references)
		
		if entropy == None:
			entropy = -1
		#if entropy_m == None:
		#	entropy_m = -1
		
		#if entropy != None:
		fh.write(str(k-region[1])+"\t"+str(entropy)+"\n")
	
	fh.close()

def get_vector_type_3(fh,region,bam,window_size):
	"""
	Exon mask
	
	11.8s
	"""
	
	n = region[2] - region[1] + 1
	hashm = {}
	
	"""
	for i in range(region[1],region[2]+1):#+1? w
		j = i+1
		for r in bam.fetch(region[0],i,j):
			covered_bases = r.get_overlap(i,j)
			if covered_bases > 0:
				if r.cigarstring.find("N") != -1:
					hashm[i-region[1]] = 0
	
	"""
	
	for r in bam.fetch(region[0],region[1],region[2]+1 ):
		if r.cigarstring.find("N") != -1:
			length = r.reference_start
			for cigar in r.cigartuples:
				if cigar[0] in [0,2,3,7,8]:
					length += cigar[1]
			
			window_start = r.reference_start
			window_end = window_start + length
			
			if window_start < region[1]:
				window_start = region[1]
			
			if window_end > region[2]:
				window_end = region[2]
			
			for i in range(window_start, window_end):
				covered_bases = r.get_overlap(i,i+1)
				
				if covered_bases > 0:
					k = i -region[1] + 1
					hashm[k] = 0
	
	
	for i in range(n):
		if i in hashm.keys():
			fh.write(str(i)+"\t0\n")
		else:
			fh.write(str(i)+"\t1\n")
	
	fh.close()





#get_vector_type_1(open("test-test-vec1_a.tabular.txt","w"),open("test-test-vec1_b.tabular.txt","w"),r1,sam_file_discordant,windowsize)
#get_vector_type_2(open("test-test-vec2.tabular.txt","w"),r1,sam_file_discordant,windowsize)
#get_vector_type_3(open("test-test-vec3.tabular.txt","w"),r1,sam_file_discordant,windowsize)





samples = ['7046-004-041','7046-004-043']

genes = {}
genes['erg'] =     ['chr21',39737183,40035618]
genes['tmprss2'] = ['chr21',42834678,42882085]
genes['ar'] =      ['chrX' ,66761874,66952461]
genes['sash1'] =   ['chr6',148661729,148875184]
genes['samd5'] =   ['chr6',147827828,147893157]
genes['klk3'] =    ['chr19',51356171,51366020]
genes['gapdh'] =   ['chr12',6641585,6649537]



for sample in samples:
	print "Running sample: "+sample
	
	bam_file_discordant = "samples/"+sample+"_discordant.fixed.bam"
	sam_file_discordant = pysam.AlignmentFile(bam_file_discordant, "rb")
	
	for gene in genes.keys():
		print " - Obtanining data from: "+gene
		
		get_vector_type_1(open("data/"+sample+"-"+gene+"-vec1_a.tabular.txt","w"),None,genes[gene],sam_file_discordant,windowsize)
		get_vector_type_2(open("data/"+sample+"-"+gene+"-vec2.tabular.txt","w"),genes[gene],sam_file_discordant,windowsize)
		get_vector_type_3(open("data/"+sample+"-"+gene+"-vec3.tabular.txt","w"),genes[gene],sam_file_discordant,windowsize)




