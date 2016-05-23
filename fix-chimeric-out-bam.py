#!/usr/bin/env python

#http://www.samformat.info/sam-format-flag

import pysam
#numpy

bam_file_discordant = "samples/7046-004-041_discordant.n.bam"

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

def fix_chain(alignments):
	"""
	aligments must come from the same MATE and reads should not be marked as DUPLICATE or MULTIMAPPING
	"""
	chains_from = {}
	chains_to = {}
	
	# Fill phase
	for alignment in alignments:
		chain_from = alignment.reference_name+":"+str(alignment.reference_start)
		chain_to = alignment.next_reference_name+":"+str(alignment.next_reference_start)
		
		print chain_from + " => " + chain_to
		
		if chain_from not in chains_from.keys():
			chains_from[chain_from] = []
		if chain_to not in chains_to.keys():
			chains_to[chain_to] = []
		
		chains_from[chain_from].append(alignment)
		chains_to[chain_to].append(alignment)
	
	chains_from = set(chains_from)
	chains_to = set(chains_to)
	
	#links = []
	#for chain_to in chains_to.keys():
		#if chain_to in chains_from.keys():
			#alignment_t = chains_from[chain_to][0]
			#chain_to_t = alignment_t.next_reference_name+":"+str(alignment_t.next_reference_start)
			#print "link: "+chain_to+"=>"+chain_to_t
	
	_from = chains_from.difference(chains_to)
	_to = chains_to.difference(chains_from)
	_linked = list(chains_to.intersection(chains_from))
	
	# situaties:
	# 1. net zoveel chains_to als chains_from -> alleen unieke start posities
	#    a. alle chains_from en chains_to zijn identiek, op 1 in beide lijsten na
	#       [a,b,c,d,e,f,g,h]
	#         [b,c,d,e,f,g,h,i]
	#       In dit geval is de chain al perfect
	#    b. de overlap van chains heeft meer dan 1 mismatch per vector - er is iets
	#       mis met de mapping van de next reads. Dit klopt niet en spuug een error
	#       uit
	#        [a,b,c,c,d,e]
	#          [b,c,c,d,e,f]
	# 2. Er zijn chains met dubbele entries in de TO en dus ook minder die overeen
	#    komen met de chains from. Updaten moet door die to's aan te passen. De
	#    entrie die het dichtste bij de 'to' zit moet er naar blijven linken, de
	#    andere(n) naar elkaar - per chromosoom?
	
	if len(list(_from)) == 1 and len(list(_to)) == 1:
		# Chain is seems to be OKAY - only add SA suffixes/tags
		pass
	else:
		print
		for alignment in alignments:
			print alignment
		print
		print "f",chains_from
		print "t",chains_to
		
		if len(_linked) == 0:
			#They're all mapping to an element not in this list
			# Just find the one closest to the chain_to and continue as if you have a start point
			pass
		else:
			#There is at least a start of linked alignments
			#print "intersect",_linked
			pass
		
		# continue as if you have a start point

"""
f set(['chr9:100772404', 'chr9:100772393'])
t set(['chr9:100772326'])

100772404-100772326=78
100772326-100772393=-67
"""



def reconstruct_alignments(alignments):
	n = len(alignments)
	r1 = []
	r2 = []
	singletons = []
	
	
	for alignment in alignments:
		if alignment.is_read1:
			r1.append(alignment)
		elif alignment.is_read2:
			r2.append(alignment)
		else:
			singletons.append(alignment)
	
	n_r1 = len(r1)
	n_r2 = len(r2)
	n_s  = len(singletons)
	
	double_disco = 0
	
	if n_r1 > 1:
		## set "SA:S:chrA:1:10;,,," tag and change the next alignment
		double_disco += 1
		fix_chain(r1)
		print "---"
	
	if n_r2 > 1:
		## set "SA:S:chrA:1:10;,,," tag and change the next alignment
		double_disco += 1
		#for alignment in r2:
		#	print alignment.next_reference_name+":"+str(alignment.next_reference_start)
		#	print
	
	#if n_r1 > 1 or n_r2 > 1:
	#	print "--------------------------------"
	
	if n_s >= 2:
		## Only set "SA:S:chrA:1:10;,,," tag and do not change the next alignment
		#print singletons[0]
		#print singletons[1]
		#print
		#import sys
		#sys.exit(1)
		pass
	
	#if double_disco >= 2:
	#	print "DOUBLE DISCO!"
	#	import sys
	#	sys.exit(1)
	
	#if n > 2:
	#print n,n_r1, n_r2, n_s
	#   9038 2 0 0 2
	# 122063 2 1 1 0
	#  12689 3 1 2 0
	#  15828 3 2 1 0


	

sam_file_discordant = pysam.AlignmentFile(bam_file_discordant, "rb")
last_read_name = False
alignments = []
for read in sam_file_discordant:
	if read.qname != last_read_name:
		if len(alignments) > 0:
			reconstruct_alignments(alignments)
		alignments = []
		last_read_name = read.qname
	alignments.append(read)
reconstruct_alignments(alignments)



# 1. Check name sorted
# 2. 


#hg19
#r1_tmprss = "chr21:42,834,478-42,882,085"
#r2_erg = "chr21:39,737,183-40,035,618"
#offset = 5000

#find_breakpoints(r1_tmprss,r2_erg,"samples/7046-004-041_concordant.bam",)
