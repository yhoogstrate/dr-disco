#!/usr/bin/env python

#http://www.samformat.info/sam-format-flag

import click,os,subprocess
import pysam

if pysam.__version__[0:4] != "0.9.":
	raise Exception("Version of pysam needs to be at least 0.9 but is: "+pysam.__version__+" instead")


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

def set_qname_to_group(all_reads_updated):
	qnames = []
	for a in all_reads_updated:
		qnames.append(a.qname)
	
	qnames = list(set(qnames))
	
	if len(qnames) != 1:
		raise Exception("Not all reads belong to the same QNAME")
	else:
		qname = qnames[0]
		for i in range(len(all_reads_updated)):
			all_reads_updated[i].set_tag('RG',qname)
	
	
	return all_reads_updated

def update_sa_tags(reads_updated,bam_file):
	sa_ids = []
	
	for a in reads_updated:
		try:
			nm = a.get_tag('nM')
		except:
			nm = -1
		
		sa_id = [bam_file.get_reference_name(a.reference_id),a.reference_start,a.cigarstring,a.mapping_quality,"-" if a.is_reverse else "+",nm]
		sa_ids.append(",".join([str(x) for x in sa_id]))
	
	for i in range(len(reads_updated)):
		sa_id = sa_ids[i]
		sa_tag = ";".join([x for x in sa_ids if x != sa_id])
		reads_updated[i].set_tag('SA',sa_tag)
	
	return reads_updated

def set_next_ref(aligned_segment,position):
	if aligned_segment.next_reference_id != position[0] or aligned_segment.next_reference_start != position[1]:
		a = pysam.AlignedSegment()
		a.query_name = aligned_segment.query_name
		a.query_sequence = aligned_segment.query_sequence
		a.flag = aligned_segment.flag
		a.reference_id = aligned_segment.reference_id
		a.reference_start = aligned_segment.reference_start
		a.mapping_quality = aligned_segment.mapping_quality
		a.cigartuples = aligned_segment.cigartuples
		a.template_length = aligned_segment.template_length
		a.query_qualities = aligned_segment.query_qualities
		a.set_tags(aligned_segment.get_tags())
		
		a.next_reference_id = position[0]
		a.next_reference_start = position[1]
		
		return a
	else:
		return aligned_segment

def get_closest(location, alignments):
	d_chr = None
	d_pos = None
	closest = None
	
	for alignment in alignments:
		is_closest = False
		if (d_chr == None and d_pos == None):
			is_closest = True
		else:
			dd_chr = abs(d_chr - location[0])
			dd_pos = abs(d_pos - location[1])
			
			if (dd_chr < d_chr) or (dd_chr == d_chr and dd_pos < d_pos):
				# Ref to itself is ackward
				#if dd_chr == 0 and dd_pos == 0:
				#	pass
				#else:
				is_closest = True
		
		if is_closest:
			closest = alignment
	
	return closest

def fix_chain(alignments,bam_file,mates):
	"""
	all segments in variable `aligments` must come from the same MATE and should not be marked as DUPLICATE or MULTIMAPPING
	"""
	chains_from = {}
	chains_to = {}
	
	k = len(alignments)
	
	# Fill phase
	for alignment in alignments:
		chain_from = str(alignment.reference_id)+":"+str(alignment.reference_start)
		chain_to = str(alignment.next_reference_id)+":"+str(alignment.next_reference_start)
		
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
	
	new_alignments = []
	
	if len(list(_from)) == 1 and len(list(_to)) == 1:
		# Chain is seems to be OKAY - only add SA suffixes/tags
		new_alignments = alignments
	else:
		if len(_linked) == 0:
			if len(chains_to) > 1:
				# Probably pick simply the link to the mate as follows:
				# mate_pos = [mates[0].reference_id,mates[0].reference_start]
				raise Exception("Error - unknown situation type 1")
			
			mate_pos = [int(x) for x in list(chains_to)[0].split(":")]
			if mate_pos[0] == -1 and mate_pos[1] == -1:
				mate_pos = [alignments[0].reference_id,alignments[0].reference_start]
			
			closest = get_closest(mate_pos,alignments)
			alignments = [a for a in alignments if a != closest]
			new_alignments.append(set_next_ref(closest,mate_pos))
			
			while len(alignments) > 0:
				seg_pos = [closest.reference_id,closest.reference_start]
				closest = get_closest(mate_pos,alignments)
				
				alignments = [a for a in alignments if a != closest]
				new_alignments.append(set_next_ref(closest,seg_pos))
			
		else:
			#There is at least a start of linked alignments
			if len(_linked) > 1:
				if len(_from) == 0 and len(_to) == 0:
					# cross reffing each other
					alignments_new = alignments
				else:
					raise Exception("Error - unknown situation type 2")
			
			if len(mates) == 1:
				seg_pos = [mates[0].reference_id,mates[0].reference_start]
			else:
				seg_pos = [int(x) for x in list(_linked)[0].split(":")]
				if seg_pos[0] == -1 and seg_pos[1] == -1:
					seg_pos = [alignments[0].reference_id,alignments[0].reference_start]
			
			closest = get_closest(seg_pos,alignments)
			alignments = [a for a in alignments if a != closest]
			new_alignments.append(set_next_ref(closest,seg_pos))
			
			while len(alignments) > 0:
				seg_pos = [closest.reference_id,closest.reference_start]
				closest = get_closest(seg_pos,alignments)
				
				alignments = [a for a in alignments if a != closest]
				new_alignments.append(set_next_ref(closest,seg_pos))
	
	if len(new_alignments) != k:
		raise Exception("Somewhere alignments got lost in this function")
	
	return new_alignments



def reconstruct_alignments(alignments,bam_file,fh_out):
	"""
	input: only reads with the same qname
	"""
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
	
	all_reads_updated = []
	
	if n_r1 > 1:
		double_disco += 1
		
		reads_updated = fix_chain(r1,bam_file,r2)
		#reads_updated = update_sa_tags(reads_updated,bam_file)
		for a in reads_updated:
			all_reads_updated.append(a)
	else:
		# No need to update - no supplementary r1 reads
		for a in r1:
			all_reads_updated.append(a)
	
	if n_r2 > 1:
		double_disco += 1
		
		reads_updated = fix_chain(r2,bam_file,r1)
		#reads_updated = update_sa_tags(reads_updated,bam_file)
		for a in reads_updated:
			all_reads_updated.append(a)
	else:
		# No need to update - no supplementary r1 reads
		for a in r2:
			all_reads_updated.append(a)
	
	if n_s >= 2:
		double_disco += 1
		
		reads_updated = fix_chain(singletons,bam_file,[])
		#reads_updated = update_sa_tags(reads_updated,bam_file)
		for a in reads_updated:
			all_reads_updated.append(a)
	else:
		# No need to update - no supplementary singleton reads
		for a in singletons:
			all_reads_updated.append(a)
	
	if double_disco >= 2:
		# Should in theory be possible but has not been observed
		raise Exception("DOUBLE DISCO! >> forward and reverse have multiple segments?")
	
	all_reads_updated = update_sa_tags(all_reads_updated,bam_file)
	if len(all_reads_updated) != n:
		raise Exception("Error - reads have been lost")
	
	all_reads_updated = set_qname_to_group(all_reads_updated)
	
	for a in all_reads_updated:
		fh_out.write(a.tostring(bam_file)+"\n")
	
	#if n > 2:
	#print n,n_r1, n_r2, n_s
	#   9038 2 0 0 2
	# 122063 2 1 1 0
	#  12689 3 1 2 0
	#  15828 3 2 1 0


@click.command(help="This tool requires the '*.Chimeric.out.sam' files of RNA STAR converted into BAM")
@click.argument('bam_file_discordant')
@click.argument('bam_file_discordant_fixed')
def main(bam_file_discordant,bam_file_discordant_fixed):
	path = os.path.dirname(bam_file_discordant)
	basename,ext = os.path.splitext(os.path.basename(bam_file_discordant))
	
	#@TODO / consider todo - start straight from sam
	#samtools view -bS samples/7046-004-041_discordant.Chimeric.out.sam > samples/7046-004-041_discordant.Chimeric.out.unsorted.bam
	
	print "Convert into a name-sorted bam file (to get all reads with the same name adjacent to each other"
	command = ["samtools",
			   "sort",
			   "-n",
			   bam_file_discordant,
			   basename+".name-sorted"]
	e_code = subprocess.call(command)
	
	if e_code != 0:
		raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
	
	
	print "--- Fixing sam file ---"
	sam_file_discordant = pysam.AlignmentFile(basename+".name-sorted.bam", "rb")
	fh = open(basename+".name-sorted.fixed.sam","w")
	fh.write(sam_file_discordant.text)
	last_read_name = False
	alignments = []
	for read in sam_file_discordant:
		if read.qname != last_read_name:
			if len(alignments) > 0:
				reconstruct_alignments(alignments,sam_file_discordant,fh)
			alignments = []
			last_read_name = read.qname
		alignments.append(read)
	reconstruct_alignments(alignments,sam_file_discordant,fh)
	fh.close()
	
	
	print "Converting fixed file into BAM"
	command = ["samtools",
			   "view",
			   "-bS",
			   "-o",
			   basename+".name-sorted.fixed.bam",
			   basename+".name-sorted.fixed.sam"]
	e_code = subprocess.call(command)
	
	if e_code != 0:
		raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
	
	
	print "Sorting position based fixed file"
	command = ["samtools",
			   "sort",
			   basename+".name-sorted.fixed.bam",
			   basename+".sorted.fixed"]
	e_code = subprocess.call(command)
	
	if e_code != 0:
		raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
	
	
	print "Indexing the position sorted bam file"
	command = ["samtools",
			   "index",
			   basename+".sorted.fixed.bam"]
	e_code = subprocess.call(command)
	
	if e_code != 0:
		raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
	
	
	print "Cleaning up temp files"
	os.remove(basename+".name-sorted.bam")
	os.remove(basename+".name-sorted.fixed.sam")
	os.remove(basename+".name-sorted.fixed.bam")
	
	
	print "Moving to final destination"
	os.rename(basename+".sorted.fixed.bam",bam_file_discordant_fixed)
	os.rename(basename+".sorted.fixed.bam"+".bai",bam_file_discordant_fixed+".bai")


if __name__ == "__main__":
	main()
