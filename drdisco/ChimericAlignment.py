#!/usr/bin/env python3

#http://www.samformat.info/sam-format-flag

import click,os,subprocess
import pysam


#@todo check samtools v >= 1.3.*

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


class ChimericAlignment:
	def __init__(self,f1,f2,f3):
		self.test_pysam_version()
	
	def test_pysam_version(self):
		if pysam.__version__[0:4] != "0.9.":
			raise Exception("Version of pysam needs to be at least 0.9 but is: "+pysam.__version__+" instead")
		else:
			return True
	
	def set_read_group(self,all_reads_updated,group):
		for a in all_reads_updated:
			a.set_tag('RG',group)

	def fix_alignment_score(self,all_reads_updated):
		for a in all_reads_updated:
			a.is_paired = True
			a.is_read2 = True

	def set_qname_to_group(self,all_reads_updated):
		qnames = []
		for a in all_reads_updated:
			qnames.append(a.qname)
		
		qnames = list(set(qnames))
		
		if len(qnames) != 1:
			raise Exception("Not all reads belong to the same QNAME")
		else:
			qname = qnames[0]
			for i in range(len(all_reads_updated)):
				all_reads_updated[i].set_tag('LB',qname.replace(":","."))
		
		return all_reads_updated

	def update_sa_tags(self,reads_updated,bam_file):
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

	def set_next_ref(self,aligned_segment,position):
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

	def get_closest(self,location, alignments):
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

	def fix_chain(self,alignments,bam_file,mates):
		"""
		all segments in variable `aligments` must come from the same MATE (unless discordant pairs without suppl. alignments) and should not be marked as DUPLICATE or MULTIMAPPING
		"""
		chains_from = {}
		chains_to = {}
		
		k = len(alignments)
		
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
		new_mates = []
		
		if len(mates) == 1 and len(alignments) >= 2:
			if str(mates[0].next_reference_id)+":"+str(mates[0].next_reference_start) in chains_from:
				next_pos = [mates[0].next_reference_id,mates[0].next_reference_start]
				last_pos = [mates[0].reference_id,mates[0].reference_start]
				
				start = get_closest(next_pos,alignments)
			
			else:
				print("Warning - mates do not correspond? - maybe empty (-1) as well?")
				
				next_pos = [alignments[0].reference_id,alignments[0].reference_start]
				last_pos = [mates[0].reference_id,mates[0].reference_start]
				
				start = get_closest(next_pos,alignments)
			
			alignments = [a for a in alignments if a != start]
			# If the mate is not exactly matched but close, fix it:
			new_mates.append(set_next_ref(mates[0],[start.reference_id,start.reference_start]))
			
			i = 0
			while len(alignments) >= 1:
				closest = get_closest(next_pos,alignments)
				next_pos = [closest.reference_id,closest.reference_start]
				
				s_fixed = set_next_ref(start,next_pos)
				s_fixed.set_tag('FI',i)
				new_alignments.append(s_fixed)
				alignments = [a for a in alignments if a != closest]
				
				start = closest
				i += 1
			
			# Map last one back to the mate again
			if len(alignments) == 0:
				start = set_next_ref(start,last_pos)
				start.set_tag('FI',i)
				new_alignments.append(start)
		
		elif len(mates) == 0:
			## Either 2 discordant mates
			## Or 2 discordant segments from one singleton
			
			if len(_linked) == len(alignments):
				# cross reffing each other - is already fine
				return alignments,new_mates
			
			"""
			for a in alignments:
				print a
			print "f",_from
			print "t",_to
			print "l",_linked
			"""
			
			seg_pos = None
			if len(_linked) > 0:
				seg_pos = [int(x) for x in list(_linked)[0].split(":")]
			
			if seg_pos == None or (seg_pos[0] == -1 and seg_pos[1] == -1):
				seg_pos = [alignments[0].reference_id,alignments[0].reference_start]
			
			closest = get_closest(seg_pos,alignments)
			alignments = [a for a in alignments if a != closest]
			new_alignments.append(set_next_ref(closest,seg_pos))
			
			while len(alignments) > 0:
				seg_pos = [closest.reference_id,closest.reference_start]
				closest = get_closest(seg_pos,alignments)
				
				alignments = [a for a in alignments if a != closest]
				new_alignments.append(set_next_ref(closest,seg_pos))
		
		else:
			raise Exception("Dunno how to handle junctions in both mates yet... not aware of STAR Fusion producing them either")
		
		if len(new_alignments) != k:
			raise Exception("Somewhere alignments got lost in this function")
		
		return new_alignments,new_mates



	def reconstruct_alignments(self,alignments,bam_file,fh_out):
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
			reads_updated,mates_updated = fix_chain(r1,bam_file,r2)
			set_read_group(reads_updated,'spanning_paired')
			set_read_group(mates_updated,'silent_mate')
			for a in reads_updated:
				all_reads_updated.append(a)
			
			for a in mates_updated:
				all_reads_updated.append(a)
		
		elif n_r2 > 1:
			reads_updated,mates_updated = fix_chain(r2,bam_file,r1)
			set_read_group(reads_updated,'spanning_paired')
			set_read_group(mates_updated,'silent_mate')
			for a in reads_updated:
				all_reads_updated.append(a)
			
			for a in mates_updated:
				all_reads_updated.append(a)
		
		elif n_s >= 2:
			reads_updated,mates_updated = fix_chain(singletons,bam_file,[])
			set_read_group(reads_updated,'spanning_singleton')
			
			fix_alignment_score(reads_updated)
			
			for a in reads_updated:
				all_reads_updated.append(a)
			
			for a in mates_updated:
				all_reads_updated.append(a)
		
		elif n_r1 == 1 and n_r2 == 1 and n_s == 0:
			reads_updated,mates_updated = fix_chain(r1 + r2,bam_file,[])
			set_read_group(reads_updated,'discordant_mates')
			
			for a in reads_updated:
				all_reads_updated.append(a)
			
			for a in mates_updated:
				all_reads_updated.append(a)
		
		else:
			if n == 1:
				print("Warning: segments of mate are missing: "+alignments[0].query_name)
				all_reads_updated.append(alignments[0])
			else:
				raise Exception("what happens here?")
		
		
		all_reads_updated = update_sa_tags(all_reads_updated,bam_file)
		if len(all_reads_updated) != n:
			raise Exception("Error - reads have been lost")
		
		all_reads_updated = set_qname_to_group(all_reads_updated)
		
		for a in all_reads_updated:
			fh_out.write(a)
		
		#if n > 2:
		#print n,n_r1, n_r2, n_s
		#   9038 2 0 0 2
		# 122063 2 1 1 0
		#  12689 3 1 2 0
		#  15828 3 2 1 0


	def convert(self,bam_file_discordant,bam_file_discordant_fixed):
		path = os.path.dirname(bam_file_discordant)
		basename,ext = os.path.splitext(os.path.basename(bam_file_discordant))
		
		#@TODO / consider todo - start straight from sam
		#samtools view -bS samples/7046-004-041_discordant.Chimeric.out.sam > samples/7046-004-041_discordant.Chimeric.out.unsorted.bam
		
		print("Convert into a name-sorted bam file, to get all reads with the same name adjacent to each other")
		command = ["samtools",
				   "sort",
				   "-o",basename+".name-sorted.bam",
				   "-n",
				   bam_file_discordant]
		e_code = subprocess.call(command)
		
		if e_code != 0:
			raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
		
		
		print("--- Fixing sam file ---")
		sam_file_discordant = pysam.AlignmentFile(basename+".name-sorted.bam", "rb")
		header = sam_file_discordant.header
		header['RG'] = [
			{'ID':'spanning_singleton','DS':'This read was aligned to two locations but no aligned mate'},
			{'ID':'discordant_mates','DS':'This read has discordant mate pair'},
			{'ID':'spanning_paired','DS':'This read was aligned to two locations and also has an aligned mate'},
			{'ID':'silent_mate','DS':'Reads of this type are not discordant while their mate is'}]
		
		fh = pysam.AlignmentFile(basename+".name-sorted.fixed.sam", "wb", header=header)
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
		
		
		print("Converting fixed file into BAM")
		command = ["samtools",
				   "view",
				   "-bS",
				   "-o",
				   basename+".name-sorted.fixed.bam",
				   basename+".name-sorted.fixed.sam"]
		e_code = subprocess.call(command)
		
		if e_code != 0:
			raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
		
		
		print("Sorting position based fixed file")
		## Samtools 1.3.1
		command = ["samtools",
				   "sort",
				   "-o",basename+".sorted.fixed.bam",
				   basename+".name-sorted.fixed.bam"]
		e_code = subprocess.call(command)
		
		if e_code != 0:
			raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
		
		
		print("Indexing the position sorted bam file")
		command = ["samtools",
				   "index",
				   basename+".sorted.fixed.bam"]
		e_code = subprocess.call(command)
		
		if e_code != 0:
			raise Exception("Abnormal termination of samtools: exit="+str(e_code)+" while running:\n"+"\t".join(command))
		
		
		print("Cleaning up temp files")
		os.remove(basename+".name-sorted.bam")
		os.remove(basename+".name-sorted.fixed.sam")
		os.remove(basename+".name-sorted.fixed.bam")
		
		
		print("Moving to final destination")
		os.rename(basename+".sorted.fixed.bam",bam_file_discordant_fixed)
		os.rename(basename+".sorted.fixed.bam"+".bai",bam_file_discordant_fixed+".bai")

