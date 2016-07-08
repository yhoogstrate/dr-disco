#!/usr/bin/env python3
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Tries to figure out within a discordant RNA-Seq read alignment:
 - Tries to explain introns and exons by fusion transcripts
"""

#http://www.samformat.info/sam-format-flag

import logging,re
import pysam

from fuma.Fusion import STRAND_FORWARD, STRAND_REVERSE, STRAND_UNDETERMINED 

class Arc:
    def __init__(self,_target):
        self._target = _target
        self._types = {}
    
    def add_type(self,_type):
        if not self._types.has_key(_type):
            self._types[_type] = 0
        
        self._types[_type] += 1

class Node:
    def __init__(self,position):
        self.position = position
        self.arcs = {}
    
    def insert_arc(self,arc):
        skey = str(arc._target)
        
        if not self.arcs.has_key(skey):
            self.arcs[skey] = Arc(skey)
        
        self.arcs[skey] = Arc
    
    def add_arc(self,node2,arc_type,do_vice_versa=True):
        if do_vice_versa:
            node2.add_arc(self,arc_type,False)
        
        arc = Arc(node2)
        arc.add_type(arc_type)
        self.insert_arc(arc)
        #print "adding arc",str(self.position)," -> ",str(node2.position)

def bam_parse_alignment_end(read):
    pos = read.reference_start
    if not read.is_reverse:
        for chunk in read.cigar:
            """ M	BAM_CMATCH	0
                I	BAM_CINS	1
                D	BAM_CDEL	2
                N	BAM_CREF_SKIP	3
                S	BAM_CSOFT_CLIP	4
                H	BAM_CHARD_CLIP	5
                P	BAM_CPAD	6
                =	BAM_CEQUAL	7
                X	BAM_CDIFF	8"""
            
            if chunk[0] in [0,2,3]:
                pos += chunk[1]
    
    return pos

pat_bam_parse_alignment_offset_using_cigar = re.compile("([0-9]+)([MIDNSHPX=])")
def bam_parse_alignment_offset_using_cigar(sa_tag):
    """Parses the offset to theend point of the reads' mate
    """
    pos = 0
    if sa_tag[4] == '+':
        for chunk in pat_bam_parse_alignment_offset_using_cigar.finditer(sa_tag[2]):
            """ M	BAM_CMATCH	0
                I	BAM_CINS	1
                D	BAM_CDEL	2
                N	BAM_CREF_SKIP	3
                S	BAM_CSOFT_CLIP	4
                H	BAM_CHARD_CLIP	5
                P	BAM_CPAD	6
                =	BAM_CEQUAL	7
                X	BAM_CDIFF	8"""
            
            if chunk.group(2) in ['M','D','N']:
                pos += int(chunk.group(1))
    
    return pos

class BreakPosition:
    """
    Unambiguous data type that determines where a break point occurs.
    
    chr1:
    | 0 | 1 | 2 | 3 | ~~~ 
    
    chr2:
     ~~~ | 6 | 7 | 8 |
     
     string function of break points:
     chr1:3/4 and chr2:5/6
    """
    def __init__(self,_chr,position_0_based,strand):
        self._chr = _chr
        self.pos = position_0_based
        self.strand = strand
    
    def __str__(self):
        if self.strand == STRAND_FORWARD:
            return str(self._chr)+":"+str(self.pos)+"/"+str(self.pos+1)+"(+)"
        elif self.strand == STRAND_REVERSE:
            return str(self._chr)+":"+str(self.pos)+"/"+str(self.pos+1)+"(-)"
        else:    
            return str(self._chr)+":"+str(self.pos)+"/"+str(self.pos+1)+"(?)"


class Chain:
    def __init__(self,pysam_fh):
        self.idx = {}
        self.pysam_fh = pysam_fh
    
    def insert_entry(self,pos1,pos2,_type):
        """
         - Checks if Node exists at pos1, otherwise creates one
         - Checks if Node exists at pos2, otherwise creates one
         - Checks if Arc exists between them
         
        """
        
        if not self.idx.has_key(pos1):
            self.idx[pos1] = Node(pos1)
        
        if not self.idx.has_key(pos2):
            self.idx[pos2] = Node(pos2)
        
        self.idx[pos1].add_arc(self.idx[pos2],_type)
    
    def insert(self,read,parsed_SA_tag):
        """Inserts a read in the Chain and determine the type of arc"""
        # Type:
        # 1. 'N' <- N alignment flag in SAM, meaning INTRON
        # 2. 'discordant_read'
        # 3. 'split_read'
        # 4. 'silent_mate'
        #
        # ... 'S' and 'H' for soft and hard clipping?
        
        #print read.query_name
        #print " -" , read.reference_start, "-" if read.is_reverse else "+"
        #print " -" , read.cigarstring
        #print " -" , read.get_tag('RG')
        #print " -" , parsed_SA_tag
        
        if read.get_tag('RG') == "discordant_mates":
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),bam_parse_alignment_end(read), not read.is_reverse)
            
            if read.mate_is_reverse:
                pos2 = BreakPosition(parsed_SA_tag[0][0], parsed_SA_tag[0][1], STRAND_FORWARD if parsed_SA_tag[0][4] == "+" else STRAND_REVERSE)
            else:
                pos2 = BreakPosition(parsed_SA_tag[0][0], parsed_SA_tag[0][1] + bam_parse_alignment_offset_using_cigar(parsed_SA_tag[0]), STRAND_FORWARD if parsed_SA_tag[0][4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,"discordant_mates")
        
        self.find_sam_SHI_arcs(read)
    
    def find_sam_SHI_arcs(self,r):
        """Tries to find ARCs introduced by:
         - Hard clipping
         - Soft clipping
         - Splicing"""
        return True
    
    def prune(self):
        pass
        # do some clever tricks to merge arcs together and reduce data points



class IntronDecomposition:
    def __init__(self,break_point):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        self.break_point = break_point
    
    def test_disco_alignment(self,alignment_file):
        # Make sure the header exists indicating that the BAM file was
        # fixed using Dr. Disco
        bam_fh = pysam.AlignmentFile(alignment_file, "rb")
        if bam_fh.header.has_key('PG'):
            for pg in bam_fh.header['PG']:
                if pg['ID'] == 'drdisco_fix_chimeric':
                    return bam_fh
        raise Exception("Invalid STAR BAM File: has to be post processed with 'dr-disco fix-chimeric ...' first")

    def annotate_genes(self,gene_set):
        pass
    
    def is_exonic(self,read,gene):
        pass
    
    def get_insert_size(self,pos1,pos2):
        if pos1[0] == pos2[0]:
            return pos2[1] - pos1[1]
        else:
            return 99999999
    
    def parse_SA(self,SA_tag):
        sa_tags = SA_tag.split(":")
        for i in range(len(sa_tags)):
            sa_tags[i] = sa_tags[i].split(",")
            sa_tags[i][1] = int(sa_tags[i][1])
        
        return sa_tags
    
    def decompose(self,alignment_file):
        pysam_fh = self.test_disco_alignment(alignment_file)
        
        lpos = self.break_point.get_left_position(True)
        rpos = self.break_point.get_right_position(True)
        
        chain_left = []
        chain_right = []
        
        c = Chain(pysam_fh)
        
        for r in pysam_fh.fetch(lpos[0],lpos[1]):
            sa = self.parse_SA(r.get_tag('SA'))
            _chr = pysam_fh.get_reference_name(r.reference_id)
            insert_size = self.get_insert_size([_chr,r.reference_start],[sa[0][0],sa[0][1]])
            
            #print r.is_reverse, r.mate_is_reverse
            #if r.is_reverse == r.mate_is_reverse:
                #print "****************"
                #print r
            
            if abs(insert_size) >= 400:
                if r.get_tag('RG') == 'discordant_mates':
                    c.insert(r,sa)
                #print r.get_aligned_pairs()
            else:
                # Figure out whether there was hard or soft clipping on one of the mates
                pass
