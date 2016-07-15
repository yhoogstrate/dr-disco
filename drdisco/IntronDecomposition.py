#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Tries to figure out within a discordant RNA-Seq read alignment:
 - Tries to explain introns and exons by fusion transcripts
"""

#http://www.samformat.info/sam-format-flag

import logging,re
import pysam
from intervaltree_bio import GenomeIntervalTree
from .CigarAlignment import *

from fuma.Fusion import STRAND_FORWARD, STRAND_REVERSE, STRAND_UNDETERMINED 



class Arc:
    def __init__(self,_target):
        if not isinstance(_target, Node):
            raise Exception("_trget must be a Node")
        
        self._target = _target
        self._types = {}
    
    def add_type(self,_type):
        if not self._types.has_key(_type):
            self._types[_type] = 0
        
        self._types[_type] += 1
    
    def get_score(self):
        """Based on this function, the start point is determined
        
        Currently it's implemented as sum(split reads)
        
        Later on spanning reads will be added, softclips will be added, etc.
        """
        
        score = 0
        #print self._types
        
        if self._types.has_key('spanning_paired'):
            score += self._types['spanning_paired']*3
        
        if self._types.has_key('spanning_singleton'):
            score += self._types['spanning_singleton']*2
        
        #if self._types.has_key('spannin'):
        #    score += self._types['spanning_singleton']*2
        
        return score



class Node:
    def __init__(self,position):
        self.position = position
        self.arcs = {}
    
    def insert_arc(self,arc,arc_type):
        skey = str(arc._target)
        
        if not self.arcs.has_key(skey):
            self.arcs[skey] = Arc(arc._target)
        
        self.arcs[skey].add_type(arc_type)
    
    def get_top_arc(self):
        maxscore = -1
        top_arc = None
        
        for arc in self:
            score = arc.get_score()
            if score > maxscore:
                maxscore = score
                top_arc = arc
        
        return (maxscore, top_arc, self, top_arc._target)
    
    def add_arc(self,node2,arc_type,do_vice_versa=True):
        if do_vice_versa:
            node2.add_arc(self,arc_type,False)
        
        arc = Arc(node2)
        self.insert_arc(arc,arc_type)
    
    def __iter__(self):
        for k in sorted(self.arcs.keys()):
            yield self.arcs[k]
    
    def __str__(self):
        return str(self.position)
    
    def str2(self):
        out = ""
        
        for sarc in self.arcs:
            arc = self.arcs[sarc]
            out += "\n\t-> "+str(arc._target)+" "+str(arc._types)
        
        return out+"\n"

def bam_parse_alignment_offset(cigartuple):
    pos = 0
    for chunk in cigartuple:
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

def bam_parse_alignment_end(read):
    pos = read.reference_start
    if not read.is_reverse:
        pos += bam_parse_alignment_offset(read.cigar)
    
    return pos

def bam_parse_alignment_pos_using_cigar(sa_tag):
    """SA tag looks like this:
        ['chr21', 42879876, '85S41M', '3', '-', '1']
    """
    
    pos = sa_tag[1]
    if sa_tag[4] == '+':
        cigartuple = cigar_to_cigartuple(sa_tag[2])
        pos += bam_parse_alignment_offset(cigartuple)
    
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
        self.idxtree = GenomeIntervalTree()
        #self.idx = {}
        self.pysam_fh = pysam_fh
    
    def create_node(self,pos):
        """
        Creates the Node but does not overwrite it
        """
        # see if position exists in idx
        if len(self.idxtree[pos._chr][pos.pos]) == 0:
            self.idxtree[pos._chr].addi(pos.pos,pos.pos+1,{})
        
        # make sure the strand is added
        if not list(self.idxtree[pos._chr][pos.pos])[0][2].has_key(pos.strand):
            list(self.idxtree[pos._chr][pos.pos])[0][2][pos.strand] = Node(pos)
    
    def get_node_reference(self,pos):
        try:
            return list(self.idxtree[pos._chr][pos.pos])[0][2][pos.strand]
        except:
            return None
    
    def insert_entry(self,pos1,pos2,_type):
        """
         - Checks if Node exists at pos1, otherwise creates one
         - Checks if Node exists at pos2, otherwise creates one
         - Checks if Arc exists between them
        """
        
        self.create_node(pos1)
        self.create_node(pos2)
        
        node1 = self.get_node_reference(pos1)
        node2 = self.get_node_reference(pos2)
        
        node1.add_arc(node2,_type)
    
    def insert(self,read,parsed_SA_tag,specific_type = None):
        """Inserts a bi-drectional arc between read and sa-tag in the Chain
        determines the type of arc by @RG tag (done by dr-disco fix-chimeric)"""
        # Type:
        # 2. 'discordant_read'
        # 3. 'split_read'
        # 4. 'silent_mate'
        
        rg = read.get_tag('RG')
        if rg in ["discordant_mates","silent_mate","spanning_singleton"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 bam_parse_alignment_end(read),
                                 not read.is_reverse)
            
            if read.mate_is_reverse:
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1],
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            else:
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg)
        
        elif rg in ["spanning_paired"]:
            ca = CigarAlignment(read.cigar, cigar_to_cigartuple(parsed_SA_tag[2]))
            
            if ca.get_order() == STRAND_FORWARD:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     bam_parse_alignment_end(read),
                                     not read.is_reverse)
                
                #if read.mate_is_reverse:
                if parsed_SA_tag[4] == "+":
                    pos2 = BreakPosition(parsed_SA_tag[0],
                                        parsed_SA_tag[1],
                                        STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
                else:
                    pos2 = BreakPosition(parsed_SA_tag[0],
                                        bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                        STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            else:
                if not read.is_reverse:
                    pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                        read.reference_start,
                                        not read.is_reverse)
                else:
                    pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                        bam_parse_alignment_end(read),
                                        not read.is_reverse)
                    
                
                ## The read is the second in pair, and the SA tag the first...
                if parsed_SA_tag[4] != "+":
                    print read
                    raise Exception("Todo")
                    # Most likely code that should do it:
                    #
                    #pos1 = BreakPosition(parsed_SA_tag[0],
                    #                    parsed_SA_tag[1],
                    #                    STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
                    #
                    #self.insert_entry(pos1,pos2,rg)
                    
                else:
                    pos1 = BreakPosition(parsed_SA_tag[0],
                                        bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                        STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
                
            self.insert_entry(pos1,pos2,rg)
            
        else:
            raise Exception("Fatal Error, RG: "+rg)
    
    def __iter__(self):
        for key in self.idxtree:
            for element in self.idxtree[key]:
                for strand in sorted(element[2].keys()):
                    yield element[2][strand]
    
    def print_chain(self):
        for node in self:
            key = node.position
            print key, node.str2()

    def prune_pos(self, insert_size, node1, node2):
        ## @todo make somehow only 
        #if 
        print insert_size
        print node1
        print node2
    
    def get_start_point(self):
        """
        Returns the chain with the higest number of counts
        """
        maxscore = -1
        node1 = None
        node2 = None
        
        for node in self:
            _score, _arc, _node1, _node2 = node.get_top_arc()
            if _arc != None and _score > maxscore:
                maxscore = _score
                node1  = _node1
                node2  = _node2
        
        if node1 != None:
            return (node1, node2)
        else:
            return None
    
    def prune(self,insert_size):
        """
        Does some 'clever' tricks to merge arcs together and reduce data points
        """
        self.print_chain()
        
        init = self.get_start_point()
        if init != None:
            self.prune_pos(insert_size, init[0], init[1])
            


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
        sa_tags = SA_tag.split(";")
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
        
        #@todo move loop into a function insert_tree() ?
        for r in pysam_fh.fetch(lpos[0],lpos[1]):
            sa = self.parse_SA(r.get_tag('SA'))
            _chr = pysam_fh.get_reference_name(r.reference_id)
            insert_size = self.get_insert_size([_chr,r.reference_start],[sa[0][0],sa[0][1]])
            
            if r.get_tag('RG') == 'discordant_mates':
                if abs(insert_size) >= 400:
                    c.insert(r,sa[0])
            
            elif r.get_tag('RG') == 'silent_mate':
                # usually in a exon?
                if r.is_read1:
                    # The complete mate is the secondary, meaning:
                    # HI:i:1       HI:i:2
                    # [====]---$---[====>-----<========]
                    # The 'break' between the mates is from spanning_paired.2 <-> read
                    # -- requires spanning_paired.1 and spanning_paired.2 to be set in the correct order --
                    
                    broken_mate = sa[1]
                    #broken_mate[1] += bam_parse_alignment_offset_using_cigar(broken_mate)

                else:# is_read2
                    # The complete mate is the primary, meaning:
                    #                HI:i:1       HI:i:2
                    # [========>-----<====]---$---[====]
                    # The 'break' between the mates is from read <-> spanning_paired.1
                    # -- requires spanning_paired.1 and spanning_paired.2 to be set in the correct order --
                    
                    broken_mate = sa[0]
                
                c.insert(r,broken_mate)
            
            elif r.get_tag('RG') in ['spanning_paired', 'spanning_singleton']:
                read = r
                c.insert(r,sa[0])

            else:
                raise Exception("Unknown type read: '"+str(r.get_tag('RG'))+"'. Was the alignment fixed with a more up to date version of Dr.Disco?")
            
            # Find introns etc:
            for internal_arc in self.find_cigar_arcs(r):
                pos1 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                     internal_arc[0],
                                     not r.is_reverse)
                pos2 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                     internal_arc[1],
                                     not r.is_reverse)
                
                c.insert_entry(pos1,pos2,internal_arc[2])
        
        # Merge arcs somehow, label nodes
        c.prune(400)

    def find_cigar_arcs(self,read):
        """Tries to find ARCs introduced by:
         - Hard clipping
         - Soft clipping
         - Splicing
         - Deletion

            M	BAM_CMATCH	0
            I	BAM_CINS	1
            D	BAM_CDEL	2
            N	BAM_CREF_SKIP	3
            S	BAM_CSOFT_CLIP	4
            H	BAM_CHARD_CLIP	5
            P	BAM_CPAD	6
            =	BAM_CEQUAL	7
            X	BAM_CDIFF	8"""
        
        tt = {
            2:'cigar_deletion',
            3:'cigar_splice_junction',
            4:'cigar_soft_clip',
            5:'cigar_hard_clip'
        }
        
        offset = read.reference_start
        for chunk in read.cigar:
                          # D N S H
            if chunk[0] in tt.keys():
                # Small softclips occur all the time - introns and 
                # deletions of that size shouldn't add weight anyway
                if chunk[1] > 3:
                    yield (offset , (offset + chunk[1]) , tt[chunk[0]])
            
            offset += chunk[1]

