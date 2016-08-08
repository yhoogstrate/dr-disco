#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Tries to figure out within a discordant RNA-Seq read alignment:
 - Tries to explain introns and exons by fusion transcripts
"""

#http://www.samformat.info/sam-format-flag

import logging,re,math,copy,sys
import pysam
from intervaltree_bio import GenomeIntervalTree, Interval
from .CigarAlignment import *
from .CircosController import *

from fuma.Fusion import STRAND_FORWARD, STRAND_REVERSE, STRAND_UNDETERMINED 


# parameters
MIN_DISCO_INS_SIZE = 400
PRUNE_INS_SIZE = 450
SPLICE_JUNC_ACC_ERR = 3 # acceptable splice junction error


# translation tables:
strand_tt = {STRAND_FORWARD:'+',STRAND_REVERSE:'-',STRAND_UNDETERMINED:'?'}


def entropy(frequency_table):
    n = sum(frequency_table.values())
    prob = [float(x)/n for x in frequency_table.values()]
    
    entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob if p > 0 ])
    if entropy <= 0:
        return 0.0
    else:
        return entropy / (math.log(n) / math.log(2.0))


class Arc:
    """Connection between two genomic locations
     - Also contains different types of evidence
    """
    
    scoring_table={
               'cigar_splice_junction': 0,
               'discordant_mates': 2,
               'silent_mate': 0,
               'spanning_paired_1': 3,
               'spanning_paired_1_r': 3,
               'spanning_paired_1_s': 3,
               'spanning_paired_1_t': 3,
               'spanning_paired_2': 3,
               'spanning_paired_2_r': 3,
               'spanning_paired_2_s': 3,
               'spanning_paired_2_t': 3,
               'spanning_singleton_1': 2, 
               'spanning_singleton_2': 2,
               'spanning_singleton_1_r': 1,
               'spanning_singleton_2_r': 1
               }
    
    def __init__(self,_origin,_target):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        if not isinstance(_target, Node) or not isinstance(_origin, Node):
            raise Exception("_origin and _target must be a Node")
        
        #self._strand = STRAND_UNDETERMINED
        
        self._origin = _origin
        self._target = _target
        self._types = {}
        
        # Used for determining entropy
        self.unique_alignments_idx = {}
    
    """
    def remove_split_reads(self):
        for key in ['spanning_paired_1',
               'spanning_paired_2',
               'spanning_singleton_1',
               'spanning_singleton_2',
               'spanning_singleton_1_r',
               'spanning_singleton_2_r']:
                """
    
    def get_entropy(self):
        """Entropy is an important metric/ propery of an arc
        The assumption is that all reads should be 'different', i.e.
        have a different start position, different end position,
        different soft/hardclips etc.
        
        It is more probable to find the following alignment in a true
        fusion gene:
        
            <==========]                              (split reads)
              <========]
              <========]
               <=======]
                <======]
                               <==========]
                               <========]
                               <========]
                               <=======]
                               <======]
       <===========]------------------<===========]    (discordant)
     <===========]------------------<===========]
  <===========]------------------<===========]
        
        than:

              <========]
              <========]
              <========]
              <========]
              <========]
                               <========]
                               <========]
                               <========]
                               <========]
                               <========]
        <===========]------------------<===========]    (discordant)
        <===========]------------------<===========]
       
       
        The easiest way to calculate this is based upon the position plus
        the cigar strings. However, if some kind of merge strategy
        will be added, we may find 'M126' for different positions while
        they have a different meaning. Hence, the combination of the
        position + the CIGAR string would be the unique key we would like
        to use for the fequency table.
        
        Of course there is some level of saturation; if you have 1.000.000
        reads there will be reads that are identical because there is simply
        a limited spanning region (determined by insert size + read length)
        but calculating the true possiblities has too many degrees of freedom.
        """
        
        return entropy(self.unique_alignments_idx)
    
    def get_complement(self):
        """
        complement = self._target.arcs[str(self._origin.position)]
        
        if complement == self:
            raise Exception("err")
        
        return complement
        """

        return self._target.arcs[str(self._origin.position)]
    
    def merge_arc(self,arc):
        """Merges discordant_mates arcs
        """
        
        for alignment_key in arc.unique_alignments_idx:
            self.add_alignment_key(alignment_key)
        
        for _type in arc._types:
            if _type in [
                    "discordant_mates",
                    "spanning_paired_1",
                    "spanning_paired_1_r",
                    "spanning_paired_1_s",
                    "spanning_paired_1_t",
                    "spanning_paired_2",
                    "spanning_paired_2_r",
                    "spanning_paired_2_s",
                    "spanning_paired_2_t",
                    "spanning_singleton_1",
                    "spanning_singleton_2",
                    "spanning_singleton_1_r",
                    "spanning_singleton_2_r",
                    "cigar_splice_junction",
                    ]:
                self.add_type(_type)
                return True
            elif _type not in ['silent_mate']:
                raise Exception("Not sure what to do here with type: %s", _type)
        return False
    
    def add_type(self,_type):
        if _type in ["cigar_soft_clip", 'cigar_hard_clip']:
            raise Exception("Clips shouldn't be added as arcs, but as properties of Nodes")
        
        if not self._types.has_key(_type):
            self._types[_type] = 0
        
        self._types[_type] += 1
    
    def add_alignment_key(self,alignment_key):
        if not self.unique_alignments_idx.has_key(alignment_key):
            self.unique_alignments_idx[alignment_key] = 0
        
        self.unique_alignments_idx[alignment_key] += 1
    
    def get_count(self, _type):
        
        if not self._types.has_key(_type):
            return 0
        else:
            return self._types[_type]
        
    def get_score(self,_type):  
        if not self.scoring_table.has_key(_type):
            raise Exception("Not implemented _type: %s", _type)
        else:
            return self.get_count(_type)*self.scoring_table[_type]
    
    def get_scores(self):
        """Based on this function, the start point is determined
        
        Currently it's implemented as sum(split reads)
        
        Later on spanning reads will be added, softclips will be added, etc.
        """
        
        score = 0

        for _type in self.scoring_table.keys():
            score += self.get_score(_type)
        
        return score
    
    def get_splice_score(self):
        #@todo use soft/hardclips or the nodes
        return (self.get_count('cigar_splice_junction'),self.get_clips())
    
    def get_clips(self):
        return self._origin.clips + self._target.clips
    
    def target_in_range(self,_range):
        if _range[0] > _range[1]:
            _max = _range[0]
            _min = _range[1]
        else:
            _min = _range[0]
            _max = _range[1]
            
        return (
                #(self._origin.position.pos >= _min) and (self._origin.position.pos <= _max)
                 #or 
                         (self._target.position.pos >= _min) and (self._target.position.pos <= _max)
                )
    
    def __str__(self):
        typestring = []
        for _t in sorted(self._types.keys()):
            typestring.append(str(_t)+":"+str(self._types[_t]))
        
        out = str(self._origin.position) + "->" + str(self._target.position) + ":("+','.join(typestring)+")"
        
        #spacer = " "*len(str(self._origin.position))
        
        #for _k in self._origin.splice_arcs.keys():
        #    out += "\n"+spacer+"=>"+str(_k.position)+":score="+str(self._origin.splice_arcs[_k][1].get_splice_score())
        
        return out


class Node:
    """A Simple node, just a chromomal location( and strand)"""
    
    def __init__(self,position):
        self.position = position
        self.clips = 0
        self.arcs = {}
        self.splice_arcs = {}
    
    def rfind_connected_sjuncs(self,left_nodes):
        """Recursively finds all nodes that are connected by splice junctions"""
        
        results = []
        for node in self.splice_arcs.keys():
            if node not in left_nodes:
                results.append(node)
        
        results_all = [x for x in left_nodes]
        for x in results:
            if x not in results_all:
                results_all.append(x)
        
        for node in results:
            for child in node.rfind_connected_sjuncs(results_all):
                if child not in results_all:
                    results_all.append(child)
        
        return results_all
    
    def add_clip(self):
        self.clips += 1
    
    def insert_arc(self,arc,arc_type,alignment_key):
        skey = str(arc._target.position)
        
        if not self.arcs.has_key(skey):
            self.set_arc(arc)
        
        self.arcs[skey].add_type(arc_type)
        
        if alignment_key != None:# Splice junctions should be skipped for entropy
            self.arcs[skey].add_alignment_key(alignment_key)
    
    def set_arc(self, arc):
        self.arcs[str(arc._target.position)] = arc
    
    def get_top_arc(self):
        maxscore = -1
        top_arc = None
        
        for arc in self:
            score = arc.get_scores()
            if score > maxscore:
                maxscore = score
                top_arc = arc
        
        if top_arc != None:
            return (maxscore, top_arc, self, top_arc._target)
        else:
            return (None, None, None, None)
    
    def new_arc(self,node2,arc_type,alignment_key,do_vice_versa):
        if do_vice_versa:
            node2.new_arc(self,arc_type,alignment_key,False)
        
        arc = Arc(self,node2)
        self.insert_arc(arc,arc_type,alignment_key)
    
    def remove_arc(self,arc,idx):
        if idx == "by-target":
            skey = str(arc._target.position)
        elif idx == "by-origin":
            skey = str(arc._origin.position)
        else:
            raise Exception("Invalid usage of function")
        
        if not self.arcs.has_key(skey):
            raise Exception("Unknown key: %s", skey)
        
        self.arcs[skey] = None
        del(self.arcs[skey])
    
    def __iter__(self):
        for k in sorted(self.arcs.keys()):
            yield self.arcs[k]
    
    def __str__(self):
        out  = str(self.position)
        
        a = 0
        for sarc in self.arcs:
            arc = self.arcs[sarc]
            filtered_arcs = {x:arc._types[x] for x in sorted(arc._types.keys()) if x not in ['cigar_soft_clip','cigar_hard_clip']}
            len_arcs = len(filtered_arcs)
            a += len_arcs
            if len_arcs > 0:
                out += "\n\t-> "+str(arc._target.position)+" "+str(filtered_arcs)

        
        if a > 0:
            out += "\n\t-> soft/hard clips: "+str(self.clips)
            
            return out+"\n"
        else:
            return ""

    def is_connected_to(self, nearby_nodes, complementary_nodes):
        """The settings in this function may be very important to tune
        in order to optimize the algorithm
        
        @todo's in future, add entropy of nearby nodes using e.g. the
        cigar strings"""
        
        rmse = self.position.get_rmse(complementary_nodes)
        
        scores = []
        for nearby_node in nearby_nodes:
            p = nearby_node.arcs[str(self.position)]
            scores.append(p.get_scores())
        
        k = len(scores)
        if k < 2:
            raise Exception("nearby_nodes contains less than 2 elements with shared Arcs")
        
        product = reduce(lambda x,y: x*y, scores)
        
        # based on rmse, product and k, determine a True or False
        if product > (8*k):
            # Average gene size = 10-15kb
            # Asymptotic func:
            # rmse <= 125000 - 110123/2^( product / 620 )
            # rmse <= 125000 - 120000/2^( product / 1200)
            # if rmse < max_rmse: valid data point
            
            max_rmse = 125000.0 - (120000/pow(2, float(product) / 1200.0))
            return (rmse <= max_rmse)
        
        return False

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
            X	BAM_CDIFF	8
            """
        
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
    
    def get_dist(self, other_bp, strand_specific):
        if not isinstance(other_bp, BreakPosition):
            raise Exception("Wrong data type used")
        
        if (not strand_specific or self.strand == other_bp.strand) and self._chr == other_bp._chr:
            return other_bp.pos - self.pos
        else:
            # chr1 has 249 000 000 bp - by replacing all digits by 9 is
            # must be larger than any Hs chr reflecting natural distances
            # in a way that interchromosomal breaks are 'larger' than
            # intrachromosomal ones
            return 999999999
    
    def get_rmse(self, pos_vec):
        err = 0
        for node in pos_vec:
            pos = node.position
            err += pow(self.get_dist(pos, True) , 2)
        return math.sqrt(err)


class Chain:
    def __init__(self,pysam_fh):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        self.idxtree = GenomeIntervalTree()
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
        if not isinstance(pos, BreakPosition):
            raise Exception("Wrong data type used")
        
        try:
            return list(self.idxtree[pos._chr][pos.pos])[0][2][pos.strand]
        except:
            return None
    
    def insert_entry(self,pos1,pos2,_type,cigarstrs,do_vice_versa):
        """
         - Checks if Node exists at pos1, otherwise creates one
         - Checks if Node exists at pos2, otherwise creates one
         - Checks if Arc exists between them
         
         cigarstrs must be something like ("126M","126M") or ("25S50M2000N","25M50S")
        """
        
        self.create_node(pos1)
        self.create_node(pos2)
        
        node1 = self.get_node_reference(pos1)
        node2 = self.get_node_reference(pos2)
        
        if cigarstrs != None:
            # Hexadec saves more mem
            short_pos1 = "%0.2X" % pos1.pos#str(pos1.pos)
            short_pos2 = "%0.2X" % pos2.pos#str(pos2.pos)
        
            alignment_key  = short_pos1+strand_tt[pos1.strand]+cigarstrs[0]+"|"
            alignment_key += short_pos2+strand_tt[pos2.strand]+cigarstrs[1]
            
            node1.new_arc(node2,_type,alignment_key,do_vice_versa)
        else:
            node1.new_arc(node2,_type,None,do_vice_versa)
    
    def insert(self,read,parsed_SA_tag,specific_type = None):
        """Inserts a bi-drectional arc between read and sa-tag in the Chain
        determines the type of arc by @RG tag (done by dr-disco fix-chimeric)
        
        
        Strands of reads... magical...
        -------------------------------
        spanning_paired_1	71M41S		+:

             [==========>
                 [======>
                    [===>

        spanning_paired_2	71S41M		+:
                               [=======>
                               [=====>
                               [===>

        (silent_mate is 2nd in pair)



        spanning_paired_1	95M31S		-:
            <==========]
              <========]
              <========]
               <=======]
                <======]

        spanning_paired_2	95S31M		-:
                               <==========]
                               <========]
                               <========]
                               <=======]
                               <======]

        spanning_pair is second in pair


        This means that for:
         - spanning_paired_1 the breakpoint is at: start + offset
         - spanning_paired_2 the breakpoint is at: start
        regardless of the strand (+/-)"""
        # Type:
        # 2. 'discordant_read'
        # 3. 'split_read'
        # 4. 'silent_mate'
        
        rg = read.get_tag('RG')
        if rg in ["discordant_mates"]:
            """
            Normally reads are:
            [==R1==> ... <==R2==]
            
            For transcription in negative strand, I expect:
            [==R2==> ... <==R1==]
            -- confirm --
            -- if this is true, strand + fip/sip are determinant --
            
            
            If first in pair:
            [========> ... <---] junction should be position + offset
            
            
            If second in pair:
              if strand is - (normal):
                [--> ... <========] junction should be start position
              if strand is + (rev):
                [--> ... [========> 
            


            spanning_singleton_2_r
            spanning_paired_2_t (left-junc):
            =======>
             ======>
               ====>


            spanning_paired_2_s:
            ---
            spanning_singleton_1_r (right-junc):
            spanning_paired_1_t    (right-junc):
                      =======>
                      ======>
                      ====>
            
            spanning_paired_1_s:
            ---
            spanning_singleton_1:
            spanning_paired_1:
            <=======
             <======
               <====

            spanning_singleton_2 (right junc):
            spanning_paired_2 (right junc):
                      <=======
                      <======
                      <====

        """

            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 bam_parse_alignment_end(read),
                                 STRAND_FORWARD if read.is_reverse else STRAND_REVERSE)
            
            if read.is_read1:
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1],
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            else:
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        elif rg in ["spanning_singleton_1", "spanning_paired_1"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)

        elif rg in ["spanning_singleton_1_r", "spanning_paired_1_t"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 not read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        elif rg in ["spanning_paired_1_r"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 not read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        elif rg in ["spanning_singleton_2", "spanning_paired_2"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        elif rg in ["spanning_singleton_2_r", "spanning_paired_2_t"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        elif rg in ["spanning_paired_2_r"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        elif rg not in ["silent_mate"]:
            raise Exception("Fatal Error, RG: "+rg)
        
        return (pos1, pos2)
    
    def reinsert_arcs(self, arcs):
        """Only works for Arcs of which the _origin and _target Node
        still exists
        """
        for arc_t in arcs:
            arc   = arc_t[0]
            arc_c = arc_t[1]
            
            node1 = arc._origin
            node2 = arc._target
            
            node1.set_arc(arc)
            node2.set_arc(arc_c)
    
    def remove_arc(self, arc):
        node1 = arc._origin
        node2 = arc._target
        
        node1.remove_arc(arc,"by-target")
        node2.remove_arc(arc,"by-origin")

        # Not necessary
        #del(arc)
        
        if node1.get_top_arc()[0] == None:
            self.remove_node(node1)

        if node2.get_top_arc()[0] == None:
            self.remove_node(node2)
    
    def remove_node(self,node):
        pos = node.position
        element = self.idxtree[pos._chr][pos.pos]
        
        root_element = list(element)[0]
        new_element = Interval(root_element[0],root_element[1],{x:root_element[2][x] for x in root_element[2].keys() if x != pos.strand})
        
        # First remove the root element
        self.idxtree[pos._chr].remove(root_element)
        
        # Delete node
        del(root_element[2][pos.strand])
        
        # Weird error
        #del(root_element)
        
        ## Re-insert if necessary
        if len(new_element[2]) > 0:
            self.idxtree[pos._chr].add(new_element)
            
    def __iter__(self):
        for key in self.idxtree:
            for element in sorted(self.idxtree[key]):
                for strand in sorted(element[2].keys()):
                    yield element[2][strand]
    
    def get_range(self,pos,insert_size,inverse_strand):
        if inverse_strand:
            if pos.strand == STRAND_FORWARD:
                return (pos.pos,(pos.pos-insert_size)-1,-1)
            else:
                return (pos.pos,(pos.pos+insert_size)+1,1)
        else:
            if pos.strand == STRAND_REVERSE:
                return (pos.pos,(pos.pos-insert_size)-1,-1)
            else:
                return (pos.pos,(pos.pos+insert_size)+1,1)
    

    def pos_to_range(self,pos,_range,exclude_itself):
        for _pos in range(_range[0], _range[1], _range[2]):
            if exclude_itself and _pos == pos.pos:
                continue
            else:
                yield BreakPosition(pos._chr, _pos, pos.strand)
        
    
    def search_arcs_between(self,pos1, pos2, insert_size):
        """Searches for reads inbetween two regions (e.g. break + ins. size)
        
        @todo: correct for splice junction width
        """
        
        range1 = self.get_range(pos1,insert_size,False)
        range2 = self.get_range(pos2,insert_size,False)
        for pos_i in self.pos_to_range(pos1,range1,True):
            node_i = self.get_node_reference(pos_i)
            if node_i != None:
                for arc in node_i.arcs.values():
                    """Although I still don't understand why this happens
                    I found about ~25% of the split reads to be in reverse
                    strand:
                    
                    <========]   $ ... $   <====]
                      <======]   $ ... $   <======]
                        <====]   $ ... $   <========]
                    <========]   $ ... $   <====]
                      <======]   $ ... $   <======]
                        <====]   $ ... $   <========]
                        
                       [=====>   $ ... $   [========>
                      <======]   $ ... $   <========]
                    
                    but this almost only seem to happen for 'spanning_singleton_[1/2]' type reads.
                    The most likely reason is that for spanning singletons the direction
                    depends on whether it was the first or second mate.
                    """
                    if arc.target_in_range(range2):
                        yield (arc)

    
        #for pos_i in self.pos_to_range(pos2,range2,True):
            #node_i = self.get_node_reference(pos_i)
            #if node_i != None:
                #for arc in node_i.arcs.values():
                    #if not arc.target_in_range(range1):
                        #yield('pos2',arc)
                     
                    #if arc.target_in_range(range1):
                    #    yield ('both',arc)
                    #else:
                    #    yield('pos2',arc)

    
    def arcs_ratio_between(self,pos1, pos2, insert_size):
        """
        Searches for reads inbetween two regions (e.g. break + ins. size)
        and returns the ratio of uniq reads that fall within this range
        
        included = [
            'discordant_mates',
            'spanning_paired_1',
            'spanning_paired_2',
            'spanning_singleton_1',
            'spanning_singleton_2',
            'spanning_singleton_1_r',
            'spanning_singleton_2_r'
            ]
        excluded = ['cigar_soft_clip', 'silent_mate', 'cigar_splice_junction']
        
        total_arcs = 0
        arcs_inbetween = 0

        range1 = self.get_range(pos1,insert_size,False)
        range2 = self.get_range(pos2,insert_size,False)
        for pos_i in self.pos_to_range(pos1,range1,False):
            node_i = self.get_node_reference(pos_i)
            if node_i != None:
                for arc in node_i.arcs.values():
                    for _type in arc._types.keys():
                        if _type in included:
                            if arc.target_in_range(range2):
                                arcs_inbetween += arc._types[_type]
                                total_arcs += arc._types[_type]
                            else:
                                # The in-range ones are counted twice
                                total_arcs += 2*arc._types[_type]
                        
                        elif _type in excluded:# weird types; cigar types
                            continue
                        else:
                            raise Exception("To be implemented: %s", _type)
                    
        
        range2 = self.get_range(pos2,insert_size,False)
        range1 = self.get_range(pos1,insert_size,False)
        for pos_i in self.pos_to_range(pos2,range2,False):
            node_i = self.get_node_reference(pos_i)
            if node_i != None:
                for arc in node_i.arcs.values():
                    for _type in arc._types.keys():
                        if _type in included:
                            if arc.target_in_range(range1):
                                arcs_inbetween += arc._types[_type]
                                total_arcs += arc._types[_type]
                            else:
                                # The in-range ones are counted twice
                                total_arcs += 2*arc._types[_type]
                        
                        elif _type in excluded:# weird types; cigar types
                            continue
                        else:
                            raise Exception("To be implemented: %s", _type)

        return (1.0 * arcs_inbetween / total_arcs)
        """

    def print_chain(self):
        print "**************************************************************"
        for node in self:
            _str = str(node)
            if len(_str.strip()) > 0:
                print _str
        print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    def get_start_point(self):
        """Returns the chain with the higest number of counts
        """
        maxscore = 0
        arc = None
        
        for node in self:
            _score, _arc, _node1, _node2 = node.get_top_arc()
            if _arc != None and _score > maxscore:
                maxscore = _score
                arc = _arc
        
        return arc
    
    def get_start_splice_junc(self):
        """Returns the chain with the higest number of counts
        """
        maxscore = (0,0)
        arc = None
        
        for node in self:
            for _arc in node:
                score = (_arc.get_count('cigar_splice_junction'),_arc.get_count('cigar_soft_clip'))
                if score[0] > maxscore[0] or (score[0] == maxscore[0] and score[1] > maxscore[1]):
                    arc = _arc
                    maxscore = score
        
        return arc
    
    def prune(self,insert_size):
        """Does some 'clever' tricks to merge arcs together and reduce data points
        """
        self.logger.info("Find and merge other arcs in close proximity (insert size)")
        
        candidates = []
        candidate = self.get_start_point()
        
        self.print_chain()
        
        while candidate != None:
            self.prune_arc(candidate, insert_size)
            candidates.append((candidate, candidate.get_complement()))
            
            # do not remove if splice junc exists?
            self.remove_arc(candidate)
            
            candidate = self.get_start_point()
        
        self.print_chain()
        
        for candidate in candidates:
            print candidate[0]
        
        return candidates
    
    #@todo ADD RECURSION DEPTH LIMIT
    def prune_arc(self, arc, insert_size):
        ## These types of arcs fall within the same space / inst. size
        ## they make the arcs heavier
        
        node1 = arc._origin
        node2 = arc._target
        
        arc_complement = arc.get_complement()
        
        for arc_m in self.search_arcs_between(node1.position, node2.position, insert_size):
            arc_mc = arc_m.get_complement()
            
            arc.merge_arc(arc_m)
            arc_complement.merge_arc(arc_mc)
            
            self.remove_arc(arc_m)
            # complement is automatically removed after removing the fwd
            #self.remove_arc(arc_mc)
        
        return None
    
    def merge_splice_juncs(self,uncertainty):
        self.logger.info("Merging splice juncs")
        
        init = self.get_start_splice_junc()
        if init != None:
            init_c = init.get_complement()
            
            for node in self:
                for junc in node:
                    if junc.get_count('cigar_splice_junction') > 0 and junc not in [init, init_c]:
                        d_origin = abs(junc._origin.position.pos - init._origin.position.pos)
                        d_target = abs(junc._target.position.pos - init._target.position.pos)
                        if d_origin <= uncertainty and d_target <= uncertainty:
                            init.merge_arc(junc)
                            init_c.merge_arc(junc.get_complement())
                            
                            self.remove_arc(junc)
                            
                            #self.print_chain()
            
            init = self.get_start_splice_junc()
    
    def rejoin_splice_juncs(self, thicker_arcs, insert_size):
        """thicker arcs go across the break point:
        
 splice juncs:
               -------  -------                                            -----
             ||       ||       ||        |          $ ... $       |      ||     ||

thick arcs:
                                ------------------------------------------
                       ---------------------------------------------------
                       ----------------------------------------------------------
              -------------------------------------------------------------------
        
        the goal is to add the splice juncs between the nodes
        """
        ## 01 collect all left and right nodes
        
        #added_splice_junctions = set()
        
        left_nodes = set()
        right_nodes = set()
        
        for arc in thicker_arcs:
            left_nodes.add(arc[0]._origin)
            right_nodes.add(arc[0]._target)
        
        ## 02 look for all left nodes if there is any set (i < j) where
        ## i and j span a splice junction
        i = -1
        for node1 in left_nodes:
            i += 1
            j = -1
            
            for node2 in left_nodes:
                j += 1
                
                if j > i:# Avoid unnecessary comparisons
                    if node1.position.strand == node2.position.strand:
                        left_junc = (999999999, None)
                        
                        for node in self:
                            for splice_junc in node:
                                if splice_junc.get_count('cigar_splice_junction') > 0:#@todo and dist splice junction > ?insert_size?
                                    dist_origin1 = abs(splice_junc._origin.position.get_dist(node1.position, False))
                                    dist_origin2 = abs(splice_junc._target.position.get_dist(node2.position, False))
                                    sq_dist_origin = pow(dist_origin1, 2) +pow(dist_origin2, 2)
                                    
                                    if dist_origin1 < insert_size and dist_origin2 < insert_size and sq_dist_origin < left_junc[0]:
                                        left_junc = (sq_dist_origin, splice_junc)
                        
                        if left_junc[1] != None:
                            node1.splice_arcs[node2] = left_junc
                            node2.splice_arcs[node1] = left_junc
                            
                            #added_splice_junctions.add(left_junc[1])

        i = -1
        for node1 in right_nodes:
            i += 1
            j = -1
            
            for node2 in right_nodes:
                j += 1
                
                if j > i:# Avoid unnecessary comparisons
                    if node1.position.strand == node2.position.strand:
                        right_junc = (999999999, None)
                        
                        for node in self:
                            for splice_junc in node:
                                if splice_junc.get_count('cigar_splice_junction') > 0:#@todo and dist splice junction > ?insert_size?
                                    dist_target1 = abs(splice_junc._origin.position.get_dist(node1.position, False))
                                    dist_target2 = abs(splice_junc._target.position.get_dist(node2.position, False))
                                    sq_dist_target = pow(dist_target1, 2) +pow(dist_target2, 2)

                                    if dist_target1 < insert_size and dist_target2 < insert_size and sq_dist_target < right_junc[0]:
                                        right_junc = (sq_dist_target, splice_junc)
                        
                        if right_junc[1] != None:
                            node1.splice_arcs[node2] = right_junc
                            node2.splice_arcs[node1] = right_junc
                            
                            #added_splice_junctions.add(right_junc[1])
        
        return thicker_arcs
    
    def extract_subnetworks(self,thicker_arcs):
        """Make sure this does not suffer from endless recursion
        """
        q = 0
        subnetworks = []
        while len(thicker_arcs) > 0:
            start_point = thicker_arcs[0][0]
            
            left_nodes = [start_point._origin]
            right_nodes = [start_point._target]
            
            ## The original nodes have been emptied, so the most important
            ## arc's are now separated.
            left_nodes = start_point._origin.rfind_connected_sjuncs(left_nodes)
            right_nodes = start_point._target.rfind_connected_sjuncs(right_nodes)
            
            subarcs = []
            
            # All arcs between any of the left and right nodes are valid arcs and have to be extracted
            # Find all direct arcs joined by splice junctions
            for left_node in left_nodes:
                for right_node in right_nodes:
                    if left_node.arcs.has_key(str(right_node.position)):
                        subarcs.append( (left_node.arcs[str(right_node.position)], right_node.arcs[str(left_node.position)]) )
            
            
            ## Find all with one indirect step - these might be alternative junctions / exons
            """
            Here we want to add new nodes to `left_nodes` or `right_nodes`
            using the `guilt-by-association` principle. Sometimes nodes are
            having arcs to the same nodes of the already existing network,
            but lack splice junction(s) to those existing network. They're
            still connected to the network they might be exons that are not
            taken into account by the aligner (classical example: exon-0 in
            TMPRSS2). We have to be careful on the other hand, as multimap
            locations may influence this.
            
            Therefore we ideally need a function that adds nodes based on:
             - Genomic distance to all entries in either `left_nodes` or `right_nodes`
               * Currently implemented as RMSE on Node::get_dist()
             - The score of all arcs going to the node
               * It also needs to correct for imbalance; 4 and 36 reads is less
                 likely than 20 and 20. Solution used multiplication.
             - Number of arcs relative to the number of nodes.
               *In case of 3 nodes, it is more likely to have 3 arcs (3/3) instead of 2
               (2/3). This weight is not yet implemented (01/08/16).


            Detection of possible candidates
            --------------------------------
            
            nodes:
            A     B         $ ... $          Y    Z
    splice:  -----

    left_nodes:  A, B
    right_nodes: Y, Z

    arcs:
            A-Z (9)
            B-Z (100)
            
            A-Y (2)
            B-Y (20)
            
            For all each left node i, compared with each other left node j,
            look for nodes they have in common:
            Y, Z (exclude Z, because it was already in `right_nodes`
            
            Then we know for sure that both A and B are also connected to
            Y, while there is no relation found between Y-Z. However,
            if the distance between Y-Z is reasonable and there are many
            reads, it is pretty likely that Y is part of the fusion event
            as well.
            
            Reverse
            -------
            Then do this in reverse(d) order, from `right_nodes` to
            `left_nodes`.
            """
            i = -1
            for left_node_i in left_nodes:
                j = -1
                i += 1
                
                for left_node_j in left_nodes:
                    j += 1
                    
                    if i < j:
                        ## Find similar destinations
                        mutual_targets = list(set(left_node_i.arcs.keys()).intersection(left_node_j.arcs.keys()))
                        mutual_targets = [left_node_i.arcs[mt]._target for mt in mutual_targets]
                        
                        for mt in mutual_targets:
                            if mt.is_connected_to((left_node_i, left_node_j), right_nodes):
                                # Add node
                                right_nodes.append(mt)
                                
                                for arc in mt.arcs.keys():
                                    for l in left_nodes:
                                        if str(l.position) == arc:
                                            arc = mt.arcs[arc]
                                            arc_c = arc.get_complement()
                                            
                                            # Make sure order is correct:
                                            subarcs.append((arc_c,arc))
            
            ## Do inverse:
            
            #tmp = left_nodes
            #left_nodes = right_nodes
            #right_nodes = tmp

            ### @todo Redo code
            
            
            # remove all the links to the arcs in each of the nodes
            for node in left_nodes:
                for arc_u in subarcs:
                    for arc in arc_u:
                        key = str(arc._target.position)
                        if node.arcs.has_key(key):
                            del(node.arcs[key])

            # pop subarcs from thicker arcs and redo until thicker arcs is empty
            popme = set()
            for arc in subarcs:
                for arc2 in thicker_arcs:
                    if arc[0] == arc2[0] or arc[1] == arc2[0]:
                        popme.add(arc2)
            
            for pop in popme:
                thicker_arcs.remove(pop)
            
            if q == 2:
                break
            
            subnetworks.append(subarcs)
        
        #for sn in subnetworks:
        #    print "sn:"
        #    for s in sn:
        #        print "  ",s[0]
        
        return subnetworks


class Subnet(Chain):
    def __init__(self,_id,arcs):
        self._id = _id
        self.arcs = arcs
        self.total_clips = 0
        self.total_score = 0
    
    def __str__(self):
        """Make tabular output
        
        print "chr-A\t"
        print "pos-A\t"
        print "direction-A\t"
        
        print "chr-B\t"
        print "pos-B\t"
        print "direction-B\t"
        
        print "score\t"
        print "soft+hardclips\t"
        print "n-split-reads\t"
        print "n-discordant-reads"
        
        print "n-arcs\t"
        print "n-nodes-A\t"
        print "n-nodes-B\t"
        
        fh.write("entropy-bp-arc\t")
        fh.write("entropy-all-arcs\t")
        """
        out = ""
        node_a = self.arcs[0][0]._origin
        node_b = self.arcs[0][0]._target
        
        out += str(node_a.position._chr)+"\t"
        out += str(node_a.position.pos)+"\t"
        out += strand_tt[node_a.position.strand]+"\t"
        
        out += str(node_b.position._chr)+"\t"
        out += str(node_b.position.pos)+"\t"
        out += strand_tt[node_b.position.strand]+"\t"
        
        out += str(self.total_score)+"\t"
        out += str(self.total_clips)+"\t"
        
        out += str(self.get_n_split_reads())+"\t"
        out += str(self.get_n_discordant_reads())+"\t"
        
        out += str(len(self.arcs))+"\t"
        nodes_a, nodes_b = self.get_n_nodes()
        out += str(nodes_a)+"\t"
        out += str(nodes_b)+"\t"
        
        out += str(self.arcs[0][0].get_entropy())+"\t"
        out += str(self.get_overall_entropy())+"\t"
        
        out += "&".join([str(arc[0]) for arc in self.arcs])
        
        return out+"\n"
    
    def get_overall_entropy(self):
        frequency_table = {}
        for arc in self.arcs:
            for key in arc[0].unique_alignments_idx:
                if not frequency_table.has_key(key):
                    frequency_table[key] = 0
                frequency_table[key] += arc[0].unique_alignments_idx[key]
        
        return entropy(frequency_table)
    
    def get_n_splice_junctions(self):
        pass
    
    def get_n_nodes(self):
        nodes_a = set()
        nodes_b = set()
        
        for arc in self.arcs:
            nodes_a.add(arc[0]._origin)
            nodes_b.add(arc[0]._target)
        
        return len(nodes_a), len(nodes_b)
    
    def get_n_split_reads(self):
        n = 0
        
        for _type in ["spanning_paired_1",
                      "spanning_paired_1_r",
                      "spanning_paired_1_s",
                      "spanning_paired_1_t",
                      "spanning_paired_2",
                      "spanning_paired_2_r",
                      "spanning_paired_2_s",
                      "spanning_paired_2_t",
                      "spanning_singleton_1",
                      "spanning_singleton_1_r",
                      "spanning_singleton_2",
                      "spanning_singleton_2_r"]:
            for arc in self.arcs:
                n += arc[0].get_count(_type)
        
        return n
    
    def get_n_discordant_reads(self):
        return sum([arc[0].get_count("discordant_mates") for arc in self.arcs])
    

class IntronDecomposition:
    def __init__(self,break_point):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        self.break_point = break_point
        self.chain = None
    
    def decompose(self,alignment_file):
        pysam_fh = self.test_disco_alignment(alignment_file)
        
        lpos = self.break_point.get_left_position(True)
        rpos = self.break_point.get_right_position(True)
        
        self.chain = Chain(pysam_fh)
        
        # Solve it somehow like this:
        #self.insert_tree(pysam_fh, lpos())
        tmprss2 = ['chr21',42834478,42882085]
        erg = ['chr21',39737183,40035618]
        self.insert_chain(pysam_fh, tmprss2)
        self.insert_chain(pysam_fh, erg)
        
        # Merge arcs somehow, label nodes
        
        # emperical evidence showed ~230bp? look into this by picking a few examples
        #c.prune(400+126-12)
        # max obs = 418 for now
        thicker_arcs = self.chain.prune(PRUNE_INS_SIZE) # Makes arc thicker by lookin in the ins. size
        self.chain.merge_splice_juncs(SPLICE_JUNC_ACC_ERR)
        thicker_arcs = self.chain.rejoin_splice_juncs(thicker_arcs, PRUNE_INS_SIZE) # Merges arcs by splice junctions and other junctions
        self.chain.reinsert_arcs(thicker_arcs)
        subnets = self.chain.extract_subnetworks(thicker_arcs)
        
        # If circos:
        #s = 1
        #for subnet in subnets:
        #    c = CircosController(str(s), subnet, "tmp/circos.conf","tmp/select-coordinates.conf", "tmp/circos-data.txt")
        #    c.draw_network("tmp/test.png","tmp/test.svg")
        #    s += 1
        
        subnets = self.filter_subnets_on_identical_nodes(subnets)
        self.results = subnets
        
        return len(self.results)
    
    def test_disco_alignment(self,alignment_file):
        # Make sure the header exists indicating that the BAM file was
        # fixed using Dr. Disco
        bam_fh = pysam.AlignmentFile(alignment_file, "rb")
        if bam_fh.header.has_key('PG'):
            for pg in bam_fh.header['PG']:
                if pg['ID'] == 'drdisco_fix_chimeric':
                    return bam_fh
        
        raise Exception("Invalid STAR BAM File: has to be post processed with 'dr-disco fix-chimeric ...' first")
    
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
    
    
    def insert_chain(self,pysam_fh, region):
        for r in pysam_fh.fetch(region[0],region[1],region[2]):
        #for r in pysam_fh.fetch(region[0],region[1]):
            sa = self.parse_SA(r.get_tag('SA'))
            _chr = pysam_fh.get_reference_name(r.reference_id)
            insert_size = self.get_insert_size([_chr,r.reference_start],[sa[0][0],sa[0][1]])
            
            pos1 = None
            pos2 = None
            
            if r.get_tag('RG') == 'discordant_mates':
                if abs(insert_size) >= MIN_DISCO_INS_SIZE:
                    pos1, pos2 = self.chain.insert(r,sa[0])
            
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
                
                #@todo silent mates do not make pairs but are one directional
                # in their input - has to be fixed in order to allow arc_merging
                #self.chain.insert(r,broken_mate)
            
            elif r.get_tag('RG') in [
                'spanning_paired_1',
                'spanning_paired_1_r',
                'spanning_paired_1_s',
                'spanning_paired_1_t',
                'spanning_paired_2',
                'spanning_paired_2_r',
                'spanning_paired_2_s',
                'spanning_paired_2_t',
                'spanning_singleton_1',
                'spanning_singleton_1_r',
                'spanning_singleton_2',
                'spanning_singleton_2_r']:
                read = r
                pos1, pos2 = self.chain.insert(r,sa[0])

            else:
                raise Exception("Unknown type read: '"+str(r.get_tag('RG'))+"'. Was the alignment fixed with a more up to date version of Dr.Disco?")
            
            """Find introns etc:
            
            Example 1: 
            
            5S10M15N10M2S
            
            S S S S S | | | | | | | | |---------------| | | | | | | | | S S
            
 soft-clip: <========]
                                                                       [==>
            Softclip are usually one-directional, from the sequenced base until
            the end. If it is before the first M/= flag, it's direction should
            be '-', otherwise '+' as it clips from the end. In principle the
            soft/hard clips are not arcs but properties of nodes. One particular
            base may have several soft/hard clips (from different locations but
            adding weigt to the same node).
            
            Splice juncs are bi-directional, and real arcs.
 
splice-junc:                           <=============>

            """
            if pos1 != None and pos2 != None:
                for internal_arc in self.find_cigar_arcs(r):
                    #@todo return Arc object instead of tuple
                    if internal_arc[2] in ['cigar_splice_junction']:
                        i_pos1 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                             internal_arc[0],
                                             STRAND_FORWARD)
                        i_pos2 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                             internal_arc[1],
                                             STRAND_REVERSE)
                        
                    elif internal_arc[2] in ['cigar_soft_clip']:
                        # @todo add all _r's also
                        if r.get_tag('RG') in [ 'spanning_paired_1',
                                                'spanning_paired_1_r',
                                                'spanning_paired_1_s',
                                                'spanning_paired_1_t',
                                                'spanning_paired_2',
                                                'spanning_paired_2_r',
                                                'spanning_paired_2_s',
                                                'spanning_paired_2_t',
                                                'spanning_singleton_1',
                                                'spanning_singleton_1_r',
                                                'spanning_singleton_2',
                                                'spanning_singleton_2_r']:
                            i_pos1 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                                 internal_arc[0],
                                                 pos2.strand)
                            i_pos2 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                                 internal_arc[1],
                                                 pos1.strand)
                        
                        else:
                            raise Exception("what todo here?")
                    else:
                        raise Exception("Arc type not implemented: %s", internal_arc)
                    
                    if internal_arc[2] in ['cigar_soft_clip', 'cigar_hard_clip']:
                        try:
                            self.chain.get_node_reference(i_pos2).add_clip()
                        except:
                            # This happens with some weird reads
                            #if r.get_tag('RG') in ['spanning_paired_2','spanning_singleton_2_r']:
                            #self.chain.create_node(i_pos2)
                            #self.chain.get_node_reference(i_pos2).add_clip()
                            self.logger.warn("softclip of "+r.qname+" ("+r.get_tag('RG')+") of dist="+str(i_pos2.pos - i_pos1.pos)+" is not at the side of the break point.")
                            #else:
                            #    raise Exception("\n\nNode was not found, cigar correctly parsed?\n----------\ninternal arc: "+str(internal_arc)+"\n----------\nqname: "+r.qname+"\ncigar: "+r.cigarstring+"\ntype:  "+r.get_tag('RG')+"\npos:   "+str(pos2))
                    else:
                        self.chain.insert_entry(i_pos1,i_pos2,internal_arc[2],None,True)
    
    
    def export(self, fh):
        fh.write(str(self))
    
    def __str__(self):
        out = "chr-A\t"
        out += "pos-A\t"
        out += "direction-A\t"
        
        out += "chr-B\t"
        out += "pos-B\t"
        out += "direction-B\t"
        
        out += "score\t"
        out += "soft+hardclips\t"
        out += "n-split-reads\t"
        out += "n-discordant-reads\t"
        
        out += "n-arcs\t"
        out += "n-nodes-A\t"
        out += "n-nodes-B\t"
        
        out += "entropy-bp-arc\t"
        out += "entropy-all-arcs\t"
        
        out += "data-structure\n"
        
        for subnet in self.results:
            out += str(subnet)
        
        return out
    
    def filter_subnets_on_identical_nodes(self, subnets):
        new_subnets = []
        _id = 0
        #1. filter based on nodes - if there are shared nodes, exclude sn
        all_nodes = set()
        for i in range(len(subnets)):
            subnet = subnets[i]
            rmme = False
            score = 0
            nodes = set()
            for dp in subnet:
                score += dp[0].get_scores()
                
                nodes.add(dp[0]._origin)
                nodes.add(dp[0]._target)
            
            clips = 0
            for n in nodes:
                if n in all_nodes:
                    rmme = True
                else:
                    all_nodes.add(n)
                
                clips += n.clips
            
            if rmme:
                subnets[i] = None
            else:
                i += 1
                sn = Subnet(_id, subnet)
                sn.total_clips = clips
                sn.total_score = score
                new_subnets.append(sn)
        
        return new_subnets

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
        solid = False
        
        for chunk in read.cigar:
            if chunk[0] in tt.keys():# D N S H
                # Small softclips occur all the time - introns and 
                # deletions of that size shouldn't add weight anyway
                
                """
                1S 5M means:
                startpoint-1 , startpoint => Softclip 
                
                5M 2S means:
                startpoint +5, startpoint +5+2 => softclip
                
                In the first case, which I call left_clipping, the junction
                occurs before the start position and it needs to be corrected
                for
                
                @todo it's still a buggy implementation because the following
                is theoretically possible too:
                
                5H 10S 5M
                
                in that case, the first arc should be -15,-10 and the
                second -10,0. Maybe soft- and hard clipping should
                be merged together?
                """
                
                if chunk[0] in [4,5] and not solid:
                    offset -= chunk[1]

                if chunk[1] > 3:
                    if solid:
                        """Clips to the first node:
                        M M M M M M M M M M S S S S S
                                           <========]
                        """
                        #yield (offset , (offset + chunk[1]) , tt[chunk[0]])
                        yield ((offset + chunk[1]), offset , tt[chunk[0]])
                    else:
                        """Clips to the second node:
                        S S S S S M M M M M M M M M M
                        [========>
                        """
                        yield (offset , (offset + chunk[1]) , tt[chunk[0]])

            
            if chunk[0] not in [4,5]:
                solid = True
            
            offset += chunk[1]
