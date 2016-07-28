#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Tries to figure out within a discordant RNA-Seq read alignment:
 - Tries to explain introns and exons by fusion transcripts
"""

#http://www.samformat.info/sam-format-flag

import logging,re,math
import pysam,copy
from intervaltree_bio import GenomeIntervalTree, Interval
from .CigarAlignment import *

from fuma.Fusion import STRAND_FORWARD, STRAND_REVERSE, STRAND_UNDETERMINED 


class Arc:
    """Connection between two genomic locations
     - Also contains different types of evidence
    """
    scoring_table={
               'cigar_splice_junction': 0,
               'discordant_mates': 2,
               'silent_mate': 0,
               'spanning_paired_1': 3,
               'spanning_paired_2': 3,
               'spanning_singleton_1': 2, 
               'spanning_singleton_2': 2,
               'spanning_singleton_1_r': 1,
               'spanning_singleton_2_r': 1
               }
    
    def __init__(self,_origin,_target):
        if not isinstance(_target, Node) or not isinstance(_origin, Node):
            raise Exception("_origin and _target must be a Node")
        
        #self._strand = STRAND_UNDETERMINED
        
        self._origin = _origin
        self._target = _target
        self._types = {}
    
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
        for _type in arc._types:
            if _type in [
                    "discordant_mates",
                    "spanning_paired_2",
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
        if not self._types.has_key(_type):
            self._types[_type] = 0
        
        self._types[_type] += 1
    
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
        return (self.get_count('cigar_splice_junction'),0)
    
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
        
        spacer = " "*len(str(self._origin.position))
        
        for _k in self._origin.splice_arcs.keys():
            out += "\n"+spacer+"=>"+str(_k)+",score:"+str(self._origin.splice_arcs[_k][1].get_splice_score())
        
        return out


class Node:
    """A Simple node, just a chromomal location( and strand)"""
    
    def __init__(self,position):
        self.position = position
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
    
    def insert_arc(self,arc,arc_type):
        skey = str(arc._target.position)
        
        if not self.arcs.has_key(skey):
            self.set_arc(arc)
        
        self.arcs[skey].add_type(arc_type)
    
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
    
    def new_arc(self,node2,arc_type,do_vice_versa):
        if do_vice_versa:
            node2.new_arc(self,arc_type,False)
        
        arc = Arc(self,node2)
        self.insert_arc(arc,arc_type)
    
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
        sc = 0
        hc = 0
        
        for sarc in self.arcs:
            arc = self.arcs[sarc]
            filtered_arcs = {x:arc._types[x] for x in sorted(arc._types.keys()) if x not in ['cigar_soft_clip','cigar_hard_clip']}
            len_arcs = len(filtered_arcs)
            a += len_arcs
            if len_arcs > 0:
                out += "\n\t-> "+str(arc._target.position)+" "+str(filtered_arcs)
            
            if arc._types.has_key('cigar_soft_clip'):
                sc += arc._types['cigar_soft_clip']

            if arc._types.has_key('cigar_hard_clip'):
                hc += arc._types['cigar_hard_clip']
        
        #if (a+sc+hc) > 0:
        if a > 0:
            out += "\n\t-> incoming soft-clips: "+str(sc)
            out += "\n\t-> incoming hard-clips: "+str(hc)
            
            return out+"\n"
        else:
            return ""

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


class Chain:
    def __init__(self,pysam_fh):
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
        try:
            return list(self.idxtree[pos._chr][pos.pos])[0][2][pos.strand]
        except:
            return None
    
    def insert_entry(self,pos1,pos2,_type,do_vice_versa):
        """
         - Checks if Node exists at pos1, otherwise creates one
         - Checks if Node exists at pos2, otherwise creates one
         - Checks if Arc exists between them
        """
        
        self.create_node(pos1)
        self.create_node(pos2)
        
        node1 = self.get_node_reference(pos1)
        node2 = self.get_node_reference(pos2)
        
        node1.new_arc(node2,_type,do_vice_versa)
    
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
        if rg in ["discordant_mates","silent_mate"]:
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
            
            self.insert_entry(pos1,pos2,rg,False)
        
        elif rg in ["spanning_paired_1","spanning_singleton_1"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 not read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,False)

        elif rg in ["spanning_singleton_1_r"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 not read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,False)
        
        elif rg in ["spanning_paired_2","spanning_singleton_2"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,False)
        
        elif rg in ["spanning_singleton_2_r"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            self.insert_entry(pos1,pos2,rg,False)
        
            #else:
            #pos2 = BreakPosition(parsed_SA_tag[0],
            #                    bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
            #                    STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
        
            #else:
                #if not read.is_reverse:
                    #pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                        #read.reference_start,
                                        #not read.is_reverse)
                #else:
                    #pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                        #bam_parse_alignment_end(read),
                                        #not read.is_reverse)
                    
                
                ### The read is the second in pair, and the SA tag the first...
                #if parsed_SA_tag[4] != "+":
                    #print read
                    #raise Exception("Todo")
                    ## Most likely code that should do it:
                    ##
                    ##pos1 = BreakPosition(parsed_SA_tag[0],
                    ##                    parsed_SA_tag[1],
                    ##                    STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
                    ##
                    ##self.insert_entry(pos1,pos2,rg)
                    
                #else:
                    #pos1 = BreakPosition(parsed_SA_tag[0],
                                        #bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                        #STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
        
        else:
            raise Exception("Fatal Error, RG: "+rg)
    
    def reinsert_arcs(self, arcs):
        """Only works for Arcs of which the _origin and _target Node
        still exists
        """
        for arc_t in arcs:
            arc   = arc_t[0]
            arc_c = arc_t[2]
            
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
            comp_strand = STRAND_REVERSE
        else:
            comp_strand = STRAND_FORWARD
        
        if pos.strand == comp_strand:
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
        """
        
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
        candidates = []
        candidate = self.get_start_point()
        
        #print candidate
        self.print_chain()
        
        i = 1
        while candidate != None:
            if i > 15:
                raise Exception("Recusion depth errr")
            
            ratio = self.prune_arc(candidate, insert_size)
            candidates.append((candidate, ratio, candidate.get_complement()))
            
            self.remove_arc(candidate)
            
        #    self.print_chain()
            
            candidate = self.get_start_point()
             
            i += 1
        
        
        self.print_chain()
        
        #for candidate in candidates:
        #    print candidate
        
        return candidates
    
    #@todo ADD RECURSION DEPTH LIMIT
    def prune_arc(self, arc, insert_size):
        ## @todo make somehow only 
        #if 
        node1 = arc._origin
        node2 = arc._target
        
        ratio = self.arcs_ratio_between(node1.position, node2.position, insert_size)
        
        for c_arc in self.search_arcs_between(node1.position, node2.position, insert_size):
            ## These types of arcs fall within the same space / inst. size
            ## they make the arcs heavier
            arc_complement = node2.arcs[str(node1.position)]
            arc.merge_arc(c_arc)
            arc.merge_arc(arc_complement)
            
            self.remove_arc(c_arc)
            
            #self.remove_arc(arc_complement)
            #if arc.merge_arc(c_arc) or arc_complement.merge_arc(c_arc):
            #    try:
            #        self.remove_arc(c_arc)
            #    except:
            #        self.print_chain()
            #        raise Exception("Error")
        
        """
        for c_arc in self.search_arcs_adjacent(node1.position, node2.position, insert_size):
            # @todo see if we can merge by other types of Arcs too..
            if c_arc[1].get_count('cigar_splice_junction') > 0:
                #if c_arc[0] == 'pos1':
                #    print "merge it with node1"
                #elif c_arc[0] == 'pos2':
                #    print "merge it with node2"
                #else:
                #    raise Exception("Unknown type: %s", c_arc[0])
                
                #complement_arc = self.get_node_reference(c_arc[1]._target.position).arcs[str(c_arc[1]._origin.position)]
                
                arc.new_arc(c_arc[1])
                #arc.new_arc(complement_arc)
                
                #arc_complement.new_arc(c_arc[1])
                
                self.remove_arc(c_arc[1])# Removes bi-directional
                
                #self.prune_arc(c_arc[1], insert_size)
                #self.prune_arc(complement_arc, insert_size)
                
                # @todo pruning needs to be done also to node2 instead
                # of to the 2 splice junctions?
                
                
                # if recursion depth not reached:
                # self.prune_arc(c_arc[1],isize)
                # else:
                # raise exception max rec. depth reached
        """
        
        return ratio
    
    def merge_splice_juncs(self,uncertainty):
        #self.print_chain()
        
        init = self.get_start_splice_junc()
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
        
        #added_splice_junctions = set([])
        
        left_nodes = set([])
        right_nodes = set([])
        
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
        subnetworks = []
        
        
        start_point = thicker_arcs[0][0]
        
        left_nodes = [start_point._origin]
        right_nodes = [start_point._target]
        
        
        print left_nodes
        ## The original nodes have been emptied, so the most important
        ## arc's are now separated.
        
        ## Let's add Node.insert() inserting Arcs directly
        
        left_nodes = start_point._origin.rfind_connected_sjuncs(left_nodes)
        right_nodes = start_point._target.rfind_connected_sjuncs(right_nodes)
        
        # All arcs between any of the left and right nodes are valid arcs and have to be extracted
        print left_nodes
        
        # Find all direct arcs joined by splice junctions
        for left_node in left_nodes:
            for right_node in right_nodes:
                if left_node.arcs.has_key(str(right_node.position)):
                    print "Found arc:", (left_node.arcs[str(right_node.position)],right_node.arcs[str(left_node.position)])
        
        
        
        ## Find all with one indirect step - these might be alternative junctions / exons
        i = -1
        for left_node_i in left_nodes:
            j = -1
            i += 1
            
            for left_node_j in left_nodes:
                j += 1
                
                if i < j:
                    #print i,j,[left_node_i],"*",[left_node_j]
                    ## Find similar destinations
                    #print left_node_i.arcs.keys()
                    #print left_node_j.arcs.keys()
                    mutual = set(left_node_i.arcs.keys()).intersection(left_node_j.arcs.keys())
                    print 
        
        return subnetworks


class IntronDecomposition:
    def __init__(self,break_point):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        self.break_point = break_point
        self.chain = None
    
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
            
            if r.get_tag('RG') == 'discordant_mates':
                if abs(insert_size) >= 400:
                    self.chain.insert(r,sa[0])
            
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
            
            elif r.get_tag('RG') in ['spanning_paired_1', 'spanning_paired_2', 'spanning_singleton_1', 'spanning_singleton_2', 'spanning_singleton_1_r', 'spanning_singleton_2_r']:
                read = r
                self.chain.insert(r,sa[0])

            else:
                raise Exception("Unknown type read: '"+str(r.get_tag('RG'))+"'. Was the alignment fixed with a more up to date version of Dr.Disco?")
            
            """Find introns etc:
            
            Example 1: 
            
            5S10M15N10M2S
            
            S S S S S | | | | | | | | |---------------| | | | | | | | | S S
            
 soft-clip: <========]
                                                                       [==>
            Softclip are usually one directional, from the sequenced base until
            the end. If it is before the first M/= flag, it's direction should
            be '-', otherwise '+'
 
splice-junc:                           <=============>

            """
            for internal_arc in self.find_cigar_arcs(r):
                #@todo return Arc object instead of tuple
                if internal_arc[2] in ['cigar_splice_junction']:
                    pos1 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                         internal_arc[0],
                                         STRAND_FORWARD)
                    pos2 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                         internal_arc[1],
                                         STRAND_REVERSE)
                    
                elif internal_arc[2] in ['cigar_soft_clip']:
                    if r.get_tag('RG') in ['spanning_paired_2', 'spanning_singleton_2']:
                        pos1 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                             internal_arc[0],
                                             r.is_reverse)
                        pos2 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                             internal_arc[1],
                                             r.is_reverse)
                    else:
                        pos1 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                             internal_arc[0],
                                             not r.is_reverse)
                        pos2 = BreakPosition(pysam_fh.get_reference_name(r.reference_id),
                                             internal_arc[1],
                                             not r.is_reverse)
                else:
                    raise Exception("Arc type not implemented: %s", internal_arc)
                
                self.chain.insert_entry(pos1,pos2,internal_arc[2],True)
    
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
        thicker_arcs = self.chain.prune(450) # Makes arc thicker by lookin in the ins. size
        self.chain.merge_splice_juncs(3)
        thicker_arcs = self.chain.rejoin_splice_juncs(thicker_arcs, 450) # Merges arcs by splice junctions and other junctions
        self.chain.reinsert_arcs(thicker_arcs)
        subnets = self.chain.extract_subnetworks(thicker_arcs)
        
        return thicker_arcs

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
        left_clipping = True
        
        for chunk in read.cigar:
                          # D N S H
            if chunk[0] in tt.keys():
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
                if chunk[0] in [4,5] and left_clipping:
                    offset -= chunk[1]
                
                if chunk[1] > 3:
                    yield (offset , (offset + chunk[1]) , tt[chunk[0]])
            
            if chunk[0] not in [4,5]:
                left_clipping = False
            
            offset += chunk[1]
