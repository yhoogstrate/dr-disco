#!/usr/bin/env python2
# *- coding: utf-8 -*-
# -- vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4
# https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

"""
Tries to find genomic and exon to exon break points within a discordant
RNA-Seq read alignment.
"""

#http://www.samformat.info/sam-format-flag
import logging,re,math,copy,sys,operator
import pysam
import HTSeq
from .CigarAlignment import *
from .CircosController import *

from fuma.Fusion import STRAND_FORWARD, STRAND_REVERSE, STRAND_UNDETERMINED 
strand_tt = {STRAND_FORWARD:'+',STRAND_REVERSE:'-',STRAND_UNDETERMINED:'?'}

# load cfg
from . import *

def entropy(frequency_table):
    n = sum(frequency_table.values())
    prob = [float(x)/n for x in frequency_table.values()]
    
    entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob if p > 0 ])
    if entropy <= 0:
        return 0.0
    else:
        return entropy / (math.log(n) / math.log(2.0))

def merge_frequency_tables(frequency_tables):
    """
    Merges different frequency tables.
    
    Parameters
    ----------
    frequency_tables: list
        [{'key1': int}, {'key1': int, 'key2: int}]
    
    Returns
    -------
    list
        {'key1': int, 'key2: int}
    """
    new_frequency_table = {}
    for t in frequency_tables:
        for key in t.keys():
            if not new_frequency_table.has_key(key):
                new_frequency_table[key] = 0
            
            new_frequency_table[key] += t[key]
    
    return new_frequency_table

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
        self._hash = self._chr.replace("chr","") + ("%0.2X" % self.pos).zfill(8) + strand_tt[self.strand]#http://stackoverflow.com/questions/38430277/python-class-hash-method-and-set
    
    def __lt__(self, other_bp):
        return self._hash < other_bp._hash
    
    def __str__(self):
        return str(self._chr)+":"+str(self.pos)+"/"+str(self.pos+1)+"("+strand_tt[self.strand]+")"
    
    def get_dist(self, other_bp, strand_specific):
        if not isinstance(other_bp, BreakPosition):# pragma: no cover
            raise Exception("Wrong data type used")
        
        if (not strand_specific or self.strand == other_bp.strand) and self._chr == other_bp._chr:
            return other_bp.pos - self.pos
        else:
            # chr1 has 249 000 000 bp - by replacing all digits by 9 is
            # must be larger than any Hs chr reflecting natural distances
            # in a way that interchromosomal breaks are 'larger' than
            # intrachromosomal ones
            return MAX_GENOME_DISTANCE


class Node:
    """A Simple node, just a chromomal location( and strand)"""

    def __init__(self,position):
        self.position = position
        self.clips = 0
        self.edges = {}
        self.splice_edges = {STRAND_FORWARD: {}, STRAND_REVERSE: {}}
    
    def __lt__(self, other_node):
        return self.position < other_node.position
    
    def get_connected_splice_junctions(self):
        """Recursively finds all nodes that are connected by splice junctions
        
        self: the node of which it's splice junctions are traversed in order to find more nodes that may not be already in the blak
        self.splice_edges:
         -  [     ] splice edge: splice p1, splice p2
        """    
        
        linked_nodes = set([self])# set(self) iterates over self.__iter__()
        linked_edges = {}
        
        self.get_connected_splice_junctions_recursively(linked_nodes, linked_edges, MAX_ACCEPTABLE_INSERT_SIZE, STRAND_REVERSE)
        self.get_connected_splice_junctions_recursively(linked_nodes, linked_edges, MAX_ACCEPTABLE_INSERT_SIZE, STRAND_FORWARD)
        
        return linked_nodes, linked_edges
    
    def get_connected_splice_junctions_recursively(self, nodes, edges, insert_size_to_travel, direction):
        new_nodes = []
        for edge_n in self.splice_edges[direction]:
            if edge_n not in nodes:
                edge = self.splice_edges[direction][edge_n]
                dkey = insert_size_to_travel - edge[0]# Calculate new traversal size. If we start with isze=450 and the first SJ is 50 bp away for the junction, we need to continue with 450-50=400
                if dkey >= 0:
                    new_nodes.append((edge_n, dkey))
                    
                    if not edges.has_key(edge_n):
                        edges[edge_n] = set()
                    edges[edge_n].add(edge[1])#use min() to consistsently use the one with the lowest mem addr - this only works if counts are used because otherwise the order may become dependent
    
        # old results, + recursive results
        for edge_n in new_nodes:
            nodes.add(edge_n[0])
        
        # only recusively add to the new ones
        for edge_n in new_nodes:
            edge_n[0].get_connected_splice_junctions_recursively(nodes, edges, edge_n[1], direction)
    
    def add_clip(self):
        self.clips += 1
    
    def get_edge_to_node(self, target_node):
        key = target_node.position._hash
        if self.edges.has_key(key):
            return self.edges[key]
        else:
            return None
    
    def insert_edge(self, edge):
        if edge._target == self:
            self.edges[edge._origin.position._hash] = edge
        else:
            self.edges[edge._target.position._hash] = edge
    
    def remove_edge(self,edge):
        try:
            del(self.edges[edge._origin.position._hash])
        except:
            del(self.edges[edge._target.position._hash])
    
    def __iter__(self):
        for k in sorted(self.edges.keys()):
            yield self.edges[k]
    
    def __str__(self):# pragma: no cover
        out  = str(self.position)
        
        a = 0
        for sedge in self.edges:
            edge = self.edges[sedge]
            filtered_edges = {JunctionTypeUtils.str(x):edge._types[x] for x in sorted(edge._types.keys())}# if x not in ['cigar_soft_clip','cigar_hard_clip']
            len_edges = len(filtered_edges)
            a += len_edges
            if len_edges > 0:
                out += "\n\t["+str(id(edge))+":"+sedge+"] "+str(edge._origin.position)+"->"+str(edge._target.position)+" "+str(filtered_edges)
        
        if a > 0:
            out += "\n\t-> soft/hard clips: "+str(self.clips)
            
            return out+"\n"
        else:
            return ""


class JunctionTypes:# Enum definition to avoid string operations
    silent_mate            = 1
    cigar_soft_clip        = 2
    cigar_hard_clip        = 3
    cigar_deletion         = 4
    cigar_splice_junction  = 5
    discordant_mates       = 6
    spanning_paired_1      = 7
    spanning_paired_1_r    = 8
    spanning_paired_1_s    = 9
    spanning_paired_1_t    = 10
    spanning_paired_2      = 11
    spanning_paired_2_r    = 12
    spanning_paired_2_s    = 13
    spanning_paired_2_t    = 14
    spanning_singleton_1   = 15
    spanning_singleton_1_r = 16
    spanning_singleton_2   = 17
    spanning_singleton_2_r = 18


class JunctionTypeUtils:
    tt = {
        'silent_mate': JunctionTypes.silent_mate,
        'discordant_mates': JunctionTypes.discordant_mates,
        'spanning_paired_1': JunctionTypes.spanning_paired_1,
        'spanning_paired_2': JunctionTypes.spanning_paired_2,
        'spanning_paired_1_r': JunctionTypes.spanning_paired_1_r,
        'spanning_paired_2_r': JunctionTypes.spanning_paired_2_r,
        'spanning_paired_1_s': JunctionTypes.spanning_paired_1_s,
        'spanning_paired_2_s': JunctionTypes.spanning_paired_2_s,
        'spanning_paired_1_t': JunctionTypes.spanning_paired_1_t,
        'spanning_paired_2_t': JunctionTypes.spanning_paired_2_t,
        'spanning_singleton_1': JunctionTypes.spanning_singleton_1,
        'spanning_singleton_2': JunctionTypes.spanning_singleton_2,
        'spanning_singleton_1_r': JunctionTypes.spanning_singleton_1_r,
        'spanning_singleton_2_r': JunctionTypes.spanning_singleton_2_r }
    
    scoring_table = {
        JunctionTypes.discordant_mates: 1,
        JunctionTypes.spanning_paired_1: 3,
        JunctionTypes.spanning_paired_1_r: 3,
        JunctionTypes.spanning_paired_1_s: 1,# Very odd type of reads
        JunctionTypes.spanning_paired_1_t: 3,
        JunctionTypes.spanning_paired_2: 3,
        JunctionTypes.spanning_paired_2_r: 3,
        JunctionTypes.spanning_paired_2_s: 1,# Very odd type of reads
        JunctionTypes.spanning_paired_2_t: 3,
        JunctionTypes.spanning_singleton_1: 2, 
        JunctionTypes.spanning_singleton_2: 2,
        JunctionTypes.spanning_singleton_1_r: 2,
        JunctionTypes.spanning_singleton_2_r: 2 }
    
    @staticmethod
    def enum(enumstring):
        return JunctionTypeUtils.tt[enumstring]
    
    @staticmethod
    def str(enumcode):
        for key, value in JunctionTypeUtils.tt.items():
            if value == enumcode:
                return key
    
    @staticmethod
    def get_score(junction_type):
        return JunctionTypeUtils.scoring_table[junction_type]
    
    @staticmethod
    def is_fusion_junction(junction_type):
        return junction_type >= 6


class Edge:
    """Connection between two genomic locations
     - the connection can by of different types
    """
    
    def __init__(self,_origin,_target):
        if not isinstance(_target, Node) or not isinstance(_origin, Node):# pragma: no cover
            raise Exception("_origin and _target must be a Node")
        
        self._origin = _origin
        self._target = _target
        self._types = {}
        self._unique_alignment_hashes = {}# Used for determining entropy
    
    def get_entropy(self):
        """Entropy is an important metric/ propery of an edge
        The assumption is that all reads should be 'different', i.e.
        have a different start position, different end position,
        different soft/hardclips etc.
        
        It is more probable to find the following alignment for a true
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
        
        return entropy(self._unique_alignment_hashes)
    
    def get_complement(self):
        try:
            return self._target.edges[self._origin.position._hash]
        except KeyError as err:#todo write test for this
            raise KeyError("Could not find complement for edge:   "+str(self))
    
    def merge_edge(self,edge):
        """Merges (non splice) edges"""
        
        for alignment_key in edge._unique_alignment_hashes:
            self.add_alignment_key(alignment_key)
        
        for _type in edge._types:
            if _type != JunctionTypes.silent_mate:
                self.add_type(_type,edge._types[_type])
    
    def add_type(self, _type, weight):
        if not self._types.has_key(_type):
            self._types[_type] = weight
        else:
            self._types[_type] += weight
    
    def add_alignment_key(self,alignment_key):
        if not self._unique_alignment_hashes.has_key(alignment_key):
            self._unique_alignment_hashes[alignment_key] = 1
        else:
            self._unique_alignment_hashes[alignment_key] += 1
    
    def get_count(self, _type):
        if not self._types.has_key(_type):
            return 0
        else:
            return self._types[_type]
        
    def get_score(self,_type):  
        return self.get_count(_type)*JunctionTypeUtils.scoring_table[_type]
    
    def get_scores(self):
        """Based on this function, the start point is determined
        
        Currently it's implemented as sum(split reads)
        
        Later on spanning reads will be added, softclips will be added, etc.
        """
        return sum([self.get_score(_type) for _type in JunctionTypeUtils.scoring_table.keys()])
    
    def get_splice_score(self):# pragma: no cover
        return (self.get_count(JunctionTypes.cigar_splice_junction),self.get_clips())
    
    def get_clips(self):# pragma: no cover
        return self._origin.clips + self._target.clips

    def __str__(self):
        typestring = []
        for _t in sorted(self._types.keys()):
            typestring.append(str(JunctionTypeUtils.str(_t))+":"+str(self._types[_t]))
        
        out = str(self._origin.position) + "->" + str(self._target.position) + ":("+','.join(typestring)+")"
        
        #spacer = " "*len(str(self._origin.position))
        
        #for _k in self._origin.splice_edges.keys():
        #    out += "\n"+spacer+"=>"+str(_k.position)+":score="+str(self._origin.splice_edges[_k][1].get_splice_score())
        
        return out


class Graph:
    def __init__(self):
        self.idxtree = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    
    def __iter__(self):
        for step in self.idxtree.steps():
            position = step[1]
            if position:
                for strand in sorted(position.keys()):
                    node = position[strand]
                    yield node
    
    def create_node(self,pos):
        position_accession = HTSeq.GenomicPosition(pos._chr,pos.pos)
        position = self.idxtree[position_accession]
        if not position:
            self.idxtree[position_accession] = {pos.strand: Node(pos)}
        else:
            if not position.has_key(pos.strand):
                position[pos.strand] = Node(pos)
    
    def get_node_reference(self,pos):
        position = self.idxtree[HTSeq.GenomicPosition(pos._chr,pos.pos)]
        if position and position.has_key(pos.strand):
            return position[pos.strand]
        return None
    
    def insert_edge(self,pos1,pos2,_type,cigarstrs):
        """ - Checks if Node exists at pos1, otherwise creates one
            - Checks if Node exists at pos2, otherwise creates one
            - Checks if Edge exists between them, otherwise inserts it into the Nodes
         
         cigarstrs must be something like ("126M","126M") or ("25S50M2000N","25M50S")
        """
        self.create_node(pos1)
        self.create_node(pos2)
        
        node1 = self.get_node_reference(pos1)
        node2 = self.get_node_reference(pos2)
        
        if cigarstrs != None:
            cigarstrs = pos1._hash+cigarstrs[0]+"|"+pos2._hash+cigarstrs[1]
        
        edge = node1.get_edge_to_node(node2)
        if edge == None:
            edge = Edge(node1, node2)
            node1.insert_edge(edge)
            node2.insert_edge(edge)
        
        edge.add_type(_type, 1)
        if node1 == edge._origin:# Avoid double insertion of all keys :) only do it if the positions don't get swapped
            edge.add_alignment_key(cigarstrs)
    
    def insert_splice_edge(self,pos1,pos2,_type,cigarstrs):
        """ - Checks if Node exists at pos1, otherwise creates one
            - Checks if Node exists at pos2, otherwise creates one
            - Checks if Edge exists between them, otherwise inserts it into the Nodes
         
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
        
            cigarstrs = short_pos1+strand_tt[pos1.strand]+cigarstrs[0]+"|"+short_pos2+strand_tt[pos2.strand]+cigarstrs[1]
        
        edge1 = node1.get_edge_to_node(node2)
        
        if edge1 == None:
            edge1 = Edge(node1, node2)
            edge2 = Edge(node2, node1)
            
            node1.insert_edge(edge1)
            node2.insert_edge(edge2)
        else:
            edge2 = node2.get_edge_to_node(node1)
        
        edge1.add_type(_type, 1)
        edge2.add_type(_type, 1)
        #if node1 == edge._origin:# Avoid double insertion of all keys :) only do it if the positions don't get swapped
        #    edge.add_alignment_key(cigarstrs)
    
    def reinsert_edges(self, edges):
        """Only works for Edges of which the _origin and _target Node
        still exists
        """
        for edge in edges:
            node1 = edge._origin
            node2 = edge._target
            
            node1.insert_edge(edge)
            node2.insert_edge(edge)
    
    def remove_edge(self, edge):
        node1 = edge._origin
        node2 = edge._target
        
        node1.remove_edge(edge)
        node2.remove_edge(edge)
        
        del(edge)
    
    def remove_node(self,node):
        self.idxtree[HTSeq.GenomicPosition(node.position._chr, node.position.pos)] = set()
        del(node)
    
    def search_splice_edges_between(self,pos1,pos2):# insert size is two directional
        for step in self.idxtree[HTSeq.GenomicInterval(pos1._chr, max(0, pos1.pos - MAX_ACCEPTABLE_INSERT_SIZE), pos1.pos + MAX_ACCEPTABLE_INSERT_SIZE + 1)].steps():
            if step[1]:
                position = step[1]
                if position:
                    for strand in position.keys():
                        node1 = position[strand]
                        for edge in node1.edges.values():
                            node2 = edge._target
                            
                            d1 = abs(pos1.pos - node1.position.pos)
                            if d1 <= MAX_ACCEPTABLE_INSERT_SIZE:
                                d2 = abs(pos2.pos - node2.position.pos)
                                if d1+d2 <= MAX_ACCEPTABLE_INSERT_SIZE:
                                    yield pow(d1, 2) + pow(d2, 2) , edge
    
    def print_chain(self):# pragma: no cover
        print "**************************************************************"
        for node in self:
            _str = str(node)
            if len(_str.strip()) > 0:
                print _str
        print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    
    def generate_edge_idx(self):
        logging.info("Creating edge index before pruning")
        
        edges = set()
        edges_tuple = []
        order = 0
        
        for node in self:
            for edge in node.edges.values():
                if edge not in edges:
                    score = edge.get_scores()
                    if score > 0:
                        edges.add(edge)
                        edges_tuple.append((edge,score,order))
                        order -= 1
        
        del(edges,order)
        
        self.edge_idx = [edge[0] for edge in sorted(edges_tuple, key=operator.itemgetter(1, 2), reverse=False)]
    
    def prune(self):
        """Does some 'clever' tricks to merge edges together and reduce data points
        """
        self.generate_edge_idx()
        logging.info("Finding and merging other edges in close proximity (insert size)")
        
        candidates = []
        #self.print_chain()
        
        while self.edge_idx:
            candidate = self.edge_idx.pop()
            
            self.prune_edge(candidate)
            candidates.append(candidate)
            
            self.remove_edge(candidate)# do not remove if splice junc exists?
        
        #self.print_chain()
        logging.info("Pruned into "+str(len(candidates))+" candidate edge(s)")
        return candidates
    
    def prune_edge(self, edge):
        ## @ todo double check if this is strand specific
        
        for edge_m in self.search_edges_between(edge):
            d1 = edge._origin.position.get_dist(edge_m._origin.position, True)
            d2 = edge._target.position.get_dist(edge_m._target.position, True)
            d = abs(d1) + abs(d2)
            
            if d <= MAX_ACCEPTABLE_INSERT_SIZE:
                edge.merge_edge(edge_m)
                
                self.remove_edge(edge_m)
                self.edge_idx.remove(edge_m)
    
    def search_edges_between(self,edge_to_prune):
        """searches for other junctions in-between edge+insert size:"""
        def pos_to_range(pos):
            if pos.strand == STRAND_REVERSE:
                return (pos.pos - MAX_ACCEPTABLE_INSERT_SIZE) - 1, pos.pos + MAX_ACCEPTABLE_ALIGNMENT_ERROR
            else:
                return pos.pos - MAX_ACCEPTABLE_ALIGNMENT_ERROR, (pos.pos + MAX_ACCEPTABLE_INSERT_SIZE) + 1
        
        pos1, pos2 = edge_to_prune._origin.position, edge_to_prune._target.position
        
        pos1_min, pos1_max = pos_to_range(pos1)
        pos2_min, pos2_max = pos_to_range(pos2)
        
        for step in self.idxtree[HTSeq.GenomicInterval(pos1._chr,pos1_min,pos1_max + 1)].steps():
            if step[1]:
                position = step[1]
                if position and position.has_key(pos1.strand):
                    node_i = position[pos1.strand]
                    for edge in node_i.edges.values():
                        if edge != edge_to_prune and edge._target.position.strand == pos2.strand and edge._target.position.pos >= pos2_min and edge._target.position.pos <= pos2_max:
                            yield edge
    
    def rejoin_splice_juncs(self, splice_junctions):
        """thicker edges go across the break point:

thick edges:
                                ------------------------------------------
                       ---------------------------------------------------
                       ----------------------------------------------------------
              -------------------------------------------------------------------

 splice juncs:
               -------  -------                                            -----
             ||       ||       ||        |          $ ... $       |      ||     ||

        
        the goal is to add the splice juncs between the nodes
        """
        logging.debug("Initiated")
        k = 0
        
        def search(pos1):
            #@todo make this member of Graph and use splice_junctions.search_..._..(pos)
            nodes = []
            for step in self.idxtree[HTSeq.GenomicInterval(pos1._chr,max(0,pos1.pos - MAX_ACCEPTABLE_INSERT_SIZE), pos1.pos + MAX_ACCEPTABLE_INSERT_SIZE + 1)].steps():
                if step[1]:
                    position = step[1]
                    if position:
                        for strand in position.keys():
                            nodes.append(position[strand])
            return nodes
        
        splice_edges_had = set()
        
        for node in splice_junctions:
            for splice_junction in node.edges.values():
                if splice_junction not in splice_edges_had:
                    lnodes = search(splice_junction._origin.position)
                    rnodes = search(splice_junction._target.position)
                    
                    for lnode in lnodes:
                        for rnode in rnodes:
                            if lnode != rnode and lnode.position.strand == rnode.position.strand:
                                d1 = abs(splice_junction._origin.position.get_dist(lnode.position, False))
                                d2 = abs(splice_junction._target.position.get_dist(rnode.position, False))
                                dist = d1+d2
                                if dist <= MAX_ACCEPTABLE_INSERT_SIZE:
                                    if lnode.position < rnode.position and lnode.splice_edges[STRAND_FORWARD].has_key(rnode):
                                        insert = dist < lnode.splice_edges[STRAND_FORWARD][rnode][0]
                                    elif rnode.position < lnode.position and lnode.splice_edges[STRAND_REVERSE].has_key(rnode):
                                        insert = dist < lnode.splice_edges[STRAND_REVERSE][rnode][0]
                                    else:
                                        insert = True
                                    
                                    if insert:
                                        if rnode.position < lnode.position:# inversed
                                            rnode.splice_edges[STRAND_FORWARD][lnode] = [dist, splice_junction]
                                            lnode.splice_edges[STRAND_REVERSE][rnode] = [dist, splice_junction]
                                        else:
                                            lnode.splice_edges[STRAND_FORWARD][rnode] = [dist, splice_junction]
                                            rnode.splice_edges[STRAND_REVERSE][lnode] = [dist, splice_junction]
                                        
                                        k += 1
                    
                    splice_edges_had.add(splice_junction)
                    splice_edges_had.add(splice_junction.get_complement())
        
        logging.info("Linked "+str(k)+" splice junction(s)")
    
    def extract_subnetworks_by_splice_junctions(self,thicker_edges):
        """ Deze functie haalt recursief per edge een set van edges op die 
binnen de splice afstand voldoen. De splice afstanden zijn voorgerekend
in de rejoin_splice_junctions functie. Het enige dat hier dient te gebeuren
is het recursief aflopen aan de hand van een bepaalde afstand (insert size,
~450bp).
Vooralsnog is deze extractie alleen gebaseerd op afstand en niet op de
imbalans van de readcount. 
In de nieuwe implementatie heeft iedere node twee vectoren waarin splice
junction staan beschreven; 1 naar posities die kleiner zijn dan zichzelf
en 1 met posities die groter zijn dan zichzelf. Hierdoor is een recursief
terugloop probleem redelijk opgelost.
            """
        logging.info("Initiated")
        
        thicker_edges.reverse()
        q = 0
        subnetworks = []
        while thicker_edges:
            start_point = thicker_edges.pop()
            
            left_nodes, left_splice_junctions_ds = start_point._origin.get_connected_splice_junctions()
            right_nodes, right_splice_junctions_ds = start_point._target.get_connected_splice_junctions()
            
            left_splice_junctions = set()
            right_splice_junctions = set()
            
            subedges = []
            
            # All edges between any of the left and right nodes are valid edges and have to be extracted
            # Find all direct edges joined by splice junctions
            for left_node in left_nodes:
                for right_node in right_nodes:
                    if left_node.edges.has_key(right_node.position._hash):
                        subedge = left_node.edges[right_node.position._hash]
                        if subedge._target.position._hash == right_node.position._hash:
                            subedges.append(subedge)
                            
                            if left_node != start_point._origin:# must be merged by a splice junction
                                left_splice_junctions = left_splice_junctions.union(left_splice_junctions_ds[left_node])
                            
                            if right_node != start_point._target:
                                right_splice_junctions = right_splice_junctions.union(right_splice_junctions_ds[right_node])

                            if subedge != start_point:
                                self.remove_edge(subedge)
                                thicker_edges.remove(subedge)
            
            del(left_splice_junctions_ds,right_splice_junctions_ds)
            subnetworks.append(Subnet(q,subedges,left_splice_junctions,right_splice_junctions))
        
        logging.info("Extracted %i subnetwork(s)" % len(subnetworks))
        return subnetworks


class Subnet():
    def __init__(self,_id,edges,left_splice_junctions,right_splice_junctions):
        self._id = _id
        self.edges = edges
        self.left_splice_junctions = left_splice_junctions
        self.right_splice_junctions = right_splice_junctions
        self.total_clips = 0
        self.total_score = 0
        self.discarded = []
        
        self.calc_clips()
        self.calc_scores()
        
        self.reorder_edges()
    
    def reorder_edges(self):
        idx = {}
        for edge in self.edges:
            key1 = edge.get_scores()
            key2 = edge._origin.position._hash+"-"+edge._target.position._hash
            
            if not idx.has_key(key1):
                idx[key1] = {}
            
            idx[key1][key2] = edge
        
        ordered = []
        for key1 in sorted(idx.keys(),reverse=True):
            for key2 in sorted(idx[key1].keys()):
                ordered.append(idx[key1][key2])
        
        self.edges = ordered
    
    def calc_clips(self):
        clips = 0
        
        nodes = set()
        for edge in self.edges:
            nodes.add(edge._origin)
            nodes.add(edge._target)
        
        for node in nodes:
            clips += node.clips
        
        self.total_clips = clips
        return self.total_clips
    
    def calc_scores(self):
        self.total_score = 0
        for edge in self.edges:
            self.total_score += edge.get_scores()

        return self.total_score
    
    def __str__(self):
        """Makes tabular output"""
        node_a,  node_b  = self.edges[0]._origin, self.edges[0]._target
        nodes_a, nodes_b = self.get_n_nodes()
        
        return (
            "%s\t%i\t%s\t"
            "%s\t%i\t%s\t"
            "%s\t"
            "%i\t%i\t%i\t%i\t"
            "%i\t%i\t%i\t"
            "%i\t%i\t"
            "%s\t%s\t"#%.2f\t%.2f\t
            "%s\n" % (
                    node_a.position._chr, node_a.position.pos, strand_tt[self.edges[0]._origin.position.strand], # Pos-A
                    node_b.position._chr, node_b.position.pos, strand_tt[self.edges[0]._target.position.strand], # Pos-B
                    ("valid" if self.discarded == [] else ','.join(self.discarded)), # Classification status
                    self.total_score/2, self.total_clips, self.get_n_split_reads()/2, self.get_n_discordant_reads()/2, # Evidence stats
                    len(self.edges), nodes_a, nodes_b, # Edges and nodes stats
                    len(self.left_splice_junctions), len(self.right_splice_junctions),
                    self.edges[0].get_entropy(), self.get_overall_entropy(), # Entropy stats
                    "&".join([str(edge) for edge in self.edges]) # Data structure
            ))
    
    def get_overall_entropy(self):
        frequency_table = merge_frequency_tables([edge._unique_alignment_hashes for edge in self.edges])
        return entropy(frequency_table)
    
    def get_n_nodes(self):
        nodes_a, nodes_b = set(), set()
        
        for edge in self.edges:
            nodes_a.add(edge._origin)
            nodes_b.add(edge._target)
        
        return len(nodes_a), len(nodes_b)
    
    def get_n_split_reads(self):
        n = 0
        
        for _type in [JunctionTypes.spanning_paired_1,      JunctionTypes.spanning_paired_2,
                      JunctionTypes.spanning_paired_1_r,    JunctionTypes.spanning_paired_2_r,
                      JunctionTypes.spanning_paired_1_s,    JunctionTypes.spanning_paired_2_s,
                      JunctionTypes.spanning_paired_1_t,    JunctionTypes.spanning_paired_2_t,
                      JunctionTypes.spanning_singleton_1,   JunctionTypes.spanning_singleton_2,
                      JunctionTypes.spanning_singleton_1_r, JunctionTypes.spanning_singleton_2_r]:
            for edge in self.edges:
                n += edge.get_count(_type)
        
        return n
    
    def get_n_discordant_reads(self):
        return sum([edge.get_count(JunctionTypes.discordant_mates) for edge in self.edges])
    
    def find_distances(self, subnet_t):
        """ - Must be symmetrical i.e. (sn1.find_distances(s2) == s2.find_distances(s1)
        if the distance of any edge == infty, return False
        """
        
        l_nodes_min = set(x for x in self.get_lnodes())
        l_nodes_max = set(x for x in subnet_t.get_lnodes())
        
        r_nodes_min = set(x for x in self.get_rnodes())
        r_nodes_max = set(x for x in subnet_t.get_rnodes())
        
        if len(l_nodes_min) >  len(l_nodes_max):
            l_nodes_min, l_nodes_max = l_nodes_max, l_nodes_min
        
        if len(r_nodes_min) >  len(r_nodes_max):
            r_nodes_min, r_nodes_max = r_nodes_max, r_nodes_min
        
        l_dists = []
        r_dists = []
        
        for l_node in l_nodes_min:# makes it symmetrical
            dist = MAX_GENOME_DISTANCE
            
            for l_node_t in l_nodes_max:
                dist = min(abs(l_node.position.get_dist(l_node_t.position, True)),dist)
            
            if dist == MAX_GENOME_DISTANCE:
                return False, False
            else:
                l_dists.append(dist)
    
        for r_node in r_nodes_min:# makes it symmetrical
            dist = MAX_GENOME_DISTANCE
        
            for r_node_t in r_nodes_max:
                dist = min(abs(r_node.position.get_dist(r_node_t.position, True)),dist)
            
            if dist == MAX_GENOME_DISTANCE:
                return False, False
            else:
                r_dists.append(dist)
        
        return l_dists, r_dists
    
    def get_lnodes(self):
        lnodes = set()
        
        for edge in self.edges:
            lnodes.add(edge._origin)
        
        for lnode in lnodes:
            yield lnode
    
    def get_rnodes(self):
        rnodes = set()
        
        for edge in self.edges:
            rnodes.add(edge._target)
        
        for rnode in rnodes:
            yield rnode
    
    def merge(self, subnet_m):
        for edge in subnet_m.edges:
            self.edges.append(edge)
        
        self.left_splice_junctions = self.left_splice_junctions.union(subnet_m.left_splice_junctions)
        self.right_splice_junctions = self.right_splice_junctions.union(subnet_m.right_splice_junctions)
        
        self.calc_clips()
        self.calc_scores()
        self.reorder_edges()


class BAMExtract(object):
    def __init__(self, bam_file):
        self.pysam_fh = self.test_disco_alignment(bam_file)
    
    @staticmethod
    def test_disco_alignment(alignment_file):
        """Ensures by reading the BAM header whether the BAM file was
        indeed fixed using Dr. Disco
        """
        bam_fh = pysam.AlignmentFile(alignment_file, "rb")
        if bam_fh.header.has_key('PG'):
            for pg in bam_fh.header['PG']:
                if pg['ID'] == 'drdisco_fix_chimeric':
                    try:# pragma: no cover
                        bam_fh.fetch()
                    except:# pragma: no cover
                        logging.info('Indexing BAM file with pysam: '+bam_fh.filename)# create index if it does not exist
                        pysam.index(bam_fh.filename)
                        bam_fh = pysam.AlignmentFile(bam_fh.filename)
                    
                    try:
                        bam_fh.fetch()
                    except:
                        raise Exception('Could not indexing BAM file: '+bam_fh.filename)
                    
                    return bam_fh
        
        #@todo write simple test
        raise Exception("Invalid STAR BAM File: has to be post processed with 'dr-disco fix-chimeric ...' first")
    
    def extract_junctions(self, fusion_junctions, splice_junctions):
        def read_to_junction(read, rg, parsed_SA_tag, specific_type = None):
            pos1, pos2 = None, None
            
            if rg in [JunctionTypes.discordant_mates]:
                # How to distinguish between:
                #
                # Type 1:
                # ===1===>|   |<===2===
                
                # Type 2:
                # ===2===>|   |<===1===
                
                # Hypothetical:
                #|<===1===    |<===2===
                
                if read.is_reverse:
                    pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                         read.reference_start,
                                         STRAND_FORWARD)
                else:
                    pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                         read.reference_start+bam_parse_alignment_offset(read.cigar),
                                         STRAND_REVERSE)
                
                if parsed_SA_tag[4] == "-":
                    pos2 = BreakPosition(parsed_SA_tag[0],
                                         parsed_SA_tag[1],
                                         STRAND_FORWARD)
                else:
                    pos2 = BreakPosition(parsed_SA_tag[0],
                                         bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                         STRAND_REVERSE)
            
            elif rg in [JunctionTypes.spanning_singleton_1, JunctionTypes.spanning_paired_1]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1],
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

            elif rg in [JunctionTypes.spanning_singleton_1_r, JunctionTypes.spanning_paired_1_t]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start,
                                     not read.is_reverse)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            elif rg in [JunctionTypes.spanning_paired_1_r]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                     not read.is_reverse)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1],
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
            
            elif rg in [JunctionTypes.spanning_singleton_2, JunctionTypes.spanning_paired_2]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start,
                                     read.is_reverse)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            elif rg in [JunctionTypes.spanning_singleton_2_r, JunctionTypes.spanning_paired_2_t]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                    (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                     read.is_reverse)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1],
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            elif rg in [JunctionTypes.spanning_paired_2_r]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start,
                                     read.is_reverse)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

            elif rg in [JunctionTypes.spanning_paired_1_s]:
                # Very clear example in S054 @ chr21:40,064,610-40,064,831  and implemented in test case 14
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)
                
                pos2 = BreakPosition(parsed_SA_tag[0],
                         parsed_SA_tag[1],
                         STRAND_REVERSE if parsed_SA_tag[4] == "+" else STRAND_FORWARD)
             
            elif rg in [JunctionTypes.spanning_paired_2_s]:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 STRAND_FORWARD if read.is_reverse else STRAND_REVERSE)
            
                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
            
            elif rg not in [JunctionTypes.silent_mate]:# pragma: no cover
                raise Exception("Fatal Error, RG: %s" % JunctionTypeUtils.str(rg))
            
            else:
                raise Exception("Unnknown read group: %s", rg)
            
            if abs(pos1.get_dist(pos2,False)) >= MAX_ACCEPTABLE_INSERT_SIZE:
                return (pos1, pos2, (read.cigarstring,parsed_SA_tag[2]))
            
            return (None,None,None)
        
        logging.debug("Parsing reads to obtain fusion gene and splice junctions")
        for read in self.pysam_fh.fetch():
            sa = self.parse_SA(read.get_tag('SA'))
            _chr = self.pysam_fh.get_reference_name(read.reference_id)
            rg = JunctionTypeUtils.enum(read.get_tag('RG'))
            
            pos1, pos2 = None, None
            
            if JunctionTypeUtils.is_fusion_junction(rg):
                pos1, pos2, junction = read_to_junction(read, rg, sa[0])
                if pos1 != None:
                    fusion_junctions.insert_edge(pos1, pos2, rg, junction)
            
            elif rg == JunctionTypes.silent_mate:# Not yet implemented, may be useful for determining type of junction (exonic / intronic)
                pass

            else:# pragma: no cover
                raise Exception("Unknown type read: %s. Was the alignment fixed with a more up to date version of Dr.Disco?" % JunctionTypeUtils.str(rg))
            
            # Detect introns, clipping and insertions/deletions by SAM flags
            for internal_edge in self.find_cigar_edges(read):
                i_pos1, i_pos2 = None, None
                
                if internal_edge[2] in [JunctionTypes.cigar_splice_junction, JunctionTypes.cigar_deletion]:#, JunctionType.cigar_deletion
                    i_pos1 = BreakPosition(_chr, internal_edge[0], STRAND_FORWARD)
                    i_pos2 = BreakPosition(_chr, internal_edge[1], STRAND_REVERSE)
                
                    if internal_edge[2] in [JunctionTypes.cigar_deletion] and i_pos1.get_dist(i_pos2,False) < MAX_ACCEPTABLE_INSERT_SIZE:
                        i_pos1, i_pos2 = None, None
                
                elif internal_edge[2] in [JunctionTypes.cigar_soft_clip]:
                    if pos1 == None or rg in [JunctionTypes.spanning_paired_1_s, JunctionTypes.spanning_paired_2_s]:
                        pass
                    elif JunctionTypeUtils.is_fusion_junction(rg):
                        i_pos1 = BreakPosition(_chr, internal_edge[0], pos2.strand)
                        i_pos2 = BreakPosition(_chr, internal_edge[1], pos1.strand)
                    else:# pragma: no cover
                        raise Exception("what todo here %s "+JunctionTypeUtils.str(rg))
                else:# pragma: no cover
                    raise Exception("Edge type not implemented: %s\n\n%s", internal_edge,str(read))
                
                if internal_edge[2] in [JunctionTypes.cigar_soft_clip, JunctionTypes.cigar_hard_clip]:
                    if i_pos1 != None:
                        node = fusion_junctions.get_node_reference(i_pos2)
                        if node != None:
                            node.add_clip()
                        # else condition happens when clips are found on positions that are not junctions
                else:
                    if i_pos1 != None:
                        if internal_edge[2] == JunctionTypes.cigar_splice_junction:
                            #splice_junctions.insert_splice_edge(i_pos1,i_pos2,internal_edge[2],None)
                            splice_junctions.insert_edge(i_pos1,i_pos2,internal_edge[2],None)
                        else:
                            fusion_junctions.insert_edge(i_pos1,i_pos2,internal_edge[2],None)
        
        logging.debug("alignment data loaded")
    
    def parse_pos(self, str_pos):
        _chr,_poss = str_pos.split(":",2)
        _poss = _poss.replace(",","").split("-",2)
        
        return str(_chr), int(_poss[0]), int(_poss[1])
    
    def extract(self, pos1, pos2, bamfile_out):
        pos1 = self.parse_pos(pos1)
        pos2 = self.parse_pos(pos2)
        
        ids = set()
        
        fh = pysam.AlignmentFile(bamfile_out, "wb", header=self.pysam_fh.header)
        
        for r in self.pysam_fh.fetch(pos1[0],pos1[1],pos1[2]):#pysam.view(b,"chr21:1-500")
            ids.add(r.query_name)
        
        for r in self.pysam_fh.fetch(pos2[0],pos2[1],pos2[2]):#pysam.view(b,"chr21:1-500")
            ids.add(r.query_name)
        
        for r in self.pysam_fh.fetch():
            if r.query_name in ids:
                fh.write(r)
        
        
        fh.close()
        
        pysam.index(str(bamfile_out))
    
    # static is elegant for functions that do not use class properties
    # http://stackoverflow.com/questions/18679803/python-calling-method-without-self
    @staticmethod
    def parse_SA(SA_tag):
        sa_tags = SA_tag.split(";")
        for i in xrange(len(sa_tags)):
            sa_tags[i] = sa_tags[i].split(",")
            sa_tags[i][1] = int(sa_tags[i][1])
        
        return sa_tags
    
    @staticmethod
    def find_cigar_edges(read):
        """Tries to find Edges introduced by:
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
            X	BAM_CDIFF	8
            
            Example 1: 
            
            5S10M15N10M2S
            
            S S S S S | | | | | | | | |---------------| | | | | | | | | S S
            
 soft-clip: <========]
                                                                       [==>
            Softclip are usually one-directional, from the sequenced base until
            the end. If it is before the first M/= flag, it's direction should
            be '-', otherwise '+' as it clips from the end. In principle the
            soft/hard clips are not edges but properties of nodes. One particular
            base may have several soft/hard clips (from different locations but
            adding weigt to the same node).
            
            Splice juncs are bi-directional, and real edges.
 
splice-junc:                           <=============>
            """
        
        tt = {
            2:JunctionTypes.cigar_deletion,
            3:JunctionTypes.cigar_splice_junction,
            4:JunctionTypes.cigar_soft_clip,
            5:JunctionTypes.cigar_hard_clip
        }
        
        offset = read.reference_start
        solid = False
        
        for chunk in read.cigar:
            if chunk[0] in tt.keys():# D N S H
                """Small softclips occur all the time - introns and 
                deletions of that size shouldn't add weight anyway
                
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
                
                in that case, the first edge should be -15,-10 and the
                second -10,0. Maybe soft- and hard clipping should
                be merged together?
                """
                
                if chunk[0] in [4,5] and not solid:
                    offset -= chunk[1]

                if chunk[1] > MAX_ACCEPTABLE_ALIGNMENT_ERROR:
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


class IntronDecomposition:
    def __init__(self,alignment_file):
        self.alignment_file = alignment_file
    
    def decompose(self):
        alignment = BAMExtract(self.alignment_file)
        
        fusion_junctions = Graph()
        splice_junctions = Graph()
        
        alignment.extract_junctions(fusion_junctions, splice_junctions)
        
        thicker_edges = fusion_junctions.prune() # Makes edge thicker by lookin in the ins. size - make a sorted data structure for quicker access - i.e. sorted list
        fusion_junctions.rejoin_splice_juncs(splice_junctions) # Merges edges by splice junctions and other junctions
        del(splice_junctions)
        fusion_junctions.reinsert_edges(thicker_edges)#@todo move function into statuc function in this class
        
        subnets = fusion_junctions.extract_subnetworks_by_splice_junctions(thicker_edges)
        subnets = self.merge_overlapping_subnets(subnets)
        self.results = self.filter_subnets(subnets)# Filters based on three rules: entropy, score and background
        
        # If circos:
        #s = 1
        #for subnet in self.results:
        #    c = CircosController(str(s), subnet, "tmp/circos.conf","tmp/select-coordinates.conf", "tmp/circos-data.txt")
        #    c.draw_network("tmp/test.png","tmp/test.svg")
        #    s += 1
        
        return len(self.results)
    
    def export(self, fh):
        fh.write(str(self))
    
    def __str__(self):
        ordered = []
        order = 0 # Order of top-edge in subnet
        for subnet in self.results:
            ordered.append((subnet, subnet.total_score, subnet.get_overall_entropy()))
        ordered = [subnet[0] for subnet in sorted(ordered, key=operator.itemgetter(1, 2), reverse=True)]
        
        return (
            "chr-A"           "\t" "pos-A"             "\t" "direction-A""\t"
            "chr-B"           "\t" "pos-B"             "\t" "direction-B""\t"
            "filter-status"   "\t"
            "score"           "\t" "soft+hardclips"    "\t" "n-split-reads" "\t" "n-discordant-reads" "\t"
            "n-edges"         "\t" "n-nodes-A"         "\t" "n-nodes-B"     "\t"
            "n-splice-junc-A" "\t" "n-splice-junc-B"   "\t"
            "entropy-bp-edge" "\t" "entropy-all-edges" "\t"
            "data-structure"  "\n"
            "%s" % (''.join([str(subnet) for subnet in ordered]) )
            )

    def merge_overlapping_subnets(self, subnets):
        """Merges very closely adjacent subnets based on the smallest
        internal distance. E.g. if we have a subnet having 1 and one
        having 2 edges:
        
  snA:  |       |            ~             |
        |        --------------------------  edgeA1
         ----------------------------------  edgeA2

  snB         |                        |
               ------------------------      edgeB1

        We would like to filter based on the distance between edgeA1 and
        edgeB1.
        
        We also don't want to have a growth pattern, i.e. that based on
        merging snB to snA, snC becomes part of it. Although it may be
        true it can become problematc as I've seen with FuMa.
        
        So, idea is:
        loop over all subnets i;
            loop over all subnets j > i
                if subnet j should be merged with i, remeber it into M
            
            merge all subnets in M into i, and remove the former subnets
        """
        logging.info("initiated")
        
        def sq_dist(vec):
            sum_of_squares = sum(pow(x,2) for x in vec)
            
            if sum_of_squares > 0:
                avg_sq_d = float(sum_of_squares) / len(vec)
            else:
                avg_sq_d = 0.0
            
            return math.sqrt(avg_sq_d)
        
        def tree_insert(genometree, pos, entry):
            """inserts entry into a set at pos in genometree"""
            position_accession = HTSeq.GenomicPosition(pos._chr, pos.pos)
            position_tree = genometree[position_accession]
            if not position_tree:# position was not found in tree, insert completely
                genometree[position_accession] = {pos.strand: set([entry])}
            else:
                if not position_tree.has_key(pos.strand):
                    position_tree[pos.strand] = set([entry])
                else:
                    position_tree[pos.strand].add(entry)
        
        def tree_remove(genometree, subnet):
            positions = set()
            for edge in subnet.edges:
                positions.add(edge._origin.position)
                positions.add(edge._target.position)
            
            for pos in positions:
                position_accession = HTSeq.GenomicPosition(pos._chr, pos.pos)
                position_tree = genometree[position_accession]
                if position_tree and position_tree.has_key(pos.strand):
                    sset = position_tree[pos.strand]
                    if subnet in sset:
                        sset.remove(subnet)
        
        n = len(subnets)
        
        idx_l = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        idx_r = HTSeq.GenomicArrayOfSets("auto", stranded=False)

        for subnet in subnets:
            for edge in subnet.edges:
                tree_insert(idx_l, edge._origin.position, subnet)
                tree_insert(idx_r, edge._target.position, subnet)
        
        k = 0
        new_subnets = []
        subnets.reverse()
        while subnets:
            subnet = subnets.pop()
            tree_remove(idx_l, subnet)
            tree_remove(idx_r, subnet)
        
            to_be_merged_with = set()
            
            for edge in subnet.edges:
                
                candidates = set()
                
                # based on l-node to be close
                pos = edge._origin.position
                for step in idx_l[HTSeq.GenomicInterval(pos._chr, max(0, pos.pos - MAX_ACCEPTABLE_INSERT_SIZE), pos.pos + MAX_ACCEPTABLE_INSERT_SIZE)].steps():
                    if step[1] and step[1].has_key(pos.strand):
                        for candidate_subnet in step[1][pos.strand]:
                            if subnet != candidate_subnet:
                                candidates.add(candidate_subnet)
                
                pos = edge._target.position
                for step in idx_r[HTSeq.GenomicInterval(pos._chr, pos.pos - MAX_ACCEPTABLE_INSERT_SIZE, pos.pos + MAX_ACCEPTABLE_INSERT_SIZE)].steps():
                    if step[1] and step[1].has_key(pos.strand):
                        for candidate_subnet in step[1][pos.strand]:
                            if subnet != candidate_subnet:
                                candidates.add(candidate_subnet)
                
                for candidate_subnet in candidates:
                    new_merged = False
                    
                    l_dists, r_dists = subnet.find_distances(candidate_subnet)
                    
                    if l_dists != False:
                        
                        n_l_dist = sum([1 for x in l_dists if x < MAX_ACCEPTABLE_INSERT_SIZE])
                        n_r_dist = sum([1 for x in r_dists if x < MAX_ACCEPTABLE_INSERT_SIZE])
                        
                        rmsq_l_dist = sq_dist(l_dists)
                        rmsq_r_dist = sq_dist(r_dists)
                        
                        """ @todo work with polynomial asymptotic equasion based on rmse, product and k, determine a True or False
                            if product > (8*k):
                                # Average gene size = 10-15kb
                                # rmse <= 125000 - 120000/2^( product / 1200)
                                # if rmse < max_rmse: valid data point
                                max_rmse = 125000.0 - (120000/pow(2, float(product) / 1200.0))
                                return (rmse <= max_rmse)
                        """
                        
                        if n_l_dist > 0 and n_r_dist > 0:
                            new_merged = True
                        else:#@todo revise if / else / elif logic, and use linear regression or sth like that
                            if n_l_dist > 0:
                                l_dist_ins_ratio = 1.0 * n_l_dist/len(l_dists)
                                if l_dist_ins_ratio == 1.0 and rmsq_r_dist < 15000:
                                    new_merged = True
                                elif l_dist_ins_ratio > 0.7 and rmsq_r_dist < 10000:
                                    new_merged = True
                                elif l_dist_ins_ratio > 0.3 and rmsq_r_dist < 5000:
                                    new_merged = True
                                        
                            if n_r_dist > 0:
                                r_dist_ins_ratio = 1.0 * n_r_dist/len(r_dists)
                                if r_dist_ins_ratio == 1.0 and rmsq_l_dist < 15000:
                                    new_merged = True
                                elif r_dist_ins_ratio > 0.7 and rmsq_l_dist < 10000:
                                    new_merged = True
                                elif r_dist_ins_ratio > 0.3 and rmsq_l_dist < 5000:
                                    new_merged = True
                        
                        if new_merged:
                            to_be_merged_with.add(candidate_subnet)
                    
            for candidate_subnet in to_be_merged_with:
                subnet.merge(candidate_subnet)
                tree_remove(idx_l, candidate_subnet)
                tree_remove(idx_r, candidate_subnet)
                subnets.remove(candidate_subnet)#[subnets.index(candidate_subnet)] = None# @todo inverse subnets and use .pop() and .remove()
                del(candidate_subnet)
                k += 1
                
            new_subnets.append(subnet)
        
        logging.info("Merged "+str(k)+" of the "+str(n)+" into "+str(len(new_subnets))+" merged subnetwork(s)")
        
        return new_subnets

    def filter_subnets(self, subnets):
        logging.debug("init")
        k = 0
        for subnet in subnets:
            """Total of 8 reads is minimum, of which 2 must be
            discordant and the entropy must be above 0.55"""
            
            entropy = subnet.get_overall_entropy()
            if entropy < MIN_SUBNET_ENTROPY:
                subnet.discarded.append("entropy="+str(entropy))
            
            n_disco = subnet.get_n_discordant_reads() / 2
            n_disco_min = MIN_DISCO_PER_SUBNET_PER_NODE * sum(subnet.get_n_nodes())
            if n_disco < n_disco_min:
                subnet.discarded.append("n_discordant_reads="+str(n_disco)+"/"+str(n_disco_min))
            
            n_support = (subnet.get_n_discordant_reads() + subnet.get_n_split_reads() ) / 2
            n_support_min = (MIN_SUPPORTING_READS_PER_SUBNET_PER_NODE * sum(subnet.get_n_nodes()))
            if n_support < n_support_min:
                subnet.discarded.append("n_support="+str(n_support)+"/"+str(n_support_min))
            
            if len(subnet.discarded) > 0:
                k += 1
        
        logging.info("Filtered "+str(k)+" of the "+str(len(subnets))+" subnetwork(s)")
        return subnets
