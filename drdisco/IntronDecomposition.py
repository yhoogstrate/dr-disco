#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:
# https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

"""
Tries to find genomic and exon to exon break points within a discordant
RNA-Seq read alignment.
"""

#http://www.samformat.info/sam-format-flag

import logging,re,math,copy,sys
import pysam
from intervaltree_bio import GenomeIntervalTree, Interval
from .CigarAlignment import *
from .CircosController import *
from .BAMExtract import BAMExtract


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


class Edge:
    """Connection between two genomic locations
     - the connection can by of different types
    """
    
    scoring_table={
               'cigar_splice_junction': 0,
               'discordant_mates': 1,
               'silent_mate': 0,
               'spanning_paired_1': 3,
               'spanning_paired_1_r': 3,
               'spanning_paired_1_s': 1,# Very odd type of reads
               'spanning_paired_1_t': 3,
               'spanning_paired_2': 3,
               'spanning_paired_2_r': 3,
               'spanning_paired_2_s': 1,# Very odd type of reads
               'spanning_paired_2_t': 3,
               'spanning_singleton_1': 2, 
               'spanning_singleton_2': 2,
               'spanning_singleton_1_r': 2,
               'spanning_singleton_2_r': 2
               }
    
    def __init__(self,_origin,_target):
        if not isinstance(_target, Node) or not isinstance(_origin, Node):
            raise Exception("_origin and _target must be a Node")
        
        #self._strand = STRAND_UNDETERMINED
        
        self._origin = _origin
        self._target = _target
        self._types = {}
        
        # Used for determining entropy
        self.unique_alignments_idx = {}
    
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
        
        return entropy(self.unique_alignments_idx)
    
    def get_complement(self):
        """complement = self._target.edges[str(self._origin.position)]
        
        if complement == self:
            raise Exception("err")
        
        return complement
        """
        try:
            return self._target.edges[str(self._origin.position)]
        except KeyError as err:
            raise KeyError("Could not find complement for edge:   "+str(self))
    
    def merge_edge(self,edge):
        """Merges discordant_mates edges
        """
        
        for alignment_key in edge.unique_alignments_idx:
            self.add_alignment_key(alignment_key)
        
        for _type in edge._types:
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
            raise Exception("Clips shouldn't be added as edges, but as properties of Nodes")
        
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
        
        #for _k in self._origin.splice_edges.keys():
        #    out += "\n"+spacer+"=>"+str(_k.position)+":score="+str(self._origin.splice_edges[_k][1].get_splice_score())
        
        return out


class Node:
    """A Simple node, just a chromomal location( and strand)"""
    
    def __init__(self,position):
        self.position = position
        self.clips = 0
        self.edges = {}
        self.splice_edges = {}
    
    def rfind_connected_sjuncs(self,left_nodes,insert_size_to_travel):
        """Recursively finds all nodes that are connected by splice junctions
        
        self: the node of which it's splice junctions are traversed in order to find more nodes that may not be already in the blak
        self.splice_edges:
         -  [     ] splice edge: splice p1, splice p2
         
        """
        
        #@todo only lnodes are taken into account? / not sure about this; name left_nodes is confusing
        print "Rfind d=",insert_size_to_travel,left_nodes
        print 'self:',[self]
        results_new2 = {}
        for edge_n in self.splice_edges.keys():
            if edge_n not in left_nodes:
                edge = self.splice_edges[edge_n]
                print '  * edge_n:',[edge_n]
                
                # distance between the node self and the nodes of the exon junction
                ds1 = abs(self.position.get_dist(edge[1]._target.position,False))
                ds2 = abs(self.position.get_dist(edge[1]._origin.position,False))# Consider strand specificness... possible with current graph model?
                # distance between the target node and the nodes of the exon junction
                dt1 = abs(edge_n.position.get_dist(edge[1]._target.position,False))
                dt2 = abs(edge_n.position.get_dist(edge[1]._origin.position,False))# Consider strand specificness... possible with current graph model?
                # this has to be crossed i.e. if the one node is close to the left part of the SJ the other node must be close to the right part
                d1 = ds1 + dt2
                d2 = ds2 + dt1
                d = min(d1,d2)
                print "  (ds1,ds2:",ds1,ds2,")(dt1,dt2:",dt1,dt2,"): [d1 ,d2]",d1, d2
                
                if d <= insert_size_to_travel:
                    print "  YES, we can!"
                    dkey = insert_size_to_travel - d# Calculate new traversal size. If we start with isze=450 and the first SJ is 50 bp away for the junction, we need to continue with 450-50=400
                    if not results_new2.has_key(dkey):
                        results_new2[dkey] = set()
                    results_new2[dkey].add(edge_n)
                else:
                    print "NO: ",d,">",insert_size_to_travel
        
        # old results, + recursive results
        results_all = set(left_nodes)
        for depth in results_new2.keys():
            results_all = results_all.union(results_new2[depth])
        
        print ">>",results_all
        print
        
        # only recusively add to the new ones
        for depth in results_new2.keys():
            if depth > 0:
                for node in results_new2[depth]:
                    for edge in node.rfind_connected_sjuncs(results_all, depth):
                        results_all.add(edge)
        
        return list(results_all)
    
    def add_clip(self):
        self.clips += 1
    
    def insert_edge(self,edge,edge_type,alignment_key):
        skey = str(edge._target.position)
        
        if not self.edges.has_key(skey):
            self.set_edge(edge)
        
        self.edges[skey].add_type(edge_type)
        
        if alignment_key != None:# Splice junctions should be skipped for entropy
            self.edges[skey].add_alignment_key(alignment_key)
    
    def set_edge(self, edge):
        self.edges[str(edge._target.position)] = edge
    
    def get_top_edge(self):
        maxscore = -1
        top_edge = None
        
        for edge in self:
            score = edge.get_scores()
            if score > maxscore:
                maxscore = score
                top_edge = edge
        
        if top_edge != None:
            return (maxscore, top_edge, self, top_edge._target)
        else:
            return (None, None, None, None)
    
    def new_edge(self,node2,edge_type,alignment_key,do_vice_versa):
        if do_vice_versa:
            node2.new_edge(self,edge_type,alignment_key,False)
        
        edge = Edge(self,node2)
        self.insert_edge(edge,edge_type,alignment_key)
    
    def remove_edge(self,edge,idx):
        if idx == "by-target":
            skey = str(edge._target.position)
        elif idx == "by-origin":
            skey = str(edge._origin.position)
        else:
            raise Exception("Invalid usage of function")
        
        if not self.edges.has_key(skey):
            raise Exception("Unknown key: %s", skey)
        
        self.edges[skey] = None
        del(self.edges[skey])
    
    def __iter__(self):
        for k in sorted(self.edges.keys()):
            yield self.edges[k]
    
    def __str__(self):
        out  = str(self.position)
        
        a = 0
        for sedge in self.edges:
            edge = self.edges[sedge]
            filtered_edges = {x:edge._types[x] for x in sorted(edge._types.keys()) if x not in ['cigar_soft_clip','cigar_hard_clip']}
            len_edges = len(filtered_edges)
            a += len_edges
            if len_edges > 0:
                out += "\n\t-> "+str(edge._target.position)+" "+str(filtered_edges)

        
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
            p = nearby_node.edges[str(self.position)]
            scores.append(p.get_scores())
        
        k = len(scores)
        if k < 2:
            raise Exception("nearby_nodes contains less than 2 elements with shared Edges")
        
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
            return MAX_GENOMIC_DIST
    
    def get_rmse(self, pos_vec):
        err = 0
        for node in pos_vec:
            pos = node.position
            err += pow(self.get_dist(pos, True) , 2)
        return math.sqrt(err)


class Chain:
    def __init__(self,pysam_fh):
        self.idxtree = GenomeIntervalTree()
        #self.scoringtree = []
        self.pysam_fh = pysam_fh
    
    def __iter__(self):
        for key in self.idxtree:
            for element in sorted(self.idxtree[key]):
                for strand in sorted(element[2].keys()):
                    yield element[2][strand]

    def __iter_chr(self,key):
        for element in sorted(self.idxtree[key]):
            for strand in sorted(element[2].keys()):
                yield element[2][strand]
    
    def create_node(self,pos):
        """Creates the Node but does not overwrite it
        """
        # see if position exists in idx
        if len(self.idxtree[pos._chr][pos.pos]) == 0:
            self.idxtree[pos._chr].addi(pos.pos,pos.pos+1,{})
        
        # make sure the strand is added
        if not list(self.idxtree[pos._chr][pos.pos])[0][2].has_key(pos.strand):
            list(self.idxtree[pos._chr][pos.pos])[0][2][pos.strand] = Node(pos)
    
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
        
    def get_node_reference(self,pos):
        try:
            return list(self.idxtree[pos._chr][pos.pos])[0][2][pos.strand]
        except:
            return None
    
    def insert_alignment(self):
        logging.debug("starting inserting alignment data")
        for read in self.pysam_fh.fetch():
            sa = BAMExtract.BAMExtract.parse_SA(read.get_tag('SA'))
            _chr = self.pysam_fh.get_reference_name(read.reference_id)
            rg = read.get_tag('RG')
            
            pos1 = None
            pos2 = None

            if rg in [
                'discordant_mates',
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
                pos1, pos2 = self.insert(read,sa[0])
            
            elif rg == 'silent_mate':
                # usually in a exon?
                if read.is_read1:
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
                # in their input - has to be fixed in order to allow edge_merging
                #self.chain.insert(r,broken_mate)

            else:
                raise Exception("Unknown type read: '"+str(rg)+"'. Was the alignment fixed with a more up to date version of Dr.Disco?")
            
            """Find introns etc:
            
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
            
            for internal_edge in BAMExtract.BAMExtract.find_cigar_edges(read):
                if internal_edge[2] in ['cigar_splice_junction']:#, 'cigar_deletion'
                    #@todo _chr?
                    i_pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                           internal_edge[0],
                                           STRAND_FORWARD)
                    #@todo _chr?
                    i_pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                           internal_edge[1],
                                           STRAND_REVERSE)
                
                elif internal_edge[2] in ['cigar_deletion']:
                    #@todo _chr?
                    i_pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                           internal_edge[0],
                                           STRAND_FORWARD)
                    #@todo _chr?
                    i_pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                           internal_edge[1],
                                           STRAND_REVERSE)
                    
                    if i_pos1.get_dist(i_pos2,False) < MIN_DISCO_INS_SIZE:
                        i_pos1 = None
                        i_pos2 = None
                
                elif internal_edge[2] in ['cigar_soft_clip']:
                    if pos1 != None and pos2 != None and rg in ['discordant_mates',
                              'spanning_paired_1',
                              'spanning_paired_1_r',
                              'spanning_paired_1_t',
                              'spanning_paired_2',
                              'spanning_paired_2_r',
                              'spanning_paired_2_t',
                              'spanning_singleton_1',
                              'spanning_singleton_1_r',
                              'spanning_singleton_2',
                              'spanning_singleton_2_r']:
                        #@todo _chr?
                        i_pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                             internal_edge[0],
                                             pos2.strand)
                        #@todo _chr?
                        i_pos2 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                             internal_edge[1],
                                             pos1.strand)
                    
                    elif rg in ['spanning_paired_1_s',
                                'spanning_paired_2_s']:
                        i_pos1 = None
                        i_pos2 = None
                    else:
                        raise Exception("what todo here - "+rg)
                else:
                    raise Exception("Edge type not implemented: %s\n\n%s", internal_edge,str(read))
                
                if internal_edge[2] in ['cigar_soft_clip', 'cigar_hard_clip']:
                    try:
                        if i_pos1 != None:
                            self.get_node_reference(i_pos2).add_clip()
                    except:
                        # This happens with some weird reads
                        #if rg in ['spanning_paired_2','spanning_singleton_2_r']:
                        #self.chain.create_node(i_pos2)
                        #self.chain.get_node_reference(i_pos2).add_clip()
                        #logging.warn("softclip of "+read.qname+" ("+rg+") of dist="+str(i_pos2.pos - i_pos1.pos)+" is not at the side of the break point.")
                        pass
                else:
                    if i_pos1 != None:
                        self.insert_entry(i_pos1,i_pos2,internal_edge[2],None,True)
        
        logging.debug("alignment data loaded")
    
    def insert_entry(self,pos1,pos2,_type,cigarstrs,do_vice_versa):
        """ - Checks if Node exists at pos1, otherwise creates one
            - Checks if Node exists at pos2, otherwise creates one
            - Checks if Edge exists between them
         
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
            
            node1.new_edge(node2,_type,alignment_key,do_vice_versa)
        else:
            node1.new_edge(node2,_type,None,do_vice_versa)
    
    def insert(self,read,parsed_SA_tag,specific_type = None):
        """Inserts a bi-drectional edge between read and sa-tag in the Chain
        determines the type of edge by @RG tag (done by dr-disco fix-chimeric)

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


        This means that for:
         - spanning_paired_1 the breakpoint is at: start + offset
         - spanning_paired_2 the breakpoint is at: start
        regardless of the strand (+/-)"""
        
        pos1 = None
        pos2 = None
        
        rg = read.get_tag('RG')
        
        if rg in ["discordant_mates"]:
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
        
        elif rg in ["spanning_singleton_1", "spanning_paired_1"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

        elif rg in ["spanning_singleton_1_r", "spanning_paired_1_t"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 not read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
        
        elif rg in ["spanning_paired_1_r"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 not read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)
        
        elif rg in ["spanning_singleton_2", "spanning_paired_2"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
        
        elif rg in ["spanning_singleton_2_r", "spanning_paired_2_t"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1],
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
        
        elif rg in ["spanning_paired_2_r"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 read.reference_start,
                                 read.is_reverse)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

        elif rg in ["spanning_paired_1_s"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                 (read.reference_start+bam_parse_alignment_offset(read.cigar)),
                                 STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)
            
            pos2 = BreakPosition(parsed_SA_tag[0],
                     parsed_SA_tag[1],
                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
         
        elif rg in ["spanning_paired_2_s"]:
            pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                             read.reference_start,
                             not read.is_reverse)
        
            pos2 = BreakPosition(parsed_SA_tag[0],
                                 parsed_SA_tag[1] + bam_parse_alignment_offset(cigar_to_cigartuple(parsed_SA_tag[2])),
                                 STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)
        
        elif rg not in ["silent_mate"]:
            raise Exception("Fatal Error, RG: "+rg)
        
        if pos1 != None:# First check if insert size makes sense anyway (only for certain types of edges)
            if rg in ['discordant_mates',
                      'spanning_paired_2',     'spanning_paired_1',
                      'spanning_paired_1_t',   'spanning_paired_2_t',
                      'spanning_paired_1_s',   'spanning_paired_2_s',
                      'spanning_singleton_1',  'spanning_singleton_2',
                      'spanning_singleton_1_r','spanning_singleton_2_r']:
                
                if abs(pos1.get_dist(pos2,False)) < MIN_DISCO_INS_SIZE:
                    pos1 = None
                    pos2 = None
            
            if pos2 != None:
                self.insert_entry(pos1,pos2,rg,(read.cigarstring,parsed_SA_tag[2]),False)
        
        return (pos1, pos2)
    
    def reinsert_edges(self, edges):
        """Only works for Edges of which the _origin and _target Node
        still exists
        """
        for edge_t in edges:
            edge   = edge_t[0]
            edge_c = edge_t[1]
            
            node1 = edge._origin
            node2 = edge._target
            
            node1.set_edge(edge)
            node2.set_edge(edge_c)
    
    def remove_edge(self, edge):
        node1 = edge._origin
        node2 = edge._target
        
        node1.remove_edge(edge,"by-target")
        node2.remove_edge(edge,"by-origin")

        # Not necessary
        #del(edge)
        
        if node1.get_top_edge()[0] == None:
            self.remove_node(node1)

        if node2.get_top_edge()[0] == None:
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
    
    def pos_to_range(self,pos,_range,exclude_itself):
        for _pos in range(_range[0], _range[1], _range[2]):
            if exclude_itself and _pos == pos.pos:
                continue
            else:
                yield BreakPosition(pos._chr, _pos, pos.strand)
    
    def search_splice_edges_between(self,pos1,pos2, insert_size):
        # insert size is two directional
        target_range = (pos2.pos - insert_size, pos2.pos+ insert_size, None)
        
        for interval in self.idxtree[pos1._chr].search(pos1.pos - insert_size, pos1.pos + insert_size + 1):
            for key in interval[2].keys():
                node1 = interval[2][key]
                for edge in node1.edges.values():
                    if edge.target_in_range(target_range):
                        yield (edge)
    
    def search_edges_between(self,pos1, pos2, insert_size):
        """searches for reads inbetween two regions (e.g. break + ins. size):
        
        [     ]                   [     ]
         |  |                       | |
         |   -----------------------  |
          ----------------------------
        
        @todo: correct for splice junction width
        @todo use idxtree.search()
        """
        range1 = self.get_range(pos1,insert_size,False)
        range2 = self.get_range(pos2,insert_size,False)
        
        if range1[2] == -1:
            _min = range1[1]
            _max = range1[0]
            interval_iterator = self.idxtree[pos1._chr].search(_min - 1, _max + 1)
        else:
            _max = range1[1]
            _min = range1[0]
            interval_iterator = self.idxtree[pos1._chr].search(_min - 1, _max + 1)
        
        for interval in interval_iterator:
            if interval[0] != pos1.pos:
                if interval[2].has_key(pos1.strand):
                    node_i = interval[2][pos1.strand]
                    for edge in node_i.edges.values():
                        if edge.target_in_range(range2) and edge._target.position.strand == pos2.strand:
                            yield (edge)

    def print_chain(self):
        print "**************************************************************"
        for node in self:
            _str = str(node)
            if len(_str.strip()) > 0:
                print _str
        print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    def get_start_point(self):
        """Returns the top scoring edges in the chain
        
        By doing this iteratively, and removing the top edge in the mean time, you get some kind of unlogical ordering function
        @todo - this trick is silly: it uses a genometree and all recursive children of the nodes
                Add a set or list of vector that keeps track of all Edges ordered by a certain score on top of the tree?
        """
        maxscore = 0
        edge = None
        
        for node in self:
            _score, _edge, _node1, _node2 = node.get_top_edge()
            if _edge != None and _score > maxscore:
                maxscore = _score
                edge = _edge
        
        return edge
    
    def get_start_splice_junc(self):
        """Returns the chain with the higest number of counts
        """
        maxscore = (0,0)
        edge = None
        
        for node in self:
            for _edge in node:
                if _edge._origin.position._chr == _edge._target.position._chr:
                    score = (_edge.get_count('cigar_splice_junction'),_edge.get_count('cigar_soft_clip'))
                    if score[0] > maxscore[0] or (score[0] == maxscore[0] and score[1] > maxscore[1]):
                        edge = _edge
                        maxscore = score
        
        return edge
    
    def prune(self,insert_size):
        """Does some 'clever' tricks to merge edges together and reduce data points
        """
        logging.info("Finding and merging other edges in close proximity (insert size)")
        
        candidates = []
        k = 0
        candidate = self.get_start_point()
        self.print_chain()
        
        while candidate != None:
            #print "- pruning candidate"
            #if str(candidate) == "chr21:42860320/42860321(+)->chr21:42861520/42861521(-):(spanning_paired_1_t:4)":
            #    print "- pruning THE candidate!"
            self.prune_edge(candidate, insert_size)
            candidates.append((candidate, candidate.get_complement()))
            
            # do not remove if splice junc exists?
            self.remove_edge(candidate)
            
            candidate = self.get_start_point()
            k += 1
        
        #self.print_chain()
        logging.info("Pruned into "+str(k)+" candidate edge(s)")
        return candidates
    
    #@todo ADD RECURSION DEPTH LIMIT
    def prune_edge(self, edge, insert_size):
        ## These types of edges fall within the same space / inst. size
        ## they make the edges heavier
        ## @ todo double check if this is strand specific
        
        node1 = edge._origin
        node2 = edge._target
        
        edge_complement = edge.get_complement()
        
        i = 0
        for edge_m in self.search_edges_between(node1.position, node2.position, insert_size):
            i += 1
            edge_mc = edge_m.get_complement()
            
            edge.merge_edge(edge_m)
            edge_complement.merge_edge(edge_mc)
            
            self.remove_edge(edge_m)
            
            # complement is automatically removed after removing the fwd
            #self.remove_edge(edge_mc)
        
        return i
    
    def merge_splice_juncs(self,uncertainty):
        logging.info("Merging splice juncs")
        """@todo while loop?"""
        
        init = self.get_start_splice_junc()
        while init != None:
            init_c = init.get_complement()
            
            for xnode in self.idxtree[init._origin.position._chr].search(init._origin.position.pos - (uncertainty+1), init._origin.position.pos + (uncertainty + 1)):
                for node in xnode[2].values():
                    for junc in node:
                        if junc.get_count('cigar_splice_junction') > 0 and junc not in [init, init_c]:
                            d_origin = abs(junc._origin.position.pos - init._origin.position.pos)
                            d_target = abs(junc._target.position.pos - init._target.position.pos)
                            if d_origin <= uncertainty and d_target <= uncertainty:
                                init.merge_edge(junc)
                                init_c.merge_edge(junc.get_complement())
                                
                                self.remove_edge(junc)
                                
                                #self.print_chain()
            
            init = self.get_start_splice_junc()
        
        logging.debug("done")
    
    def rejoin_splice_juncs(self, thicker_edges, insert_size):
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
        #@todo use a separate genometree for this?
        
        k = 0
        logging.debug("Initiated")
        
        ## 01 collect all left and right nodes, indexed per chromosome
        left_nodes = {}
        right_nodes = {}
        
        for edge in thicker_edges:
            if not left_nodes.has_key(edge[0]._origin.position._chr):#lnodes
                left_nodes[edge[0]._origin.position._chr] = set()
            left_nodes[edge[0]._origin.position._chr].add(edge[0]._origin)
            
            if not right_nodes.has_key(edge[0]._target.position._chr):#rnodes
                right_nodes[edge[0]._target.position._chr] = set()
            right_nodes[edge[0]._target.position._chr].add(edge[0]._target)
        
        ## 02 look for all left nodes if there is any set (i < j) where
        ## i and j span a splice junction
        for _chr in left_nodes.keys():
            i = -1
            for node1 in left_nodes[_chr]:
                i += 1
                j = -1                
                
                for node2 in left_nodes[_chr]:
                    j += 1
                    
                    if j > i:# Avoid unnecessary comparisons
                        if node1.position.strand == node2.position.strand:
                            left_junc = (MAX_GENOMIC_DIST, None)
                            
                            for splice_junc in self.search_splice_edges_between(node1.position, node2.position, insert_size):
                                if splice_junc.get_count('cigar_splice_junction') > 0:#@todo and dist splice junction > ?insert_size?
                                    dist_origin1 = abs(splice_junc._origin.position.get_dist(node1.position, False))
                                    dist_origin2 = abs(splice_junc._target.position.get_dist(node2.position, False))
                                    sq_dist_origin = pow(dist_origin1, 2) + pow(dist_origin2, 2)
                                    
                                    if dist_origin1 < insert_size and dist_origin2 < insert_size and sq_dist_origin < left_junc[0]:
                                        left_junc = (sq_dist_origin, splice_junc)
                                
                            if left_junc[1] != None:
                                node1.splice_edges[node2] = left_junc
                                node2.splice_edges[node1] = left_junc
                                
                                k += 1

        for _chr in right_nodes.keys():
            i = -1
            for node1 in right_nodes[_chr]:
                i += 1
                j = -1
                
                for node2 in right_nodes[_chr]:
                    j += 1
                    
                    if j > i:# Avoid unnecessary comparisons
                        if node1.position.strand == node2.position.strand:
                            right_junc = (MAX_GENOMIC_DIST, None)
                            
                            for splice_junc in self.search_splice_edges_between(node1.position, node2.position, insert_size):
                                if splice_junc.get_count('cigar_splice_junction') > 0:#@todo and dist splice junction > ?insert_size?
                                    dist_target1 = abs(splice_junc._origin.position.get_dist(node1.position, False))
                                    dist_target2 = abs(splice_junc._target.position.get_dist(node2.position, False))
                                    sq_dist_target = pow(dist_target1, 2) +pow(dist_target2, 2)

                                    if dist_target1 < insert_size and dist_target2 < insert_size and sq_dist_target < right_junc[0]:
                                        right_junc = (sq_dist_target, splice_junc)
                            
                            if right_junc[1] != None:
                                node1.splice_edges[node2] = right_junc
                                node2.splice_edges[node1] = right_junc
                                
                                k += 1
        
        logging.info("Linked "+str(k)+" splice junction(s)")
        
        return thicker_edges
    
    def extract_subnetworks(self,thicker_edges):
        """ Here we add additional nodes an edge's current `left_node` 
or `right_node` using the `guilt-by-association` principle. Sometimes nodes
have edges to the same nodes of the already existing network,
            but lack splice junction(s) to those existing network. They're
            still connected to the network they might be exons that are not
            taken into account by the aligner (classical example: exon-0 in
            TMPRSS2). We have to be careful on the other hand, as multimap
            locations may influence this.
            
            Therefore we ideally need a function that adds nodes based on:
             - Genomic distance to all entries in either `left_nodes` or `right_nodes`
               * Currently implemented as RMSE on Node::get_dist()
             - The score of all edges going to the node
               * It also needs to correct for imbalance; 4 and 36 reads is less
                 likely than 20 and 20. Solution used multiplication.
             - Number of edges relative to the number of nodes.
               *In case of 3 nodes, it is more likely to have 3 edges (3/3) instead of 2
               (2/3). This weight is not yet implemented (01/08/16).


            Detection of possible candidates
            --------------------------------
            
            nodes:
            A     B         $ ... $          Y    Z
    edges:   ....................................... (9)
            ..................................      (2)
                  ................................. (100)
                  ............................      (20)
    splice:  -----
    
    
    start point: thickest edge
    

    left_nodes:  A, B
    right_nodes: Y, Z

    edges:
            A-Z (9)
            B-Z (100)
            
            A-Y (2)
            B-Y (20)
            
            For each left node i, compared to any other left node j,
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
        logging.info("Initiated -- *** rnodes not yet implemented ***")
        
        q = 0
        subnetworks = []
        while len(thicker_edges) > 0:
            start_point = thicker_edges[0][0]
            
            left_nodes = [start_point._origin]
            right_nodes = [start_point._target]
            
            ## The original nodes have been emptied, so the most important
            ## edge's are now separated.
            print left_nodes[0]
            left_nodes = start_point._origin.rfind_connected_sjuncs(left_nodes,450)
            #right_nodes = start_point._target.rfind_connected_sjuncs(right_nodes,450)
            print "\n\n"
            
            subedges = []
            
            # All edges between any of the left and right nodes are valid edges and have to be extracted
            # Find all direct edges joined by splice junctions
            for left_node in left_nodes:
                for right_node in right_nodes:
                    if left_node.edges.has_key(str(right_node.position)):
                        subedges.append( (left_node.edges[str(right_node.position)], right_node.edges[str(left_node.position)]) )
            
            ## Find all with one indirect step - these might be alternative junctions / exons
            i = -1
            for left_node_i in left_nodes:
                j = -1
                i += 1
                
                for left_node_j in left_nodes:
                    j += 1
                    
                    if i < j:
                        ## Find similar destinations
                        mutual_targets = list(set(left_node_i.edges.keys()).intersection(left_node_j.edges.keys()))
                        mutual_targets = [left_node_i.edges[mt]._target for mt in mutual_targets]
                        
                        for mt in mutual_targets:
                            if mt.is_connected_to((left_node_i, left_node_j), right_nodes):
                                # Add node
                                right_nodes.append(mt)
                                
                                for edge in mt.edges.keys():
                                    for l in left_nodes:
                                        if str(l.position) == edge:
                                            edge = mt.edges[edge]
                                            edge_c = edge.get_complement()
                                            
                                            # Make sure order is correct:
                                            subedges.append((edge_c,edge))
            
            ## Do inverse:
            
            #tmp = left_nodes
            #left_nodes = right_nodes
            #right_nodes = tmp

            ### @todo Redo code
            
            
            # remove all the links to the edges in each of the nodes
            for node in left_nodes:
                for edge_u in subedges:
                    for edge in edge_u:
                        key = str(edge._target.position)
                        if node.edges.has_key(key):
                            del(node.edges[key])

            # pop subedges from thicker edges and redo until thicker edges is empty
            popme = set()
            for edge in subedges:
                for edge2 in thicker_edges:
                    if edge[0] == edge2[0] or edge[1] == edge2[0]:
                        popme.add(edge2)
            
            for pop in popme:
                thicker_edges.remove(pop)
            
            subnetworks.append(Subnet(q,subedges))
            
            # @todo REMOVE THIS  - dev only:
            break
        
        logging.info("Extracted "+str(len(subnetworks))+" subnetwork(s)")
        return subnetworks


class Subnet(Chain):
    def __init__(self,_id,edges):
        self._id = _id
        self.edges = edges
        self.total_clips = 0
        self.total_score = 0
        self.discarded = []
        
        self.calc_clips()
        self.calc_scores()
        
        self.reorder_edges()
    
    def reorder_edges(self):
        idx = {}
        for edge in self.edges:
            key1 = edge[0].get_scores()
            key2 = str(edge[0]._origin.position)+"-"+str(edge[0]._target.position)
            
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
            nodes.add(edge[0]._origin)
            nodes.add(edge[0]._target)
        
        for node in nodes:
            clips += node.clips
        
        self.total_clips = clips
        return self.total_clips
    
    def calc_scores(self):
        score = 0
        for edge in self.edges:
            score += edge[0].get_scores()
        
        self.total_score = score
        return self.total_score
    
    def __str__(self):
        """Make tabular output
        
        print "chr-A\t"
        print "pos-A\t"
        print "direction-A\t"
        
        print "chr-B\t"
        print "pos-B\t"
        print "direction-B\t"
        
        print "discarded\t"
        
        print "score\t"
        print "soft+hardclips\t"
        print "n-split-reads\t"
        print "n-discordant-reads"
        
        print "n-edges\t"
        print "n-nodes-A\t"
        print "n-nodes-B\t"
        
        fh.write("entropy-bp-edge\t")
        fh.write("entropy-all-edges\t")
        """
        out = ""
        node_a = self.edges[0][0]._origin
        node_b = self.edges[0][0]._target
        
        out += str(node_a.position._chr)+"\t"
        out += str(node_a.position.pos)+"\t"
        out += strand_tt[node_a.position.strand]+"\t"
        
        out += str(node_b.position._chr)+"\t"
        out += str(node_b.position.pos)+"\t"
        out += strand_tt[node_b.position.strand]+"\t"
        
        out += ("valid" if self.discarded == [] else ','.join(self.discarded)) + "\t"
        
        out += str(self.total_score)+"\t"
        out += str(self.total_clips)+"\t"
        
        out += str(self.get_n_split_reads())+"\t"
        out += str(self.get_n_discordant_reads())+"\t"
        
        out += str(len(self.edges))+"\t"
        nodes_a, nodes_b = self.get_n_nodes()
        out += str(nodes_a)+"\t"
        out += str(nodes_b)+"\t"
        
        out += str(self.edges[0][0].get_entropy())+"\t"
        out += str(self.get_overall_entropy())+"\t"
        
        out += "&".join([str(edge[0]) for edge in self.edges])
        
        return out+"\n"
    
    def get_overall_entropy(self):
        frequency_table = merge_frequency_tables([edge[0].unique_alignments_idx for edge in self.edges])
        return entropy(frequency_table)
    
    def get_n_nodes(self):
        nodes_a = set()
        nodes_b = set()
        
        for edge in self.edges:
            nodes_a.add(edge[0]._origin)
            nodes_b.add(edge[0]._target)
        
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
            for edge in self.edges:
                n += edge[0].get_count(_type)
        
        return n
    
    def get_n_discordant_reads(self):
        return sum([edge[0].get_count("discordant_mates") for edge in self.edges])
    
    def find_distance(self, subnet_t):
        """
        Later on correct for that these are closer:
        |    |     ... ~ ...   |
           |                   |
        
        than these:
        |    |     ... ~ ...   |
      |                        |
        """
        if not isinstance(subnet_t, Subnet):
            raise Exception("subnet_t must be Subnet")
        
        ldist = MAX_GENOMIC_DIST
        rdist = MAX_GENOMIC_DIST

        for lnode in self.get_lnodes():
            for lnode_t in subnet_t.get_lnodes():
                dist = abs(lnode.position.get_dist(lnode_t.position, True))
                if dist == MAX_GENOMIC_DIST:
                    return (None, None, MAX_GENOMIC_DIST)
                
                if dist < ldist:
                    ldist = dist
        
        for rnode in self.get_rnodes():
            for rnode_t in subnet_t.get_rnodes():
                dist = abs(rnode.position.get_dist(rnode_t.position, True))
                if dist == MAX_GENOMIC_DIST:
                    return (None, None, MAX_GENOMIC_DIST)
                
                if dist < rdist:
                    rdist = dist
        
        return (ldist, rdist, math.sqrt(pow(ldist,2) + pow(rdist, 2)))
    
    def find_distances(self, subnet_t):
        """
            - Must compare Edges with Edges:
            
            self.edges:
            
         s1 [-------]
         s2 [---------------------]
            
            subnet_t.edges:
            
         t1    [----------------------------------------]
         t2     [-----------------]
         t3 [------------]
            
            This should return a vector containing min(2,3) = 2 distance entries:
            
            s1-t3   (dist s1.lnode <-> t3.lnode , dist s1.rnode <-> t3.rnode , RMSQD)
            s2-t2   (dist s2.lnode <-> t2.lnode , dist s2.rnode <-> t2.rnode , RMSQD)
            
            In order to be acceptable (this should be part of the merge_overlapping_subnets), each distance should:
             - RMSQD < math.fucn
             - min(ldist,rdist) < insert size
            Such that always one node is within the insert size of a node in the other subnet and the maximum distance can still be explained by exon distsances
        
        @todo
        """
        pass
    
    def get_lnodes(self):
        lnodes = set()
        
        for edge in self.edges:
            lnodes.add(edge[0]._origin)
        
        for lnode in lnodes:
            yield lnode
    
    def get_rnodes(self):
        rnodes = set()
        
        for edge in self.edges:
            rnodes.add(edge[0]._target)
        
        for rnode in rnodes:
            yield rnode
    
    def merge(self, subnet_m):
        for edge in subnet_m.edges:
            self.edges.append(edge)

        self.calc_clips()
        self.calc_scores()
        self.reorder_edges()


class IntronDecomposition:
    def __init__(self,alignment_file):
        self.pysam_fh = self.test_disco_alignment(alignment_file)
    
    def decompose(self):
        chain = Chain(self.pysam_fh)
        chain.insert_alignment()
        
        #@todo: thicker_edges = self.index_edges() and come up with class
        thicker_edges = chain.prune(PRUNE_INS_SIZE) # Makes edge thicker by lookin in the ins. size
        thicker_edges = chain.rejoin_splice_juncs(thicker_edges, PRUNE_INS_SIZE) # Merges edges by splice junctions and other junctions
        for te in thicker_edges:
            print te[0]
        print "-----------------------------------------"
        chain.reinsert_edges(thicker_edges)
        
        #subnets = extract_subnetworks_by_splice_junctions(thicker edges)
        subnets = chain.extract_subnetworks(thicker_edges)
        ##subnets = self.filter_subnets_on_identical_nodes(subnets)
        subnets = self.merge_overlapping_subnets(subnets, PRUNE_INS_SIZE)
        self.results = self.filter_subnets(subnets)# Filters based on three rules: entropy, score and background
        
        # If circos:
        #s = 1
        #for subnet in self.results:
        #    c = CircosController(str(s), subnet, "tmp/circos.conf","tmp/select-coordinates.conf", "tmp/circos-data.txt")
        #    c.draw_network("tmp/test.png","tmp/test.svg")
        #    s += 1
        
        return len(self.results)
    
    # @todo drop this into the BAMExtract class
    def test_disco_alignment(self,alignment_file):
        """Ensures by reading the BAM header whether the BAM file was
        indeed fixed using Dr. Disco
        """
        bam_fh = pysam.AlignmentFile(alignment_file, "rb")
        if bam_fh.header.has_key('PG'):
            for pg in bam_fh.header['PG']:
                if pg['ID'] == 'drdisco_fix_chimeric':
                    # create index if it does not exist
                    try:
                        bam_fh.fetch()
                    except:
                        logging.info('Indexing BAM file with pysam: '+bam_fh.filename)
                        pysam.index(bam_fh.filename)
                        bam_fh = pysam.AlignmentFile(bam_fh.filename)
                    
                    try:
                        bam_fh.fetch()
                    except:
                        raise Exception('Could not indexing BAM file: '+bam_fh.filename)
                    
                    return bam_fh
        
        raise Exception("Invalid STAR BAM File: has to be post processed with 'dr-disco fix-chimeric ...' first")
    
    def export(self, fh):
        fh.write(str(self))
    
    def __str__(self):
        out = "chr-A\t"
        out += "pos-A\t"
        out += "direction-A\t"
        
        out += "chr-B\t"
        out += "pos-B\t"
        out += "direction-B\t"
        
        out += "filter-status\t"
        
        out += "score\t"
        out += "soft+hardclips\t"
        out += "n-split-reads\t"
        out += "n-discordant-reads\t"
        
        out += "n-edges\t"
        out += "n-nodes-A\t"
        out += "n-nodes-B\t"
        
        out += "entropy-bp-edge\t"
        out += "entropy-all-edges\t"
        
        out += "data-structure\n"
        
        for subnet in self.results:
            out += str(subnet)
        
        return out
    
    #def filter_subnets_on_identical_nodes(self, subnets):
        ## This function is deprecated because it often happens that unknown exons are included
        ##
        #new_subnets = []
        #_id = 0
        ##1. filter based on nodes - if there are shared nodes, exclude sn
        #all_nodes = set()
        #for i in range(len(subnets)):
            #subnet = subnets[i]
            #rmme = False
            #score = 0
            #nodes = set()
            #for dp in subnet:
                #score += dp[0].get_scores()
                
                #nodes.add(dp[0]._origin)
                #nodes.add(dp[0]._target)
            
            #clips = 0
            #for n in nodes:
                #if n in all_nodes:
                    #rmme = True
                #else:
                    #all_nodes.add(n)
                
                #clips += n.clips
            
            #if rmme:
                #subnets[i] = None
            #else:
                #i += 1
                #sn = Subnet(_id, subnet)
                #sn.total_clips = clips
                #sn.total_score = score
                #new_subnets.append(sn)
        
        #return new_subnets

    def merge_overlapping_subnets(self, subnets, insert_size):
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
        n = len(subnets)
        
        k = 0
        for i in range(n):
            if subnets[i] != None:
                candidates = []
                for j in range(i+1,n):
                    if subnets[j] != None:
                        dist = subnets[i].find_distance(subnets[j])
                        #@todo use this one day
                        #subnets[i].find_distances(subnets[j])
                        
                        if dist[2] <= MAX_SUBNET_MERGE_DIST and min(dist[0],dist[1]) < insert_size:
                            candidates.append(subnets[j])
                            subnets[j] = None
            
                for sn_j in candidates:
                    subnets[i].merge(sn_j)
                    del(sn_j)
                    k += 1
        
        logging.info("Merged "+str(k)+" of the "+str(n)+" subnetwork(s)")
        
        return [sn for sn in subnets if sn != None]

    def filter_subnets(self, subnets):
        logging.debug("init")
        k = 0
        for subnet in subnets:
            """Total of 8 reads is minimum, of which 2 must be
            discordant and the entropy must be above 0.55"""
            
            entropy = subnet.get_overall_entropy()
            if entropy < MIN_SUBNET_ENTROPY:
                subnet.discarded.append("entropy="+str(entropy))
            
            n_disco = subnet.get_n_discordant_reads()
            n_disco_min = MIN_DISCO_PER_SUBNET_PER_NODE * sum(subnet.get_n_nodes())
            if n_disco < n_disco_min:
                subnet.discarded.append("n_discordant_reads="+str(n_disco)+"/"+str(n_disco_min))
            
            n_support = subnet.get_n_discordant_reads() + subnet.get_n_split_reads()
            n_support_min = (MIN_SUPPORTING_READS_PER_SUBNET_PER_NODE * sum(subnet.get_n_nodes()))
            if n_support < n_support_min:
                subnet.discarded.append("n_support="+str(n_support)+"/"+str(n_support_min))
            
            if len(subnet.discarded) > 0:
                k += 1
        
        logging.info("* Filtered "+str(k)+" of the "+str(len(subnets))+" subnetwork(s)")
        return subnets
