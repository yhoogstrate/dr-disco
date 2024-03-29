#!/usr/bin/env python
# *- coding: utf-8 -*-
# -- vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4
# https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

from drdisco.__init__ import MAX_ACCEPTABLE_INSERT_SIZE, MAX_ACCEPTABLE_ALIGNMENT_ERROR, MAX_GENOME_DISTANCE, MAX_SIZE_CIRCULAR_RNA
from drdisco.__init__ import STRAND_FORWARD, STRAND_REVERSE, STRAND_UNDETERMINED

import math
import operator

from tqdm import tqdm

import pysam
import HTSeq

import scipy.stats

from drdisco import log
from .CigarAlignment import cigar_to_cigartuple


"""[License: GNU General Public License v3 (GPLv3)]

    Dr. Disco: fusion gene detection in random hexamer RNA-seq data
    Copyright (C) 2017  Youri Hoogstrate

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    You can contact me via the github repository at the following url:
    <https://github.com/yhoogstrate/dr-disco>

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
    gmail dot com
"""

"""
Tries to find genomic and exon to exon break points within a discordant
RNA-Seq read alignment.
"""


def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    return sortedLst[index]
    # return lowest value if there are two center points
    # if (lstLen % 2):
    #     return sortedLst[index]
    # else:
    #     return min(sortedLst[index], sortedLst[index + 1])


strand_tt = {STRAND_FORWARD: '+', STRAND_REVERSE: '-', STRAND_UNDETERMINED: '?'}
x_onic_tt = {0: 'unknown', 1: 'exonic', 2: 'intronic'}


def entropy(frequency_table):
    n = sum(frequency_table.values())
    prob = [float(x) / n for x in frequency_table.values()]

    entropy = -sum([p * math.log(p) / math.log(2.0) for p in prob if p > 0])
    if entropy <= 0:
        return 0.0
    else:
        return entropy / (math.log(n) / math.log(2.0))


def stddev(vector):
    """ Calculates the variance of the breakpoints post merge - low values indicate consistent results"""
    n = len(vector)

    if n > 1:
        sq_dists = 0.0
        mean = 1.0 * sum(vector) / len(vector)
        for value in vector:
            dist = value - mean
            sq_dists += dist ** 2

        return math.sqrt(sq_dists / (n - 1))
    else:
        return 0.0


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
        for key in t:
            if key not in new_frequency_table:
                new_frequency_table[key] = 0

            new_frequency_table[key] += t[key]

    return new_frequency_table


def list_to_freq_table(as_list):
    table = {}
    for val in as_list:
        if val not in table:
            table[val] = 0
        table[val] += 1
    return table


def freq_table_to_list(ft):
    ls = []
    for val in ft:
        ls = ls + ft[val] * [val]
    return ls


def lin_regres_from_fq(ft):
    y = sorted(freq_table_to_list(ft))

    if(len(y) == 1):
        y = [y[0], y[0]]

    x = range(len(y))

    fit = scipy.stats.linregress(x, y)
    return fit


def bam_parse_alignment_offset(cigartuple, skip_N=False):
    pos = 0

    if skip_N:
        charset = [0, 2]
    else:
        charset = [0, 2, 3]

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

        if chunk[0] in charset:
            pos += chunk[1]

    return pos


def bam_parse_alignment_pos_using_cigar(sa_tag, skip_N=False):
    """SA tag looks like this:
        ['chr21', 42879876, '85S41M', '3', '-', '1']
        This function calculates the offset relative the the start position
    """

    pos = sa_tag[1]
    if sa_tag[4] == '+':
        cigartuple = cigar_to_cigartuple(sa_tag[2])
        pos += bam_parse_alignment_offset(cigartuple, skip_N=False)

    return pos


class BreakPosition:
    """
    Unambiguous data type that determines where a break point occurs.

    chr1:
    | 0 | 1 | 2 | 3 | ~~~

    chr2:
     ~~~ | 6 | 7 | 8 |

     string function of break points:
     chr1:3 /4 and chr2:5 /6
    """

    def __init__(self, _chr, position_0_based, strand):
        self._chr = _chr
        self.pos = position_0_based
        self.strand = strand
        self._hash = self._chr.replace("chr", "") + ("%0.2X" % self.pos).zfill(8) + strand_tt[self.strand]  # http://stackoverflow.com /questions /38430277 /python-class-hash-method-and-set

    def __lt__(self, other_bp):
        return self._hash < other_bp._hash

    def __str__(self):
        return str(self._chr) + ":" + str(self.pos) + "/" + str(self.pos + 1) + "(" + strand_tt[self.strand] + ")"

    def get_dist(self, other_bp, strand_specific):
        if not isinstance(other_bp, BreakPosition):  # pragma: no cover
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

    def __init__(self, position):
        self.position = position
        self.clips = 0
        self.edges = {}
        self.splice_edges = {STRAND_FORWARD: {}, STRAND_REVERSE: {}}

    def idx_sj_by_dist(self):
        if len(self.edges) == 0:  # Edges have been pruned - Node should actually be removed
            # @todo there is quite some performance to gain here i guess, as extract_by_sj iterates over all nodes that exist, even if they were used only before pruning
            # possibly remove all empty nodes directly after pruning
            self.splice_edges = {STRAND_FORWARD: [], STRAND_REVERSE: []}
        else:
            for strand in self.splice_edges:
                self.splice_edges[strand] = sorted([(value[0], key, value[1]) for (key, value) in self.splice_edges[strand].items()], key=operator.itemgetter(1, 2), reverse=False)

    def __lt__(self, other_node):
        return self.position < other_node.position

    def get_connected_splice_junctions(self):
        """Recursively finds all nodes that are connected by splice junctions

        self: the node of which it's splice junctions are traversed in order to find more nodes that may not be already in the blak
        self.splice_edges:
         -  [     ] splice edge: splice p1, splice p2
        """

        edges = {}

        def func(direction):
            idx = {self: MAX_ACCEPTABLE_INSERT_SIZE}
            todo = [self]  # Node objects
            i = 1

            while len(todo) > 0:
                next_iter = []

                for t in todo:
                    dist_left = idx[t]
                    for splice_edge in t.splice_edges[direction]:
                        cumulative_dist = dist_left - splice_edge[0]
                        if cumulative_dist >= 0:  # als cumulative_dist met een kleinere value al in idx zit, sla over
                            if splice_edge[1] in idx:
                                if cumulative_dist > idx[splice_edge[1]]:
                                    idx[splice_edge[1]] = cumulative_dist
                                    next_iter.append(splice_edge[1])

                                    edges[splice_edge[1]].add(splice_edge[2])

                                    # @todo kill entire tree / linked list?
                                    # in the example with e5.5 that is inbetween ex5 and ex6

                            else:
                                idx[splice_edge[1]] = cumulative_dist
                                next_iter.append(splice_edge[1])

                                edges[splice_edge[1]] = set([splice_edge[2]])
                        else:
                            break

                i += 1

                todo = next_iter
            return set(idx.keys())

        rnodes = func(STRAND_REVERSE)
        lnodes = func(STRAND_FORWARD)

        return rnodes.union(lnodes), edges

    def add_clip(self):
        self.clips += 1

    def get_edge_to_node(self, target_node):
        key = target_node.position._hash
        if key in self.edges:
            return self.edges[key]
        else:
            return None

    def insert_edge(self, edge):
        if edge.target == self:
            self.edges[edge.origin.position._hash] = edge
        else:
            self.edges[edge.target.position._hash] = edge

    def remove_edge(self, edge):
        try:
            del(self.edges[edge.origin.position._hash])
        except Exception:
            del(self.edges[edge.target.position._hash])

    def __str__(self):  # pragma: no cover
        out = str(self.position)

        a = 0
        for k in sorted(self.edges):
            edge = self.edges[k]
            filtered_edges = {JunctionTypeUtils.str(x): edge.types[x] for x in sorted(edge.types)}  # if x not in ['cigar_soft_clip','cigar_hard_clip']

            len_edges = len(filtered_edges)
            a += len_edges
            if len_edges > 0:
                out += "\n\t[" + str(id(edge)) + "][m:" + str(edge.n_matches) + ",mm:" + str(edge.n_mismatches) + "] " + str(edge.origin.position) + "->" + str(edge.target.position) + " " + str(filtered_edges)

        if a > 0:
            return out + "\n\t-> soft /hard clips: " + str(self.clips) + "\n"
        else:
            return ""


class JunctionTypes:  # Enum definition to avoid string operations
    silent_mate = 1
    cigar_soft_clip = 2
    cigar_hard_clip = 3
    cigar_deletion = 4
    cigar_splice_junction = 5
    discordant_mates = 6
    spanning_paired_1 = 7
    spanning_paired_1_r = 8
    spanning_paired_1_s = 9
    spanning_paired_1_t = 10
    spanning_paired_1_u = 11
    spanning_paired_2 = 12
    spanning_paired_2_r = 13
    spanning_paired_2_s = 14
    spanning_paired_2_t = 15
    spanning_paired_2_u = 16
    spanning_singleton_1 = 17
    spanning_singleton_1_r = 18
    spanning_singleton_2 = 19
    spanning_singleton_2_r = 20


class JunctionTypeUtils:
    tt = {
        'silent_mate': JunctionTypes.silent_mate,
        'cigar_splice_junction': JunctionTypes.cigar_splice_junction,
        'discordant_mates': JunctionTypes.discordant_mates,
        'spanning_paired_1': JunctionTypes.spanning_paired_1,
        'spanning_paired_2': JunctionTypes.spanning_paired_2,
        'spanning_paired_1_r': JunctionTypes.spanning_paired_1_r,
        'spanning_paired_2_r': JunctionTypes.spanning_paired_2_r,
        'spanning_paired_1_s': JunctionTypes.spanning_paired_1_s,
        'spanning_paired_2_s': JunctionTypes.spanning_paired_2_s,
        'spanning_paired_1_t': JunctionTypes.spanning_paired_1_t,
        'spanning_paired_2_t': JunctionTypes.spanning_paired_2_t,
        'spanning_paired_1_u': JunctionTypes.spanning_paired_1_u,
        'spanning_paired_2_u': JunctionTypes.spanning_paired_2_u,
        'spanning_singleton_1': JunctionTypes.spanning_singleton_1,
        'spanning_singleton_2': JunctionTypes.spanning_singleton_2,
        'spanning_singleton_1_r': JunctionTypes.spanning_singleton_1_r,
        'spanning_singleton_2_r': JunctionTypes.spanning_singleton_2_r}

    scoring_table = {
        JunctionTypes.discordant_mates: 1,
        JunctionTypes.spanning_paired_1: 3,
        JunctionTypes.spanning_paired_1_r: 3,
        JunctionTypes.spanning_paired_1_s: 3,
        JunctionTypes.spanning_paired_1_t: 3,
        JunctionTypes.spanning_paired_1_u: 3,
        JunctionTypes.spanning_paired_2: 3,
        JunctionTypes.spanning_paired_2_r: 3,
        JunctionTypes.spanning_paired_2_s: 3,
        JunctionTypes.spanning_paired_2_t: 3,
        JunctionTypes.spanning_paired_2_u: 3,
        JunctionTypes.spanning_singleton_1: 2,
        JunctionTypes.spanning_singleton_2: 2,
        JunctionTypes.spanning_singleton_1_r: 2,
        JunctionTypes.spanning_singleton_2_r: 2}

    complement_table = {
        JunctionTypes.discordant_mates: JunctionTypes.discordant_mates,
        JunctionTypes.spanning_paired_1: JunctionTypes.spanning_paired_2,
        JunctionTypes.spanning_paired_1_r: JunctionTypes.spanning_paired_2_r,
        JunctionTypes.spanning_paired_1_s: JunctionTypes.spanning_paired_2_s,
        JunctionTypes.spanning_paired_1_t: JunctionTypes.spanning_paired_2_t,
        JunctionTypes.spanning_paired_1_u: JunctionTypes.spanning_paired_2_u,
        JunctionTypes.spanning_paired_2: JunctionTypes.spanning_paired_1,
        JunctionTypes.spanning_paired_2_r: JunctionTypes.spanning_paired_1_r,
        JunctionTypes.spanning_paired_2_s: JunctionTypes.spanning_paired_1_s,
        JunctionTypes.spanning_paired_2_t: JunctionTypes.spanning_paired_1_t,
        JunctionTypes.spanning_paired_2_u: JunctionTypes.spanning_paired_1_u,
        JunctionTypes.spanning_singleton_1: JunctionTypes.spanning_singleton_2,
        JunctionTypes.spanning_singleton_2: JunctionTypes.spanning_singleton_1,
        JunctionTypes.spanning_singleton_1_r: JunctionTypes.spanning_singleton_2_r,
        JunctionTypes.spanning_singleton_2_r: JunctionTypes.spanning_singleton_1_r}

    @staticmethod
    def enum(enumstring):
        return JunctionTypeUtils.tt[enumstring]

    @staticmethod
    def str(enumcode):
        for key, value in JunctionTypeUtils.tt.items():
            if value == enumcode:
                return key

#    @staticmethod
#    def get_score(junction_type):
#        return JunctionTypeUtils.scoring_table[junction_type]

    @staticmethod
    def is_fusion_junction(junction_type):
        return junction_type >= 6


class Edge:
    """Connection between two genomic locations
     - the connection can by of different types

    Attributes:
        origin: this is a genomic location (BreakPosition) - donor / acceptor is not determined and term 'origin' and 'target' is meaningless
        target: this is a genomic location (BreakPosition)
    """

    def __init__(self, origin, target):
        if not isinstance(target, Node) or not isinstance(origin, Node):  # pragma: no cover
            raise Exception("origin and target must be a Node")

        self.origin = origin
        self.target = target
        self.types = {}
        self.unique_alignment_hashes = {}  # Used for determining entropy
        self.alignment_counter_positions = ({}, {})

        #  Used for determining entropy / rmsqd:
        self.unique_breakpoints = {True: ({}, {}),  # discordant mates: (origin , target)
                                   False: ({}, {})}  # other types: (origin, target)

        self.acceptor_donor = {}

        self.n_mismatches = 0
        self.n_matches = 0

        # These are the number of matches per piece of an edge:
        # [::::::::::>  -----   <:::]
        # will add p1=12, p2=5 given the alignment length of the chunks
        self.n_matches_p1 = []
        self.n_matches_p2 = []

    def get_entropy_of_alignments(self):
        """Entropy is an important metric / propery of an edge
        The assumption is that all reads should be 'different', i.e.
        have a different start position, different end position,
        different soft /hardclips etc.

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

        return entropy(self.unique_alignment_hashes)

    def get_distances_to_breakpoints(self, distances_of_discordant_mates):
        """
         In the following case all (non-disco) breakpoints are nearly identical:

          [==========> |
            [========> |
           [=========> |
          [==========> |
           [=========> |

         and here quite variable:

          [=========>   |
            [=======>   |
           [=========>  |
          [===========> |
           [==========> |

        Here we return the dist to the median of each of the values (0,0,0,0,0) for top example, (-1, -1, 0, 2, 2) for bottom
        """

        medians = (self.origin.position.pos, self.target.position.pos)
        dists = []
        for i in [0, 1]:
            for key in self.unique_breakpoints[distances_of_discordant_mates][i]:
                dists += self.unique_breakpoints[distances_of_discordant_mates][i][key] * [medians[i] - key]

        return dists

    def get_complement(self):
        try:
            return self.target.edges[self.origin.position._hash]
        except KeyError:  # todo write test for this
            raise KeyError("Could not find complement for edge:   " + str(self))

    def merge_edge(self, edge):
        """Merges (non splice) edges"""

        def update_pos(pos, i, is_discordant_mates, n):
            if pos not in self.unique_breakpoints[is_discordant_mates][i]:
                self.unique_breakpoints[is_discordant_mates][i][pos] = 0
            self.unique_breakpoints[is_discordant_mates][i][pos] += n

        def update_ctps(edge):
            for i in [0, 1]:
                for ctp in edge.alignment_counter_positions[i]:
                    n = edge.alignment_counter_positions[i][ctp]
                    if ctp not in self.alignment_counter_positions[i]:
                        self.alignment_counter_positions[i][ctp] = 0
                    self.alignment_counter_positions[i][ctp] += n

        def update_acceptor_donor(edge):
            for key in edge.acceptor_donor:
                if key not in self.acceptor_donor:
                    self.acceptor_donor[key] = edge.acceptor_donor[key]
                else:
                    self.acceptor_donor[key][0] += edge.acceptor_donor[key][0]
                    self.acceptor_donor[key][1] += edge.acceptor_donor[key][1]

        def update_matches_mismatches(edge):
            self.n_mismatches += edge.n_mismatches
            self.n_matches += edge.n_matches

        def update_n_matches_p1_p2(edge):
            for n_match_p1 in edge.n_matches_p1:
                self.n_matches_p1.append(n_match_p1)
            for n_match_p2 in edge.n_matches_p2:
                self.n_matches_p2.append(n_match_p2)

        for alignment_key in edge.unique_alignment_hashes:
            self.add_alignment_key(alignment_key)

        for key in [True, False]:
            for i in [0, 1]:
                for breakpoint in edge.unique_breakpoints[key][i]:
                    update_pos(breakpoint, i, key, edge.unique_breakpoints[key][i][breakpoint])

        for _type in edge.types:
            if _type != JunctionTypes.silent_mate:
                self.add_type(_type, edge.types[_type])

        update_ctps(edge)
        update_acceptor_donor(edge)
        update_matches_mismatches(edge)
        update_n_matches_p1_p2(edge)

    def add_type(self, _type, weight):
        if _type in self.types:
            self.types[_type] += weight
        else:
            self.types[_type] = weight

    def add_alignment_key(self, alignment_key):
        if alignment_key in self.unique_alignment_hashes:
            self.unique_alignment_hashes[alignment_key] += 1
        else:
            self.unique_alignment_hashes[alignment_key] = 1

    def add_ctps(self, ctps):
        """
        ctps = (dist1, dist2) << remaining lengths of remaining reads
        """
        for i in [0, 1]:
            ctp = ctps[i]
            if ctp in self.alignment_counter_positions[i]:
                self.alignment_counter_positions[i][ctp] += 1
            else:
                self.alignment_counter_positions[i][ctp] = 1

    def add_break_pos(self, pos1, pos2, is_discordant_mates):
        def update_pos(pos, i):
            if pos not in self.unique_breakpoints[is_discordant_mates][i]:
                self.unique_breakpoints[is_discordant_mates][i][pos] = 1
            else:
                self.unique_breakpoints[is_discordant_mates][i][pos] += 1

        update_pos(pos1.pos, 0)
        update_pos(pos2.pos, 1)

    def get_count(self, _type):
        if _type in self.types:
            return self.types[_type]
        else:
            return 0

    def get_score(self, _type):
        return self.get_count(_type) * JunctionTypeUtils.scoring_table[_type]

    def get_scores(self):
        """Based on this function, the start point is determined

        Currently it's implemented as sum(split reads)

        Later on spanning reads will be added, softclips will be added, etc.
        """
        return sum([self.get_score(_type) for _type in JunctionTypeUtils.scoring_table])

    def get_splice_score(self):  # pragma: no cover
        return (self.get_count(JunctionTypes.cigar_splice_junction), self.get_clips())

    def get_clips(self):  # pragma: no cover
        return self.origin.clips + self.target.clips

    def is_circular(self):
        return (self.origin.position._chr == self.target.position._chr) and \
               (self.origin.position.strand == STRAND_FORWARD) and \
               (self.target.position.strand == STRAND_REVERSE) and \
               (self.origin.position.get_dist(self.target.position, False) <= MAX_SIZE_CIRCULAR_RNA)

    def __str__(self):
        typestring = []

        for _t in sorted(self.types):
            typestring.append(str(JunctionTypeUtils.str(_t)) + ":" + str(self.types[_t]))

        out = str(self.origin.position) + "->" + str(self.target.position) + ":(" + ','.join(typestring) + ")"

        # spacer = " "*len(str(self.origin.position))

        # for _k in self.origin.splice_edges:
        #    out += "\n"+spacer+"=>"+str(_k.position)+":score="+str(self.origin.splice_edges[_k][1].get_splice_score())

        return out


class Graph:
    def __init__(self):
        self.idxtree = HTSeq.GenomicArrayOfSets("auto", stranded=False)

    def __iter__(self):
        for step in self.idxtree.steps():
            if step[1]:
                for strand in sorted(step[1]):
                    yield step[1][strand]

    def reindex_sj(self):
        for node in self:
            node.idx_sj_by_dist()

    def check_symmetry(self):  # debug function
        for edge in self.edge_idx:
            for key in edge.types:
                score = edge.get_score(key)
                cscore = edge.get_score(JunctionTypeUtils.complement_table[key])

                if score != cscore:
                    raise Exception("Read types %s and %s do not have an identical score\n%s" % (JunctionTypeUtils.str(key), JunctionTypeUtils.str(JunctionTypeUtils.complement_table[key]), str(edge)))

    def create_node(self, pos):
        position_accession = HTSeq.GenomicPosition(pos._chr, pos.pos)
        position = self.idxtree[position_accession]

        if position:
            if pos.strand not in position:
                position[pos.strand] = Node(pos)
        else:
            self.idxtree[position_accession] = {pos.strand: Node(pos)}

    def get_node_reference(self, pos):
        position = self.idxtree[HTSeq.GenomicPosition(pos._chr, pos.pos)]

        if position and pos.strand in position:
            return position[pos.strand]

        return None

    def insert_edge(self, pos1, pos2, _type, cigarstrs, ctps, acceptor, donor, n_mismatches, n_matches, n_match_p1, n_match_p2):
        """ - Checks if Node exists at pos1, otherwise creates one
            - Checks if Node exists at pos2, otherwise creates one
            - Checks if Edge exists between them, otherwise inserts it into the Nodes

         cigarstrs must be something like ("126M","126M") or ("25S50M2000N","25M50S")

         @todo what a MESS - C++ would be so neat here
        """
        self.create_node(pos1)
        self.create_node(pos2)

        node1 = self.get_node_reference(pos1)
        node2 = self.get_node_reference(pos2)

        if cigarstrs is not None:
            cigarstrs = pos1._hash + cigarstrs[0] + "|" + pos2._hash + cigarstrs[1]

        edge = node1.get_edge_to_node(node2)
        if edge is None:
            edge = Edge(node1, node2)
            node1.insert_edge(edge)
            node2.insert_edge(edge)

        edge.n_mismatches += n_mismatches
        edge.n_matches += n_matches

        edge.add_type(_type, 1)
        if node1 == edge.origin and _type != JunctionTypes.cigar_splice_junction:  # Avoid double insertion of all keys :) only do it if the positions don't get swapped
            edge.add_alignment_key(cigarstrs)
            edge.add_break_pos(pos1, pos2, (_type == JunctionTypes.discordant_mates))
            if ctps is not None:
                edge.add_ctps(ctps)

        if acceptor is not None:
            if str(acceptor) == str(edge.origin.position):
                pos_key = 'pos-A'
            elif str(acceptor) == str(edge.target.position):
                pos_key = 'pos-B'
            else:
                raise Exception("invalid acceptor position")

            if pos_key not in edge.acceptor_donor:
                edge.acceptor_donor[pos_key] = [1, 0]  # acceptor score, donor score
            else:
                edge.acceptor_donor[pos_key][0] += 1

        if donor is not None:
            if str(donor) == str(edge.origin.position):
                pos_key = 'pos-A'
            elif str(donor) == str(edge.target.position):
                pos_key = 'pos-B'
            else:
                raise Exception("invalid acceptor position")

            if pos_key not in edge.acceptor_donor:
                edge.acceptor_donor[pos_key] = [0, 1]  # acceptor score, donor score
            else:
                edge.acceptor_donor[pos_key][1] += 1

        if str(pos1) == str(edge.origin.position):
            edge.n_matches_p1.append(n_match_p1)
            edge.n_matches_p2.append(n_match_p2)
        else:
            edge.n_matches_p1.append(n_match_p2)
            edge.n_matches_p2.append(n_match_p1)

    def reinsert_edges(self, edges):
        """Only works for Edges of which the origin and target Node
        still exists
        """
        for edge in edges:
            node1 = edge.origin
            node2 = edge.target

            node1.insert_edge(edge)
            node2.insert_edge(edge)

        # self.print_chain()

    def remove_edge(self, edge):
        node1 = edge.origin
        node2 = edge.target

        node1.remove_edge(edge)
        node2.remove_edge(edge)

        del(edge)

    def remove_node(self, node):
        self.idxtree[HTSeq.GenomicPosition(node.position._chr, node.position.pos)] = set()
        del(node)

    def print_chain(self):  # pragma: no cover
        print("**************************************************************")
        for node in self:
            _str = str(node)
            if len(_str.strip()) > 0:
                print(_str)
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    def generate_edge_idx(self):
        edges = set()
        edges_tuple = []
        order = 0

        for node in self:
            for edge in node.edges.values():
                if edge not in edges:
                    score = edge.get_scores()
                    if score > 0:
                        edges.add(edge)
                        edges_tuple.append((edge, score, order))
                        order -= 1

        del(edges, order)

        self.edge_idx = [edge[0] for edge in sorted(edges_tuple, key=operator.itemgetter(1, 2), reverse=False)]

    def prune(self, pbarx):
        """Does some 'clever' tricks to merge edges together and reduce data points
        """
        # self.print_chain()

        candidates = []

        previous = len(self.edge_idx)
        while self.edge_idx:
            candidate = self.edge_idx.pop()

            self.prune_edge(candidate)
            candidates.append(candidate)

            self.remove_edge(candidate)  # do not remove if splice junc exists?

            pbarx.update(previous - len(self.edge_idx))
            previous = len(self.edge_idx)

        return candidates

    def prune_edge(self, edge):
        #  @todo double check if this is strand specific

        to_merge = []

        for edge_m in self.search_edges_between(edge):
            d1 = edge.origin.position.get_dist(edge_m.origin.position, True)
            d2 = edge.target.position.get_dist(edge_m.target.position, True)
            d = abs(d1) + abs(d2)

            if d <= MAX_ACCEPTABLE_INSERT_SIZE:
                to_merge.append(edge_m)

        for edge_m in to_merge:
            edge.merge_edge(edge_m)

            self.remove_edge(edge_m)
            self.edge_idx.remove(edge_m)


    def search_edges_between(self, edge_to_prune, skip_duplicates = True):
        """searches for other junctions in-between edge+insert size:
        
        if there are small junctions with both ends located within the search area, it can be returned twice
        """
        def pos_to_range(pos):
            if pos.strand == STRAND_REVERSE:
                #  return (pos.pos - MAX_ACCEPTABLE_INSERT_SIZE) - 1, pos.pos + MAX_ACCEPTABLE_ALIGNMENT_ERROR
                return (pos.pos - MAX_ACCEPTABLE_INSERT_SIZE) - 1, pos.pos + MAX_ACCEPTABLE_INSERT_SIZE
            else:
                #  return pos.pos - MAX_ACCEPTABLE_ALIGNMENT_ERROR, (pos.pos + MAX_ACCEPTABLE_INSERT_SIZE) + 1
                return pos.pos - MAX_ACCEPTABLE_ALIGNMENT_ERROR, (pos.pos + MAX_ACCEPTABLE_INSERT_SIZE) + 1

        pos1, pos2 = edge_to_prune.origin.position, edge_to_prune.target.position

        pos1_min, pos1_max = pos_to_range(pos1)
        pos2_min, pos2_max = pos_to_range(pos2)
        
        duplicates = set([])

        for step in self.idxtree[HTSeq.GenomicInterval(pos1._chr, max(0, pos1_min), pos1_max + 1)].steps():
            if step[1]:
                position = step[1]
                if position and pos1.strand in position:
                    node_i = position[pos1.strand]
                    for edge in node_i.edges.values():
                        if edge != edge_to_prune and edge.target.position.strand == pos2.strand and edge.target.position.pos >= pos2_min and edge.target.position.pos <= pos2_max:
                            
                            if skip_duplicates:
                                if edge not in duplicates:
                                    yield edge
                                duplicates.add(edge)
                            
                            else:
                                yield edge

    def rejoin_splice_juncs(self, splice_junctions_all, chrs):
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
        k = 0

        def search(pos1):
            # @todo make this member of Graph and use splice_junctions.search_..._..(pos)
            for step in self.idxtree[HTSeq.GenomicInterval(pos1._chr, max(0, pos1.pos - MAX_ACCEPTABLE_INSERT_SIZE), max(1, pos1.pos + MAX_ACCEPTABLE_INSERT_SIZE + 1))].steps():
                if step[1]:
                    for strand in step[1]:
                        yield step[1][strand]

        splice_edges_had = set()

        for _chr in chrs:
            if _chr in splice_junctions_all:
                splice_junctions = splice_junctions_all[_chr]
                for node in splice_junctions:
                    for splice_junction in node.edges.values():
                        if splice_junction not in splice_edges_had:
                            for lnode in search(splice_junction.origin.position):
                                for rnode in search(splice_junction.target.position):
                                    if lnode.position.strand == rnode.position.strand and lnode != rnode:
                                        d1 = abs(splice_junction.origin.position.get_dist(lnode.position, False))
                                        d2 = abs(splice_junction.target.position.get_dist(rnode.position, False))
                                        dist = d1 + d2
                                        if dist <= MAX_ACCEPTABLE_INSERT_SIZE:
                                            if lnode.position < rnode.position and rnode in lnode.splice_edges[STRAND_FORWARD]:
                                                insert = dist < lnode.splice_edges[STRAND_FORWARD][rnode][0]
                                            elif rnode.position < lnode.position and rnode in lnode.splice_edges[STRAND_REVERSE]:
                                                insert = dist < lnode.splice_edges[STRAND_REVERSE][rnode][0]
                                            else:
                                                insert = True

                                            if insert:
                                                if rnode.position < lnode.position:  # inversed
                                                    rnode.splice_edges[STRAND_FORWARD][lnode] = (dist, splice_junction)
                                                    lnode.splice_edges[STRAND_REVERSE][rnode] = (dist, splice_junction)
                                                else:
                                                    lnode.splice_edges[STRAND_FORWARD][rnode] = (dist, splice_junction)
                                                    rnode.splice_edges[STRAND_REVERSE][lnode] = (dist, splice_junction)

                                                k += 1

                            splice_edges_had.add(splice_junction)
                            splice_edges_had.add(splice_junction.get_complement())

        return k

    def extract_subnetworks_by_splice_junctions(self, thicker_edges, MIN_SCORE_FOR_EXTRACTING_SUBGRAPHS):
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
        thicker_edges.reverse()
        q = 0
        subnetworks = []
        while thicker_edges:
            start_point = thicker_edges.pop()
            left_splice_junctions = set()
            right_splice_junctions = set()

            if start_point.get_scores() >= MIN_SCORE_FOR_EXTRACTING_SUBGRAPHS:
                left_nodes, left_splice_junctions_ds = start_point.origin.get_connected_splice_junctions()
                right_nodes, right_splice_junctions_ds = start_point.target.get_connected_splice_junctions()

                subedges = []  # [start_point] if code below is commented out

                # All edges between any of the left and right nodes are valid edges and have to be extracted
                # Find all direct edges joined by splice junctions
                for left_node in left_nodes:
                    for right_node in right_nodes:
                        if right_node.position._hash in left_node.edges:
                            subedge = left_node.edges[right_node.position._hash]
                            if subedge.target.position._hash == right_node.position._hash:
                                subedges.append(subedge)

                                if left_node != start_point.origin:  # must be merged by a splice junction
                                    left_splice_junctions = left_splice_junctions.union(left_splice_junctions_ds[left_node])

                                if right_node != start_point.target:
                                    right_splice_junctions = right_splice_junctions.union(right_splice_junctions_ds[right_node])

                                if subedge != start_point:
                                    self.remove_edge(subedge)

                                    thicker_edges.remove(subedge)  # goes wrong
            else:
                subedges = [start_point]

            # del(left_splice_junctions_ds,right_splice_junctions_ds)
            subnetworks.append(SubGraph(q, subedges, left_splice_junctions, right_splice_junctions))
            self.remove_edge(start_point)

        return subnetworks


class SubGraph():
    def __init__(self, _id, edges, left_splice_junctions, right_splice_junctions):
        self._id = _id
        self.edges = edges
        self.left_splice_junctions = left_splice_junctions
        self.right_splice_junctions = right_splice_junctions
        self.total_clips = 0
        self.total_score = 0

        self.xonic = 0

        self.calc_clips()
        self.calc_scores()
        self.calc_acceptor_donor()

        self.reorder_edges()

    def reorder_edges(self):
        idx = {}
        for edge in self.edges:
            key1 = edge.get_scores()
            key2 = edge.origin.position._hash + "-" + edge.target.position._hash

            if key1 not in idx:
                idx[key1] = {}

            idx[key1][key2] = edge

        ordered = []
        for key1 in sorted(idx, reverse=True):
            for key2 in sorted(idx[key1]):
                ordered.append(idx[key1][key2])

        self.edges = ordered

    def calc_clips(self):
        clips = 0

        nodes = set()
        for edge in self.edges:
            nodes.add(edge.origin)
            nodes.add(edge.target)

        for node in nodes:
            clips += node.clips

        self.total_clips = clips
        return self.total_clips

    def calc_m_mm(self):
        m = 0
        mm = 0

        for edge in self.edges:
            m += edge.n_matches
            mm += edge.n_mismatches

        return (m, mm)

    def calc_scores(self):
        self.total_score = 0
        for edge in self.edges:
            self.total_score += edge.get_scores()

        return self.total_score

    def calc_acceptor_donor(self):
        self.posA_acceptor = 0
        self.posA_donor = 0
        self.posB_acceptor = 0
        self.posB_donor = 0

        for edge in self.edges:
            for key in edge.acceptor_donor:
                if key == 'pos-A':
                    self.posA_acceptor += edge.acceptor_donor[key][0]
                    self.posA_donor += edge.acceptor_donor[key][1]
                elif key == 'pos-B':
                    self.posB_acceptor += edge.acceptor_donor[key][0]
                    self.posB_donor += edge.acceptor_donor[key][1]
                else:
                    raise Exception("error")

    def classify_intronic_exonic(self, splice_junctions_all, chrs):
        for pos in [self.edges[0].origin.position, self.edges[0].target.position]:
            for _chr in chrs:
                if _chr in splice_junctions_all:
                    splice_junctions = splice_junctions_all[_chr]
                    for step in splice_junctions.idxtree[HTSeq.GenomicInterval(pos._chr, max(0, pos.pos - MAX_ACCEPTABLE_ALIGNMENT_ERROR), max(1, pos.pos + MAX_ACCEPTABLE_ALIGNMENT_ERROR + 1))].steps():
                        if step[1]:
                            for strand in step[1]:
                                self.xonic = 1
                                return self.xonic
        self.xonic = 2

    def get_breakpoint_stddev(self):
        """ Calculates the variance of the breakpoints post merge - low values indicate consistent results"""
        sq_dists = 0.0
        n = 0
        for edge in self.edges:
            for dist in edge.get_distances_to_breakpoints(False):
                sq_dists += dist ** 2
                # dist * dist ? math.pow(dist,2)?
                n += 1
        if n > 0:
            return math.sqrt(sq_dists / float(n))
        else:
            return 0.0

    def get_breakpoint_disco_entropy(self):
        """ Calculates the entropy of the breakpoints post merge - low values indicate consistent results"""
        as_list = []

        for edge in self.edges:
            as_list += edge.get_distances_to_breakpoints(True)

        return entropy(list_to_freq_table(as_list))

    def get_unique_breakpoint_stats(self):
        set1 = merge_frequency_tables([edge.unique_breakpoints[False][0] for edge in self.edges if len(edge.unique_breakpoints[False][0]) > 0])
        set2 = merge_frequency_tables([edge.unique_breakpoints[False][1] for edge in self.edges if len(edge.unique_breakpoints[False][1]) > 0])
        return stddev(freq_table_to_list(set1)), entropy(set1), stddev(freq_table_to_list(set2)), entropy(set2)

    def __str__(self):
        """Makes tabular output"""
        node_a, node_b = self.edges[0].origin, self.edges[0].target
        nodes_a, nodes_b = self.get_n_nodes()
        dist = node_a.position.get_dist(node_b.position, False)
        if dist == MAX_GENOME_DISTANCE:
            dist = 'inf'

        lregA = lin_regres_from_fq(self.edges[0].alignment_counter_positions[0])
        lregB = lin_regres_from_fq(self.edges[0].alignment_counter_positions[1])
        m_mm = self.calc_m_mm()

        s1, s2, s3, s4 = self.get_unique_breakpoint_stats()

        return (
            "%s\t%i\t%s\t"
            "%i\t%i\t"
            "%s\t%i\t%s\t"
            "%i\t%i\t"
            "%s\tunclassified\t%s\t%s\t"
            "%i\t%i\t%i\t%i\t"
            "%i\t%i\t"
            "%i\t%i\t%i\t"
            "%i\t%i\t"
            "%.4f\t%.4f\t"
            "%.4f\t%.4f\t"
            "%.4f\t%.4f\t%.4f\t%.4f\t"  # num_unique L, num_unique R, entropy L, entropy L, self.unique_breakpoints = {True: ({}, {}),  # discordant mates: (origin , target)  {False:}other types: (origin, target)
            "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t"  # lreg-A
            "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t"  # lreg-B
            "%.4f\t%.4f\t%.4f\t"
            "%i\t%i\t%i\t%i\t"
            "%s\n" % (node_a.position._chr, node_a.position.pos, strand_tt[self.edges[0].origin.position.strand],  # Pos-A
                      self.posA_acceptor, self.posA_donor,
                      node_b.position._chr, node_b.position.pos, strand_tt[self.edges[0].target.position.strand],  # Pos-B
                      self.posB_acceptor, self.posB_donor,
                      dist, ("circular" if self.edges[0].is_circular() else "linear"), x_onic_tt[self.xonic],  # Classification status
                      self.total_score / 2, self.total_clips, self.get_n_split_reads() / 2, self.get_n_discordant_reads() / 2,  # Evidence stats
                      m_mm[0], m_mm[1],
                      len(self.edges), nodes_a, nodes_b,  # Edges and nodes stats
                      len(self.left_splice_junctions), len(self.right_splice_junctions),
                      self.edges[0].get_entropy_of_alignments(), self.get_overall_entropy(),  # Entropy stats
                      self.get_breakpoint_stddev(), self.get_breakpoint_disco_entropy(),
                      s1, s2, s3, s4,
                      lregA[0], lregA[1], lregA[2], lregA[3], lregA[4],
                      lregB[0], lregB[1], lregB[2], lregB[3], lregB[4],
                      1.0 * self.get_n_discordant_reads() / max(self.get_n_split_reads(), 0.1),
                      1.0 * self.total_clips / self.total_score,
                      1.0 * (nodes_a + nodes_b) / len(self.edges),
                      self.get_n_matches_p1_median(), self.get_n_matches_p2_median(),
                      self.get_n_matches_p1_max(), self.get_n_matches_p2_max(),
                      "&".join([str(edge) for edge in self.edges])))  # Data structure

    def get_overall_entropy(self):
        frequency_table = merge_frequency_tables([edge.unique_alignment_hashes for edge in self.edges])
        return entropy(frequency_table)

    def get_n_nodes(self):
        nodes_a, nodes_b = set(), set()

        for edge in self.edges:
            nodes_a.add(edge.origin)
            nodes_b.add(edge.target)

        return len(nodes_a), len(nodes_b)

    def get_n_split_reads(self):
        n = 0

        for _type in [JunctionTypes.spanning_paired_1, JunctionTypes.spanning_paired_2,
                      JunctionTypes.spanning_paired_1_r, JunctionTypes.spanning_paired_2_r,
                      JunctionTypes.spanning_paired_1_s, JunctionTypes.spanning_paired_2_s,
                      JunctionTypes.spanning_paired_1_t, JunctionTypes.spanning_paired_2_t,
                      JunctionTypes.spanning_singleton_1, JunctionTypes.spanning_singleton_2,
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

        if len(l_nodes_min) > len(l_nodes_max):
            l_nodes_min, l_nodes_max = l_nodes_max, l_nodes_min

        if len(r_nodes_min) > len(r_nodes_max):
            r_nodes_min, r_nodes_max = r_nodes_max, r_nodes_min

        l_dists = []
        r_dists = []

        for l_node in l_nodes_min:  # makes it symmetrical
            dist = MAX_GENOME_DISTANCE

            for l_node_t in l_nodes_max:
                dist = min(abs(l_node.position.get_dist(l_node_t.position, True)), dist)

            if dist == MAX_GENOME_DISTANCE:
                return False, False
            else:
                l_dists.append(dist)

        for r_node in r_nodes_min:  # makes it symmetrical
            dist = MAX_GENOME_DISTANCE

            for r_node_t in r_nodes_max:
                dist = min(abs(r_node.position.get_dist(r_node_t.position, True)), dist)

            if dist == MAX_GENOME_DISTANCE:
                return False, False
            else:
                r_dists.append(dist)

        return l_dists, r_dists

    def get_lnodes(self):
        lnodes = set()

        for edge in self.edges:
            lnodes.add(edge.origin)

        for lnode in lnodes:
            yield lnode

    def get_rnodes(self):
        rnodes = set()

        for edge in self.edges:
            rnodes.add(edge.target)

        for rnode in rnodes:
            yield rnode

    def get_n_matches_p1_median(self):
        a = []
        for edge in self.edges:
            for n_match_p1 in edge.n_matches_p1:
                a.append(n_match_p1)

        return median(a)

    def get_n_matches_p2_median(self):
        a = []
        for edge in self.edges:
            for n_match_p2 in edge.n_matches_p2:
                a.append(n_match_p2)

        return median(a)

    def get_n_matches_p1_max(self):
        a = []
        for edge in self.edges:
            for n_match_p1 in edge.n_matches_p1:
                a.append(n_match_p1)

        return max(a)

    def get_n_matches_p2_max(self):
        a = []
        for edge in self.edges:
            for n_match_p2 in edge.n_matches_p2:
                a.append(n_match_p2)

        return max(a)

    def merge(self, subnet_m):
        for edge in subnet_m.edges:
            self.edges.append(edge)

        self.left_splice_junctions = self.left_splice_junctions.union(subnet_m.left_splice_junctions)
        self.right_splice_junctions = self.right_splice_junctions.union(subnet_m.right_splice_junctions)

        self.calc_clips()
        self.calc_scores()
        self.reorder_edges()


class BAMExtract(object):
    def __init__(self, bam_file, require_fixed_bam_file):
        self.pysam_fh = self.test_disco_alignment(bam_file, require_fixed_bam_file)
        self.overlapping_reads = 0
        self.crosshybridization_artifacts = 0

    @staticmethod
    def test_disco_alignment(alignment_file, require_fixed_bam_file):
        """Ensures by reading the BAM header whether the BAM file was
        indeed fixed using Dr. Disco
        """
        bam_fh = pysam.AlignmentFile(alignment_file, "rb")

        if require_fixed_bam_file:
            proper_tag = False
            if 'PG' in bam_fh.header:
                for pg in bam_fh.header['PG']:
                    if pg['ID'] == 'drdisco_fix_chimeric':
                        proper_tag = True

            if not proper_tag:
                bam_fh.close()
                raise Exception("Invalid STAR BAM File: has to be post processed with 'dr-disco fix-chimeric ...' first")

        try:  # pragma: no cover
            bam_fh.fetch()
        except Exception:  # pragma: no cover
            fname = bam_fh.filename.decode("utf-8")
            log.info('Indexing BAM file with pysam: ' + fname)  # create index if it does not exist
            bam_fh.close()

            pysam.index(fname)
            bam_fh = pysam.AlignmentFile(bam_fh.filename)

        try:
            bam_fh.fetch()
        except Exception:
            raise Exception('Could not indexing BAM file: ' + bam_fh.filename)

        return bam_fh

    def extract_junctions(self, fusion_junctions, splice_junctions):
        def read_to_junction(read, rg, parsed_SA_tags, specific_type=None):
            parsed_SA_tag = parsed_SA_tags[0]
            parsed_cigar_tuple = cigar_to_cigartuple(parsed_SA_tag[2])
            pos1, pos2, acceptor, donor = None, None, None, None

            if rg == JunctionTypes.discordant_mates:
                # How to distinguish between:
                #
                # Type 1:
                # ===1===>|   |<===2===

                # Type 2:
                # ===2===>|   |<===1===

                # Hypothetical:
                # |<===1===    |<===2==

                if read.is_reverse:
                    pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                         read.reference_start,
                                         STRAND_FORWARD)
                else:
                    pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                         read.reference_start + bam_parse_alignment_offset(read.cigar),
                                         STRAND_REVERSE)

                if parsed_SA_tag[4] == "-":
                    pos2 = BreakPosition(parsed_SA_tag[0],
                                         parsed_SA_tag[1],
                                         STRAND_FORWARD)
                else:
                    pos2 = BreakPosition(parsed_SA_tag[0],
                                         bam_parse_alignment_pos_using_cigar(parsed_SA_tag),
                                         STRAND_REVERSE)

                # more or less tested
                if read.is_read1:
                    donor, acceptor = pos2, pos1
                else:
                    donor, acceptor = pos1, pos2

            elif rg == JunctionTypes.spanning_singleton_1:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + bam_parse_alignment_offset(read.cigar),
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1],
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

                # more or less tested
                donor, acceptor = pos2, pos1

            elif rg == JunctionTypes.spanning_singleton_2:
                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start,
                                     read.is_reverse)

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + bam_parse_alignment_offset(parsed_cigar_tuple),
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

                donor, acceptor = pos1, pos2

            elif rg == JunctionTypes.spanning_singleton_1_r:
                # test 24 covers these - positions are 100% correct, strands not sure ...
                if not read.is_reverse:
                    pos1_offset = 0
                else:
                    pos1_offset = bam_parse_alignment_offset(read.cigar)

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)

                if parsed_SA_tag[4] == "-":
                    pos2_offset = 0
                else:
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

                donor, acceptor = pos1, pos2

            elif rg in [JunctionTypes.spanning_singleton_2_r]:
                if read.is_reverse:
                    pos1_offset = 0
                else:
                    pos1_offset = bam_parse_alignment_offset(read.cigar)

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     read.is_reverse)

                if parsed_SA_tag[4] == "-":
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)
                else:
                    pos2_offset = 0

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

                donor, acceptor = pos2, pos1

            elif rg in [JunctionTypes.spanning_paired_1]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = bam_parse_alignment_offset(read.cigar)
                else:
                    pos1_offset = 0

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)

                if parsed_SA_tag[4] == '+':   # Tested in TERG s55 double inversion
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)
                else:
                    pos2_offset = 0

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

                donor, acceptor = pos2, pos1

            elif rg in [JunctionTypes.spanning_paired_2]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = 0
                else:
                    pos1_offset = bam_parse_alignment_offset(read.cigar)

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     STRAND_REVERSE if not read.is_reverse else STRAND_FORWARD)

                if parsed_SA_tag[4] == '-':  # Tested in TERG s55 double inversion
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)
                else:
                    pos2_offset = 0

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

                donor, acceptor = pos1, pos2

            elif rg in [JunctionTypes.spanning_paired_1_r]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = bam_parse_alignment_offset(read.cigar)
                else:
                    pos1_offset = 0

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)

                if parsed_SA_tag[4] == '+':   # Tested in TERG s55 double inversion
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)
                else:
                    pos2_offset = 0

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

                donor, acceptor = pos2, pos1

            elif rg in [JunctionTypes.spanning_paired_2_r]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = 0
                else:
                    pos1_offset = bam_parse_alignment_offset(read.cigar)

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     read.is_reverse)

                if parsed_SA_tag[4] == '-':  # Tested in TERG s55 double inversion
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)
                else:
                    pos2_offset = 0

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

                donor, acceptor = pos1, pos2

            elif rg in [JunctionTypes.spanning_paired_1_s]:
                # Very clear example in S054 @ chr21:40, 064, 610-40, 064, 831  and implemented in test case 14
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = 0
                else:
                    pos1_offset = bam_parse_alignment_offset(read.cigar)

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     STRAND_REVERSE if not read.is_reverse else STRAND_FORWARD)

                if parsed_SA_tag[4] == "-":
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)
                else:
                    pos2_offset = 0

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_REVERSE if parsed_SA_tag[4] == "-" else STRAND_FORWARD)

                donor, acceptor = pos2, pos1

            elif rg in [JunctionTypes.spanning_paired_2_s]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = bam_parse_alignment_offset(read.cigar)
                else:
                    pos1_offset = 0

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     STRAND_REVERSE if read.is_reverse else STRAND_FORWARD)

                if parsed_SA_tag[4] == "-":
                    pos2_offset = 0
                else:
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] != "+" else STRAND_REVERSE)

                donor, acceptor = pos1, pos2

            elif rg in [JunctionTypes.spanning_paired_1_t]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = bam_parse_alignment_offset(read.cigar)
                else:
                    pos1_offset = 0

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     not read.is_reverse)

                if parsed_SA_tag[4] == "-":  # Tested in TERG s55 double inversion
                    pos2_offset = 0
                else:
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "-" else STRAND_REVERSE)

                donor, acceptor = pos1, pos2

            elif rg in [JunctionTypes.spanning_paired_2_t]:
                if read.is_reverse:  # Tested in TERG s55 double inversion
                    pos1_offset = 0
                else:
                    pos1_offset = bam_parse_alignment_offset(read.cigar)

                pos1 = BreakPosition(self.pysam_fh.get_reference_name(read.reference_id),
                                     read.reference_start + pos1_offset,
                                     read.is_reverse)

                if parsed_SA_tag[4] == "+":  # Tested in TERG s55 double inversion
                    pos2_offset = 0
                else:
                    pos2_offset = bam_parse_alignment_offset(parsed_cigar_tuple)

                pos2 = BreakPosition(parsed_SA_tag[0],
                                     parsed_SA_tag[1] + pos2_offset,
                                     STRAND_FORWARD if parsed_SA_tag[4] == "+" else STRAND_REVERSE)

                donor, acceptor = pos2, pos1

            elif rg != JunctionTypes.silent_mate:  # pragma: no cover
                raise Exception("Fatal Error, RG: %s" % JunctionTypeUtils.str(rg))

            def is_hybridization_artifact(read, pos1, pos2, sa_tags):
                """
                Example 01:

                silent       hit
                 mate       index=1

                [=====>      <==:
                 [==>

                  hit
                index=2

                assumption: HI:2 is messed up AND INVERTED (relative to HI:1)
                """
                #  silent mate exists & intrachromosomal       & (-- or ++)                   &   strand is actually determined
                if len(sa_tags) > 1 and pos1._chr == pos2._chr and pos1.strand == pos2.strand and pos1.strand is not None:
                    silent_mate = sa_tags[1]
                    read_strand = '-' if read.is_reverse else '+'

                    if silent_mate[0] == pos1._chr:
                        offset1 = sa_tags[1][1]
                        for chunk1 in cigar_to_cigartuple(sa_tags[1][2]):
                            if chunk1[0] == 0:
                                range1 = offset1, chunk1[1] + offset1

                                # 01 - is silent mate (sa_tags[1]) overlappnig the other split read (sa_tags[0])
                                if sa_tags[1][4] == sa_tags[0][4]:
                                    offset2 = sa_tags[0][1]
                                    for chunk2 in cigar_to_cigartuple(sa_tags[0][2]):
                                        if chunk2[0] == 0:
                                            range2 = offset2, chunk2[1] + offset2
                                            # print range1 , ' * ' , range2
                                            # print 1,range1
                                            # print 2,range2

                                            if range1[0] <= range2[1] and range2[0] <= range1[1]:
                                                # print " overlap!?"
                                                return True

                                        if chunk2[0] in [0, 2, 3, 7, 8]:  # CIGAR letters: MX=ND
                                            offset2 += chunk2[1]

                                # 01 - is silent mate (sa_tags[1]) overlappnig the current split read (read)
                                if sa_tags[1][4] == read_strand:
                                    offset2 = read.pos
                                    for chunk2 in read.cigar:
                                        if chunk2[0] == 0:
                                            range2 = offset2, chunk2[1] + offset2
                                            # print range1 , ' * ' , range2
                                            # print 1,range1
                                            # print 2,range2

                                            if range1[0] <= range2[1] and range2[0] <= range1[1]:
                                                # print " overlap!?"
                                                return True

                                        if chunk2[0] in [0, 2, 3, 7, 8]:  # CIGAR letters: MX=ND
                                            offset2 += chunk2[1]

                            if chunk1[0] in [0, 2, 3, 7, 8]:  # CIGAR letters: MX=ND
                                offset1 += chunk1[1]

                    # als HI:i:1 en HI:i:2 in tegengestelde strand zijn EN elkaar overlappen?
                    if read_strand != sa_tags[0][4]:
                        offset1 = sa_tags[0][1]

                        for chunk1 in cigar_to_cigartuple(sa_tags[0][2]):
                            if chunk1[0] == 0:
                                range1 = offset1, chunk1[1] + offset1

                                offset2 = read.pos
                                for chunk2 in read.cigar:
                                    if chunk2[0] == 0:
                                        range2 = offset2, chunk2[1] + offset2
                                        # print range1 , ' * ' , range2

                                        if range1[0] <= range2[1] and range2[0] <= range1[1]:
                                            # print " overlap!?"
                                            return True

                                    if chunk2[0] in [0, 2, 3, 7, 8]:  # CIGAR letters: MX=ND
                                        offset2 += chunk2[1]

                            if chunk1[0] in [0, 2, 3, 7, 8]:  # CIGAR letters: MX=ND
                                offset1 += chunk1[1]
                        # print "\n ---"

                return False

            # @todo Shouldn't this return either and Edge or None
            # @todo only skip if type == DiscordantMates/Reads, if Spanning then always allow (example 072?)
            if abs(pos1.get_dist(pos2, False)) < MAX_ACCEPTABLE_INSERT_SIZE:
                self.overlapping_reads += 1
                # print pos1, pos2, pos1.get_dist(pos2, False) , '<'  , MAX_ACCEPTABLE_INSERT_SIZE
            elif is_hybridization_artifact(read, pos1, pos2, parsed_SA_tags):
                self.crosshybridization_artifacts += 1
            else:
                n_match_p1, n_match_p2 = sum([x[1] for x in read.cigar if x[0] == 0]), sum([x[1] for x in parsed_cigar_tuple if x[0] == 0])

                return (pos1, pos2,
                        (read.cigarstring, parsed_SA_tag[2]),
                        (
                            bam_parse_alignment_offset(read.cigar, True),
                            bam_parse_alignment_offset(parsed_cigar_tuple, True)
                        ),
                        acceptor, donor,
                        n_match_p1, n_match_p2
                        )

            return (None, None, None, None, None, None, None, None)

        log.debug("Parsing reads to obtain fusion gene and splice junctions")

        with tqdm(total=self.pysam_fh.mapped) as pbar:
            for read in self.pysam_fh.fetch():
                _chr = self.pysam_fh.get_reference_name(read.reference_id)
                try:
                    sa = self.parse_SA(read.get_tag('SA'), _chr)
                except Exception:
                    raise ValueError("Problems parsing SA tag for:\n\t%s", str(read))
                rg = JunctionTypeUtils.enum(read.get_tag('RG'))

                pos1, pos2 = None, None

                if JunctionTypeUtils.is_fusion_junction(rg):
                    pos1, pos2, junction, ctps, acceptor, donor, n_match_p1, n_match_p2 = read_to_junction(read, rg, sa)
                    if pos1 is not None:
                        alignment_score = read.get_tag('AS')
                        # It happens that STAR assigns incorrenct alignment scores like -2147483577
                        # If this happens we set it to 12 as some kind of wild guess
                        if alignment_score < 0:
                            alignment_score = 12

                        p1, p2 = pos1._chr, pos2._chr
                        if p1 > p2:
                            p2, p1 = p1, p2

                        if p1 not in fusion_junctions:
                            fusion_junctions[p1] = {}
                        if p2 not in fusion_junctions[p1]:
                            fusion_junctions[p1][p2] = Graph()

                        fusion_junctions[p1][p2].insert_edge(pos1, pos2, rg, junction, ctps, acceptor, donor, read.get_tag('nM'), alignment_score, n_match_p1, n_match_p2)

                elif rg == JunctionTypes.silent_mate:  # Not yet implemented, may be useful for determining type of junction (exonic / intronic)
                    pass

                else:  # pragma: no cover
                    raise Exception("Unknown type read: %s. Was the alignment fixed with a more up to date version of Dr.Disco?" % JunctionTypeUtils.str(rg))

                # Detect introns, clipping and insertions /deletions by SAM flags
                for internal_edge in self.find_cigar_edges(read):
                    i_pos1, i_pos2 = None, None

                    if internal_edge[2] in [JunctionTypes.cigar_splice_junction, JunctionTypes.cigar_deletion]:  # , JunctionType.cigar_deletion
                        i_pos1 = BreakPosition(_chr, internal_edge[0], STRAND_FORWARD)
                        i_pos2 = BreakPosition(_chr, internal_edge[1], STRAND_REVERSE)

                        if internal_edge[2] in [JunctionTypes.cigar_deletion] and i_pos1.get_dist(i_pos2, False) < MAX_ACCEPTABLE_INSERT_SIZE:
                            i_pos1, i_pos2 = None, None

                    elif internal_edge[2] in [JunctionTypes.cigar_soft_clip]:
                        if pos1 is None or rg in [JunctionTypes.spanning_paired_1_s, JunctionTypes.spanning_paired_2_s]:
                            pass
                        elif JunctionTypeUtils.is_fusion_junction(rg):
                            i_pos1 = BreakPosition(_chr, internal_edge[0], pos2.strand)
                            i_pos2 = BreakPosition(_chr, internal_edge[1], pos1.strand)
                        else:  # pragma: no cover
                            raise Exception("what todo here %s " + JunctionTypeUtils.str(rg))
                    else:  # pragma: no cover
                        raise Exception("Edge type not implemented: %s\n\n%s", internal_edge, str(read))

                    if internal_edge[2] in [JunctionTypes.cigar_soft_clip, JunctionTypes.cigar_hard_clip]:
                        if i_pos1 is not None:
                            node = fusion_junctions[p1][p2].get_node_reference(i_pos2)
                            if node is not None:
                                node.add_clip()
                            # else condition happens when clips are found on positions that are not junctions
                    else:
                        if i_pos1 is not None:
                            if internal_edge[2] == JunctionTypes.cigar_splice_junction:
                                # splice_junctions.insert_splice_edge(i_pos1, i_pos2, internal_edge[2], None)
                                if _chr not in splice_junctions:
                                    splice_junctions[_chr] = Graph()
                                splice_junctions[_chr].insert_edge(i_pos1, i_pos2, internal_edge[2], None, None, None, None, 0, 0, -1, -1)
                            else:
                                if _chr not in fusion_junctions:
                                    fusion_junctions[_chr] = {}
                                if _chr not in fusion_junctions[_chr]:
                                    fusion_junctions[_chr][_chr] = Graph()
                                fusion_junctions[_chr][_chr].insert_edge(i_pos1, i_pos2, internal_edge[2], None, None, None, None, 0, 0, -1, -1)  # By insertion tags in cigar string

                pbar.update(1)

        log.debug("Alignment data loaded - overlapping reads excluded: " + str(self.overlapping_reads) + " - Putative crosshybridization artifact reads excluded: " + str(self.crosshybridization_artifacts))

    def parse_pos(self, str_pos):
        try:
            _chr, _poss = str_pos.split(":", 2)
            _poss = _poss.replace(",", "").split("-", 2)
        except Exception:
            raise ValueError("No understandable region format (chr:123-345): %s", str_pos)

        return str(_chr), int(_poss[0]), int(_poss[1])

    def extract(self, pos1, pos2, bamfile_out, restrict_to_targeted_chromosomes=False):
        """
        The bam-extract function
        """
        pos1 = self.parse_pos(pos1)
        pos2 = self.parse_pos(pos2)
        target_chrs = set([pos1[0],pos2[0]])

        ids = set()

        fh = pysam.AlignmentFile(bamfile_out, "wb", header=self.pysam_fh.header)

        log.info('Scanning read names in pos1')
        for r in tqdm(self.pysam_fh.fetch(pos1[0], pos1[1], pos1[2])):  # pysam.view(b,"chr21:1-500")
            ids.add(r.query_name)
        log.info('  Total unique reads: ' + str(len(ids)))

        log.info('Scanning read names in pos2')
        for r in tqdm(self.pysam_fh.fetch(pos2[0], pos2[1], pos2[2])):  # pysam.view(b,"chr21:1-500")
            ids.add(r.query_name)
        log.info('  Total unique reads: ' + str(len(ids)))


        if restrict_to_targeted_chromosomes:
            ids_other_chromosomes = set() # read id's for reads that also map to other chromosomes, to be exclused afterwards
            log.info('  Scanning for reads that also map to non targeted chromosomes')

            with tqdm(total=self.pysam_fh.mapped) as pbar:
                for r in tqdm(self.pysam_fh.fetch()):
                    if r.query_name in ids and r.reference_name not in target_chrs:
                        ids_other_chromosomes.add(r.query_name)
                    pbar.update(1)
                ids = ids.difference(ids_other_chromosomes)
            log.info('  Total reads to be excluded because of also mapping to other chromosomes:' + str(len(ids_other_chromosomes)))
            log.info('  Total unique reads: ' + str(len(ids)))


        log.info('Scanning through entire bam file looking for those reads')
        with tqdm(total=self.pysam_fh.mapped) as pbar:
            for r in tqdm(self.pysam_fh.fetch()):
                if r.query_name in ids:
                    fh.write(r)
                pbar.update(1)

        fh.close()
        pysam.index(str(bamfile_out))

    # static is elegant for functions that do not use class properties
    # http://stackoverflow.com /questions /18679803 /python-calling-method-without-self
    @staticmethod
    def parse_SA(SA_tag, _chr):
        """
        old data struct: [chr, pos, cigar, mapq, strand, Nm]
        new data struct: [chr, pos, strand, cigar, mapQ, Nm]
        convert [n->o]:  [0,   1,   3,      4,     2,    5 ]

        for backwards compatibility, keep old data structure internally in analysis code

        hence, if sa_tag[2] == '+' or sa_tag[2] == '-', we have the new type and need to reconvert it
        """
        sa_tags = SA_tag.split(";")

        for i in range(len(sa_tags)):
            sa_tags[i] = sa_tags[i].split(",")

            if sa_tags[i][2] == '+' or sa_tags[i][2] == '-':
                sa_tags[i] = [sa_tags[i][0], sa_tags[i][1], sa_tags[i][3], sa_tags[i][4], sa_tags[i][2], sa_tags[i][5]]
            sa_tags[i][1] = int(sa_tags[i][1])

            if sa_tags[i][0] == '=':
                sa_tags[i][0] = _chr

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
            the end. If it is before the first M /= flag, it's direction should
            be '-', otherwise '+' as it clips from the end. In principle the
            soft /hard clips are not edges but properties of nodes. One particular
            base may have several soft /hard clips (from different locations but
            adding weigt to the same node).

            Splice juncs are bi-directional, and real edges.

splice-junc:                           <=============>
            """

        tt = {
            2: JunctionTypes.cigar_deletion,
            3: JunctionTypes.cigar_splice_junction,
            4: JunctionTypes.cigar_soft_clip,
            5: JunctionTypes.cigar_hard_clip
        }

        offset = read.reference_start
        solid = False

        for chunk in read.cigar:
            if chunk[0] in tt:  # D N S H
                """Small softclips occur all the time - introns and
                deletions of that size shouldn't add weight anyway

                1S 5M means:
                startpoint-1, startpoint => Softclip

                5M 2S means:
                startpoint +5, startpoint +5+2 => softclip

                In the first case, which I call left_clipping, the junction
                occurs before the start position and it needs to be corrected
                for

                @todo it's still a buggy implementation because the following
                is theoretically possible too:

                5H 10S 5M

                in that case, the first edge should be -15,-10 and the
                second -10, 0. Maybe soft- and hard clipping should
                be merged together?
                """

                if chunk[0] in [4, 5] and not solid:
                    offset -= chunk[1]

                if chunk[1] > MAX_ACCEPTABLE_ALIGNMENT_ERROR:
                    if solid:
                        """Clips to the first node:
                        M M M M M M M M M M S S S S S
                                           <========]
                        """
                        # yield (offset, (offset + chunk[1]), tt[chunk[0]])
                        yield ((offset + chunk[1]), offset, tt[chunk[0]])
                    else:
                        """Clips to the second node:
                        S S S S S M M M M M M M M M M
                        [========>
                        """
                        yield (offset, (offset + chunk[1]), tt[chunk[0]])

            if chunk[0] not in [4, 5]:
                solid = True

            offset += chunk[1]


class IntronDecomposition:
    def __init__(self, alignment_file):
        self.alignment_file = alignment_file

    def decompose(self, MIN_SCORE_FOR_EXTRACTING_SUBGRAPHS):
        alignment = BAMExtract(self.alignment_file, True)

        # initally worked, but slow because of same insane nodes
        fusion_junctions = {}  # fusion_junctions['chr1']['chr3'] = Graph() - where first chr < second chr
        splice_junctions = {}  # [chr1] = Graph, etc.

        alignment.extract_junctions(fusion_junctions, splice_junctions)

        # Calc progress bar
        iterations = 0
        for chr1 in fusion_junctions:
            for chr2 in fusion_junctions[chr1]:
                fusion_junctions[chr1][chr2].generate_edge_idx()
                fusion_junctions[chr1][chr2].check_symmetry()
                iterations += len(fusion_junctions[chr1][chr2].edge_idx)

        # prune
        k = 0
        log.info("Merging edges by splice junctions")
        with tqdm(total=iterations) as pbar:
            thicker_edges = {}
            for chr1 in fusion_junctions:
                thicker_edges[chr1] = {}
                for chr2 in fusion_junctions[chr1]:
                    thicker_edges[chr1][chr2] = fusion_junctions[chr1][chr2].prune(pbar)  # Makes edge thicker by lookin in the ins. size - make a sorted data structure for quicker access - i.e. sorted list
                    k += len(thicker_edges[chr1][chr2])
        log.info("Pruned into " + str(k) + " candidate edge(s)")

        iterations = 0
        log.info("Rejoining splice junctions & re-insterting edges in order to extract subgraphs [MIN_SCORE_FOR_EXTRACTING_SUBGRAPHS=%i]" % MIN_SCORE_FOR_EXTRACTING_SUBGRAPHS)
        subnets = {}
        k = 0
        for chr1 in tqdm(fusion_junctions):
            subnets[chr1] = {}
            for chr2 in tqdm(fusion_junctions[chr1], position=0):  # if position is not set to 0, it indents and overlaps with the log.info messages
                unique_chrs = list(set([chr1, chr2]))
                k += fusion_junctions[chr1][chr2].rejoin_splice_juncs(splice_junctions, unique_chrs)  # Merges edges by splice junctions and other junctions
                fusion_junctions[chr1][chr2].reinsert_edges(thicker_edges[chr1][chr2])  # @todo move function into statuc function in this class

                # fusion_junctions.print_chain()

                fusion_junctions[chr1][chr2].reindex_sj()
                subnets[chr1][chr2] = fusion_junctions[chr1][chr2].extract_subnetworks_by_splice_junctions(thicker_edges[chr1][chr2], MIN_SCORE_FOR_EXTRACTING_SUBGRAPHS)
                iterations += len(subnets[chr1][chr2])
        log.info("Rejoined " + str(k) + " splice junction(s)")
        log.info("Extracted %i subnetwork(s)" % iterations)

        k, m = 0, 0
        log.info("Merging overlapping subgraphs based on genomic distance")
        with tqdm(total=iterations) as pbar:
            self.results = []
            for chr1 in fusion_junctions:
                for chr2 in fusion_junctions[chr1]:
                    #  @todo unclear why graph is combined with thicker edges - makes no sense as this may duplice things at first glance?
                    results = self.merge_overlapping_subnets(subnets[chr1][chr2], pbar)
                    k += results[0]
                    results = results[1]
                    m += len(results)
                    for _ in results:
                        _.classify_intronic_exonic(splice_junctions, unique_chrs)
                    self.results += results
        log.info("Merged " + str(k) + " of the " + str(iterations) + " into " + str(m) + " merged subnetwork(s)")

        return len(self.results)

    def export(self, fh):
        log.info("Exporting results")
        fh.write(str(self))

    def __str__(self):
        ordered = []
        for subnet in self.results:
            ordered.append((subnet, subnet.total_score, subnet.get_overall_entropy()))
        
        # keys = total_score, overall_entropy
        ordered = [subnet[0] for subnet in sorted(ordered, key=operator.itemgetter(1, 2), reverse=True)]

        return ("chr-A"            "\t" "pos-A"             "\t" "direction-A"   "\t"
                "pos-A-acceptor"   "\t" "pos-A-donor"       "\t"
                "chr-B"            "\t" "pos-B"             "\t" "direction-B"   "\t"
                "pos-B-acceptor"   "\t" "pos-B-donor"       "\t"
                "genomic-distance" "\t" "filter-status"     "\t" "circRNA"       "\t" "intronic/exonic"    "\t"
                "score"            "\t" "soft+hardclips"    "\t" "n-split-reads" "\t" "n-discordant-reads" "\t"
                "alignment-score"  "\t" "mismatches"        "\t"
                "n-edges"          "\t" "n-nodes-A"         "\t" "n-nodes-B"     "\t"
                "n-splice-junc-A"  "\t" "n-splice-junc-B"   "\t"
                "entropy-bp-edge"  "\t" "entropy-all-edges" "\t"
                "bp-pos-stddev"    "\t" "entropy-disco-bps" "\t"
                "unique-bp-stddev-A" "\t" "entropy-disco-bps-A" "\t" "unique-bp-stddev-B" "\t" "entropy-disco-bps-B" "\t"
                "lr-A-slope"       "\t" "lr-A-intercept"    "\t" "lr-A-rvalue"   "\t" "lr-A-pvalue"        "\t" "lr-A-stderr" "\t"
                "lr-B-slope"       "\t" "lr-B-intercept"    "\t" "lr-B-rvalue"   "\t" "lr-B-pvalue"        "\t" "lr-B-stderr" "\t"
                "disco/split"      "\t" "clips/score"       "\t" "nodes/edge"    "\t"
                "median-AS-A"      "\t" "median-AS-B"       "\t" "max-AS-A"      "\t" "max-AS-B"           "\t"
                "data-structure"   "\n"
                "%s" % (''.join([str(subnet) for subnet in ordered])))

    def merge_overlapping_subnets(self, subnets, pbar):
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

        def sq_dist(vec):
            sum_of_squares = sum(pow(x, 2) for x in vec)

            if sum_of_squares > 0:
                avg_sq_d = float(sum_of_squares) / len(vec)
            else:
                avg_sq_d = 0.0

            return math.sqrt(avg_sq_d)

        def tree_insert(genometree, pos, entry):
            """inserts entry into a set at pos in genometree"""
            position_accession = HTSeq.GenomicPosition(pos._chr, pos.pos)
            position_tree = genometree[position_accession]
            if position_tree:  # position was not found in tree, insert completely
                if pos.strand in position_tree:
                    position_tree[pos.strand].add(entry)
                else:
                    position_tree[pos.strand] = set([entry])
            else:
                genometree[position_accession] = {pos.strand: set([entry])}

        def tree_remove(genometree, subnet):
            positions = set()
            for edge in subnet.edges:
                positions.add(edge.origin.position)
                positions.add(edge.target.position)

            for pos in positions:
                position_accession = HTSeq.GenomicPosition(pos._chr, pos.pos)
                position_tree = genometree[position_accession]
                if position_tree and pos.strand in position_tree:
                    sset = position_tree[pos.strand]
                    if subnet in sset:
                        sset.remove(subnet)

        idx_l = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        idx_r = HTSeq.GenomicArrayOfSets("auto", stranded=False)

        for subnet in subnets:
            for edge in subnet.edges:
                tree_insert(idx_l, edge.origin.position, subnet)
                tree_insert(idx_r, edge.target.position, subnet)

        k = 0
        new_subnets = []
        subnets.reverse()

        previous = len(subnets)
        while subnets:
            subnet = subnets.pop()
            tree_remove(idx_l, subnet)
            tree_remove(idx_r, subnet)

            to_be_merged_with = set()
            candidates = set()

            for edge in subnet.edges:
                # based on l-node to be close
                pos = edge.origin.position
                for step in idx_l[HTSeq.GenomicInterval(pos._chr, max(0, pos.pos - MAX_ACCEPTABLE_INSERT_SIZE), pos.pos + MAX_ACCEPTABLE_INSERT_SIZE)].steps():
                    if step[1] and pos.strand in step[1]:
                        for candidate_subnet in step[1][pos.strand]:
                            if subnet != candidate_subnet:
                                candidates.add(candidate_subnet)

                pos = edge.target.position
                for step in idx_r[HTSeq.GenomicInterval(pos._chr, max(0, pos.pos - MAX_ACCEPTABLE_INSERT_SIZE), pos.pos + MAX_ACCEPTABLE_INSERT_SIZE)].steps():
                    if step[1] and pos.strand in step[1]:
                        for candidate_subnet in step[1][pos.strand]:
                            if subnet != candidate_subnet:
                                candidates.add(candidate_subnet)

            for candidate_subnet in candidates:
                new_merged = False
                l_dists, r_dists = subnet.find_distances(candidate_subnet)

                if l_dists:
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
                    else:  # @todo revise if / else / elif logic, and use linear regression or sth like that
                        if n_l_dist > 0:
                            l_dist_ins_ratio = 1.0 * n_l_dist / len(l_dists)
                            if l_dist_ins_ratio == 1.0 and rmsq_r_dist < 15000:
                                new_merged = True
                            elif l_dist_ins_ratio > 0.7 and rmsq_r_dist < 10000:
                                new_merged = True
                            elif l_dist_ins_ratio > 0.3 and rmsq_r_dist < 5000:
                                new_merged = True

                        if n_r_dist > 0:
                            r_dist_ins_ratio = 1.0 * n_r_dist / len(r_dists)
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
                subnets.remove(candidate_subnet)  # [subnets.index(candidate_subnet)] = None# @todo inverse subnets and use .pop() and .remove()
                del(candidate_subnet)
                k += 1

            new_subnets.append(subnet)

            pbar.update(previous - len(subnets))
            previous = len(subnets)

        return (k, new_subnets)
