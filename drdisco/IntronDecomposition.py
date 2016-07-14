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

pat_bam_parse_alignment_offset_using_cigar = re.compile("([0-9]+)([MIDNSHPX=])")
def cigar_to_cigartuple(cigar_str):
    tt = {'M':0,#	BAM_CMATCH	0
          'I':1,#	BAM_CINS	1
          'D':2,#	BAM_CDEL	2
          'N':3,#	BAM_CREF_SKIP	3
          'S':4,#	BAM_CSOFT_CLIP	4
          'H':5,#	BAM_CHARD_CLIP	5
          'P':6,#	BAM_CPAD	6
          '=':7,#	BAM_CEQUAL	7
          'X':8#	BAM_CDIFF	8
          }
    
    cigartup = []
    
    for chunk in pat_bam_parse_alignment_offset_using_cigar.finditer(cigar_str):
        flag = chunk.group(2)
        length = chunk.group(1)
        cigartup.append((tt[flag],int(length)))
    
    return cigartup



class CigarAlignment:
    """
    STAR Fusion is not clear about which chunk is the one that started
    the read. We can calculate this back based on the CIGAR string.
    
    Usually you find a discordant pair as follows:
    
    part1: 46S80M
    part2: 80M46s
    
    Here the M and S are exectly complementary to each other. It does
    however occur that additional soft/hard clips are added. Hence we
    have to align them to each other and find the one with the least
    distance (in cases above r^2 = 0).
    
    We use matrices in a alignment context:
    
    tup1: 2S55M69S
    tup2: 57S69M
    
    Intiation:
    
                   2S     55M     69S
        ||======||=====||======||======||
        || 0    || 2^2 || 55^2 || 69^2 ||
        ||======||=====||======||======||
    57S || 57^2 ||      |       |       |
        ||======||------+-------+-------+
    69M || 69^2 ||      |       |       |
        ||======||------+-------+-------+
    
    fill function:
    i,j = min(
       (tup1[i] - tup2[j])   +   matrix[i-1,j-1]
       matrix[i-1]              , i >= 1
       matrix[j-1]              , j >= 1
    )
    
                   2S     55M     69S
        ||======||=====||======||======||
        || 0    || 2^2 || 55^2 || 69^2 ||
        ||======||=====||======||======||
    57S || 57^2 || i1   |  i2   |  i3   |
        ||======||------+-------+-------+
    69M || 69^2 || i2   |  i3   |  i4   |
        ||======||------+-------+-------+
        
    """
    def __init__(self,cigtup1,cigtup2):
        self.cigtup1 = cigtup1
        self.cigtup2 = cigtup2
        
        self.m = len(self.cigtup1)
        self.n = len(self.cigtup2)
        
        self.init_matrix()
    
    def init_matrix(self):
        #self.matrix = [[0] * (self.n+1) ] * (self.m + 1) << copies references of lists >,<
        self.matrix = [[0 for i in range(self.n+1)] for i in range(self.m+1)]
        
        for i in range(1,self.m+1):
            self.matrix[i][0] = pow(self.cigtup1[i-1][1],2)
        
        for j in range(1,self.n+1):
            self.matrix[0][j] = pow(self.cigtup2[j-1][1],2)
        
        
        
        self.tb_matrix = [["?" for i in range(self.n+1)] for i in range(self.m+1)]
        self.tb_matrix[0][0] = 't' # Terminate traceback; finished
        
        for i in range(1,self.m+1):
            self.tb_matrix[i][0] = "|"
        
        for j in range(1,self.n+1):
            self.tb_matrix[0][j] = "-"
    
    def print_tb_matrix(self):
        for j in range(self.n+1):
            for i in range(self.m+1):
                print self.tb_matrix[i][j],
            print
        print
    
    def print_matrix(self):
        for j in range(self.n+1):
            for i in range(self.m+1):
                print str(self.matrix[i][j])+"\t",
            print
        print


    def get_diagonal(self,diagonal):
            # 0,0
            # 1,0   0,1
            # 2,0   1,1   0,2

            i = diagonal
            j = 0

            for block in range(diagonal+1):
                if i < self.m and j < self.n:
                    yield (i,j)
                
                i -= 1
                j += 1
    
    def cigar_diff(self,cig_chunk1,cig_chunk2):
        """
        Maybe sort on cigar type, to reduce if statements
        """
        if (cig_chunk1[0] == 0 and cig_chunk2[0] == 0) or (cig_chunk1[0] == 4 and cig_chunk2[0] == 4):
            # M * M or S*S => sqaure 
            return pow(cig_chunk1[1] + cig_chunk2[1],2)
        
        elif (cig_chunk1[0] == 0 and cig_chunk2[0] == 4) or (cig_chunk1[0] == 4 and cig_chunk2[0] == 0):
            return pow((cig_chunk1[1] - cig_chunk2[1]),2)
        
        else:
            raise Exception("Not yet implemented")
        
    
    def calc_diff(self,i,j):
#                   2S     55M     69S
#        ||======||=====||======||======||
#        || 0    || 2^2 || 55^2 || 69^2 ||
#        ||======||=====||======||======||
#    57S || 57^2 ||      |       |       |
#        ||======||------+-------+-------+
#    69M || 69^2 ||      |       |       |
#        ||======||------+-------+-------+
        
        #c_ins_1 = pow(self.cigtup1[i][1],2)
        #c_ins_2 = pow(self.cigtup2[j][1],2)
        c_ins_1 = self.matrix[i-1][j]     + pow(self.cigtup1[i][1],2)# Insertion vertical, in tup1
        c_ins_2 = self.matrix[i][j-1]     + pow(self.cigtup2[j][1],2)# Insertion horizontal, in tup2
        c_diag = self.matrix[i][j] + self.cigar_diff(self.cigtup1[i],self.cigtup2[j])
        
        if c_ins_1 < c_diag:
            _type = "|"
            _min = c_ins_1
        else:
            _type = "m"
            _min = c_diag
        
        if c_ins_2 < _min:
            _type = "-"
            _min = c_ins_2
        
        return (_min, _type)
    
    def fill_matrix(self):
        n_diagonals = (self.n + self.m) - 1
        for diagonal in range(n_diagonals):
            for cell in self.get_diagonal(diagonal):
               
                diff, _type = self.calc_diff(cell[0],cell[1])
                self.tb_matrix[cell[0]+1][cell[1]+1] = _type
                self.matrix[cell[0]+1][cell[1]+1] = diff
                
        #self.print_tb_matrix()
        #self.print_matrix()
    
    def traceback_matrix(self):
        tup1_rev = []
        tup2_rev = []
        
        i = self.m + 1
        j = self.n + 1
        
        action = "s"
        while action != "t":
            if i < 0 and j < 0:
                raise Exception("Unpredicted out of bound")
            
            if action == "s":
                i -= 1 
                j -= 1
            
            elif action == "m":
                tup1_rev.append([action,self.cigtup1[i-1]])
                tup2_rev.append([action,self.cigtup2[j-1]])
                
                i -= 1 
                j -= 1
            
            elif action == "-":
                tup1_rev.append([action,(-1,0)])
                tup2_rev.append([action,self.cigtup2[j-1]])
                
                j -= 1
            
            elif action == "|":
                tup1_rev.append([action,self.cigtup1[i-1]])
                tup2_rev.append([action,(-1,0)])
                
                i -= 1
            
            action = self.tb_matrix[i][j]
        
        return tup1_rev[::-1], tup2_rev[::-1]

    def get_order(self):
        # 1. align
        #   a. fill
        self.fill_matrix()
        #   b. trace back
        self.traceback_matrix()
        
        # 2. calculate M / S only on those chunks that are aligned (M to S flags and vice versa)


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
    
    def insert_arc(self,arc,arc_type):
        skey = str(arc._target)
        
        if not self.arcs.has_key(skey):
            self.arcs[skey] = Arc(skey)
        
        self.arcs[skey].add_type(arc_type)
    
    def add_arc(self,node2,arc_type,do_vice_versa=True):
        if do_vice_versa:
            node2.add_arc(self,arc_type,False)
        
        arc = Arc(node2)
        self.insert_arc(arc,arc_type)
    
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
        self.idx = {}
        self.pysam_fh = pysam_fh
    
    def insert_entry(self,pos1,pos2,_type):
        """
         - Checks if Node exists at pos1, otherwise creates one
         - Checks if Node exists at pos2, otherwise creates one
         - Checks if Arc exists between them
        """
        
        key1 = str(pos1)
        if not self.idx.has_key(key1):
            self.idx[key1] = Node(pos1)
        
        key2 = str(pos2)
        if not self.idx.has_key(key2):
            self.idx[key2] = Node(pos2)
        
        self.idx[key1].add_arc(self.idx[key2],_type)
    
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
            
            #self.insert_entry(pos1,pos2,rg)
        elif rg in ["spanning_paired"]:
            #cig1 = cigar_to_cigartuple("2S55M69S")
            #cig2 = cigar_to_cigartuple("57S69M")
            
            #ca = CigarAlignment(cig1, cig2)
            #ca.get_order()
            
            
            #cig1 = cigar_to_cigartuple("57S69M")
            #cig2 = cigar_to_cigartuple("2S55M69S")
            
            #ca = CigarAlignment(cig1, cig2)
            #ca.get_order()
            
            cig1 = cigar_to_cigartuple("57S69M")
            cig2 = cigar_to_cigartuple("55M69S2M")
            
            #cig1 = cigar_to_cigartuple("1M57S69M")
            #cig2 = cigar_to_cigartuple("55M69S2M")
            
            
            ca = CigarAlignment(cig1, cig2)
            ca.get_order()
            
            import sys
            sys.exit()
        else:
            raise Exception("Fatal Error, RG: "+rg)
    
    def prune(self):
        for key in sorted(self.idx):
            print key, self.idx[key].str2()
            
            #print self.idx[key].arcs
            #print str(self.idx[key].arcs)
            #print [str(x) for x in self.idx[key].arcs.values()]
            
            #print "k:",key, " | ".join([str(arc) for arc in self.idx[key].arcs])
            #str(self.idx[key].arcs)
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
        c.prune()

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

