#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import logging
import pysam

# load cfg
from . import *

class BAMExtract:
    def __init__(self, bam_file):
        self.samfile = self.test_idx(bam_file)
    
    def test_idx(self, bam_file):
        pysam_fh = pysam.AlignmentFile(bam_file, "rb")
        
        try:
            pysam_fh.fetch()
        except:
            logging.info('Indexing BAM file with pysam: '+pysam_fh.filename)
            pysam.index(pysam_fh.filename)
            pysam_fh = pysam.AlignmentFile(pysam_fh.filename)
        
        try:
            pysam_fh.fetch()
        except:
            raise Exception('Could not indexing BAM file: '+pysam_fh.filename)
        
        return pysam_fh
    
    def parse_pos(self, str_pos):
        _chr,_poss = str_pos.split(":",2)
        _poss = _poss.replace(",","").split("-",2)
        
        return str(_chr), int(_poss[0]), int(_poss[1])

    def extract(self, pos1, pos2, bamfile_out):
        pos1 = self.parse_pos(pos1)
        pos2 = self.parse_pos(pos2)
        
        ids = set()
        
        fh = pysam.AlignmentFile(bamfile_out, "wb", header=self.samfile.header)
        
        for r in self.samfile.fetch(pos1[0],pos1[1],pos1[2]):#pysam.view(b,"chr21:1-500")
            ids.add(r.query_name)
        
        for r in self.samfile.fetch(pos2[0],pos2[1],pos2[2]):#pysam.view(b,"chr21:1-500")
            ids.add(r.query_name)
        
        for r in self.samfile.fetch():
            if r.query_name in ids:
                fh.write(r)
        
        
        fh.close()
        
        pysam.index(str(bamfile_out))
    
    # static is elegant for functions that do not use class properties
    # http://stackoverflow.com/questions/18679803/python-calling-method-without-self
    @staticmethod
    def parse_SA(SA_tag):
        sa_tags = SA_tag.split(";")
        for i in range(len(sa_tags)):
            sa_tags[i] = sa_tags[i].split(",")
            sa_tags[i][1] = int(sa_tags[i][1])
        
        return sa_tags
    
    @staticmethod
    def find_cigar_arcs(read):
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
                
                in that case, the first arc should be -15,-10 and the
                second -10,0. Maybe soft- and hard clipping should
                be merged together?
                """
                
                if chunk[0] in [4,5] and not solid:
                    offset -= chunk[1]

                if chunk[1] > SPLICE_JUNC_ACC_ERR:
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
