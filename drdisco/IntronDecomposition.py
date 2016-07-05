#!/usr/bin/env python3
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

#http://www.samformat.info/sam-format-flag

import logging
import pysam


class Element:
    def __init__(self,_chr,_pos):
        self._chr = _chr
        self._pos = _pos
        arcs_fwd = []
        arcs_rev = []


class Chain:
    def __init__(self):
        idx = {}
    
    def insert(self,pos1,pos2,_type):
        # Type:
        # 1. 'N' <- N alignment flag in SAM, meaning INTRON
        # 2. 'discordant_read'
        # 3. 'split_read'
        # 4. 'silent_mate'
        #
        # ... 'S' and 'H' for soft and hard clipping?
        pass
    
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
            #sa_tags[i][2] = int(sa_tags[i][2])
        
        return sa_tags
    
    def decompose(self,alignment_file):
        pysam_fh = self.test_disco_alignment(alignment_file)
        
        lpos = self.break_point.get_left_position(True)
        rpos = self.break_point.get_right_position(True)
        
        chain_left = []
        chain_right = []
        
        c = Chain()
        
        for r in pysam_fh.fetch(lpos[0],lpos[1]):
            sa = self.parse_SA(r.get_tag('SA'))
            _chr = pysam_fh.get_reference_name(r.reference_id)
            insert_size = self.get_insert_size([_chr,r.reference_start],[sa[0][0],sa[0][1]])
            
            if abs(insert_size) >= 1.5*126:
                print r.query_name
                print " -" , r.reference_start, "-" if r.is_reverse else "+"
                print " -" , r.cigarstring
                print " -" , r.get_tag('RG')
                print " -" , sa
                
                if r.get_tag('RG') == 'discordant_mates':
                    print ' *', self.get_insert_size([_chr,r.reference_start],[sa[0][0],sa[0][1]])
                #print r.get_aligned_pairs()
            else:
                # Figure out whether there was hard or soft clipping on one of the mates
                pass
