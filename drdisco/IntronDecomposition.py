#!/usr/bin/env python3
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

#http://www.samformat.info/sam-format-flag

import logging
import pysam

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
    
    def decompose(self,alignment_file):
        pysam_fh = self.test_disco_alignment(alignment_file)
        #self.break_point
