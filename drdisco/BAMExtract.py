#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import logging
import pysam
from .BAMExtract import *


class BAMExtract:
    def __init__(self, bam_file):
        self.samfile = pysam.AlignmentFile(bam_file, "rb")
    
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

