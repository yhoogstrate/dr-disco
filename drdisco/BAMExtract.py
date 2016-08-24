#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import logging
import pysam

def parse_pos(str_pos):
    _chr,_poss = str_pos.plit(":",2)
    _poss = _poss.replace(",","").split("-",2)
    
    return _chr, int(_poss[0]), int(_poss[1])

def get_ids(bam,pos1,pos2):
    pos1 = parse_pos(pos1)
    pos2 = parse_pos(pos2)
    
    ids = set()
    
    samfile = pysam.AlignmentFile(bam, "rb")
    for r in samfile.fetch(pos1[0],pos1[1],pos1[2]):#pysam.view(b,"chr21:1-500")
        ids.add(r.query_name)
    
    for r in samfile.fetch(pos2[0],pos2[1],pos2[2]):#pysam.view(b,"chr21:1-500")
        ids.add(r.query_name)
    
    for r in samfile.fetch():
        if r.query_name in ids:
            print q

