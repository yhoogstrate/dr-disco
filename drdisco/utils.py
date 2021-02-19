#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import gzip
import re
import pysam
from drdisco import log


alt_map = {'ins': '0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def reverse_complement(seq):
    seq = seq.upper()
    for k, v in alt_map.items():
        seq = seq.replace(k, v)
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    for k, v in alt_map.items():
        bases = bases.replace(v, k)

    return bases


def is_gzip(filename):
    try:
        f = gzip.GzipFile(filename, 'rb')
        f.read()
        return True
    except Exception:
        return False


prog1 = re.compile('^([0-9]+)$')
prog2 = re.compile('^([0-9]+)\(([+\-\.\?])\)$')


def parse_pos(strpos):
    """
    chr1:12345 -> ['chr1', 12345, '.']
    chr1:12345(-) -> ['chr1', 12345, '-']
    """
    val = [None, None, '.']
    strpos = strpos.replace(',', '')

    params = strpos.split(':', 2)
    if len(params) != 2:
        raise ValueError("Invalid pos: '" + strpos + "' needs to be formatted as 'chr1:12,345' or 'chr1:12345(+)'")

    m1 = prog1.match(params[1])
    m2 = prog2.match(params[1])
    if m1:
        val[0] = str(params[0])
        val[1] = int(m1.group(1))
    elif m2:
        val[0] = str(params[0])
        val[1] = int(m2.group(1))
        val[2] = str(m2.group(2).replace('?', '.'))
    else:
        raise ValueError("Invalid pos: '" + strpos + "' needs to be formatted as 'chr1:12,345' or 'chr1:12345(+)'")

    return val


def str_to_bytearray(s):
    b = bytearray()
    b.extend(map(ord, s))
    return b


def get_drdisco_version(alignment_file):
    """
    dr-disco fix <= v0.18.1: 0-based SA tags
    dr-disco fix >= v0.18.2: 1-based SA tags
    """

    drdisco_version = None # the version of dr-disco used to perform `dr-disco fix`

    with pysam.AlignmentFile(alignment_file, "rb") as bam_fh:
        if 'PG' in bam_fh.header:
            for pg in bam_fh.header['PG']:
                if 'VN' in pg and 'ID' in pg and pg['ID'] == 'drdisco_fix_chimeric':
                    drdisco_version = pg['VN'].split(".")
                    
                    for i in range(len(drdisco_version)):
                        try:
                            drdisco_version[i] = int(drdisco_version[i])
                        except:
                            raise ValueError("Inconsistent Dr. Disco version: " + str(pg) )
    
    if drdisco_version == None:
        log.warning("No Dr. Disco version in fixed bam file detected? Using v0.0.0 as fallback.")
    else:
        log.info("Dr. Disco version of fixed bam file: v" + ".".join([str(_) for _ in drdisco_version]))
    
    return drdisco_version


def drdisco_version_leq(actual_version, requirement_version):
    """
    returns true if the actual verions is larger or equal to the requirement version
    [0,18,2] >= [0,0,0] = T
    [1,18,2] >= [0,0,0] = T
    [0,18,2] >= [1,0,0] = F
    [0, 0,0] >= [0,0,0] = T
    """

    if len(actual_version) != len(requirement_version):
        raise ValueError("Inconsistent Dr. Disco version: " + str(pg) )

    for i in range(len(actual_version)):
        if actual_version[i] > requirement_version[i]:
            return True
        elif actual_version[i] < requirement_version[i]:
            return False
        # else : equal, skip to next
    
    return True
