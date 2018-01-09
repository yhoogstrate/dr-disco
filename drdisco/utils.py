#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import gzip


alt_map = {'ins': '0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def reverse_complement(seq):
    seq = seq.upper()
    for k, v in alt_map.iteritems():
        seq = seq.replace(k, v)
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    for k, v in alt_map.iteritems():
        bases = bases.replace(v, k)
    return bases


def is_gzip(filename):
    try:
        f = gzip.GzipFile(filename, 'rb')
        f.read()
        return True
    except Exception:
        return False
