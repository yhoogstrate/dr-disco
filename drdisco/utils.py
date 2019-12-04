#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import gzip
import re


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


