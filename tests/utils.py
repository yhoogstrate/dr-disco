#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Dr. Disco

[License: GNU General Public License v3 (GPLv3)]

 This file is part of Dr. Disco.

 FuMa is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Dr. Disco is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
"""


import unittest
import filecmp
import os
import pysam
import subprocess
from drdisco.ChimericAlignment import ChimericAlignment


def sam_to_fixed_bam(sam, fixed_bam, T_TEST_DIR):
    basename, ext = os.path.splitext(os.path.basename(sam))
    bam_file = T_TEST_DIR + '/' + basename + '.bam'

    fhq = open(bam_file, "wb")
    fhq.write(pysam.view('-b', sam))
    fhq.close()

    alignment_handle = ChimericAlignment(bam_file)
    return alignment_handle.convert(fixed_bam, T_TEST_DIR)


def bam_diff(f1, f2, T_TEST_DIR):
    basename, ext = os.path.splitext(os.path.basename(f1))

    f1sorted = T_TEST_DIR + basename + '.f1.sorted.bam'
    f2sorted = T_TEST_DIR + basename + '.f2.sorted.bam'

    pysam.sort(f1, '-n', '-o', f1sorted)
    pysam.sort(f2, '-n', '-o', f2sorted)

    f1sam = T_TEST_DIR + basename + '.f1.sam'
    f2sam = T_TEST_DIR + basename + '.f2.sam'

    fhq = open(f1sam, "w")
    fhq.write(pysam.view('-h', f1sorted))
    fhq.close()

    fhq = open(f2sam, "w")
    fhq.write(pysam.view('-h', f2sorted))
    fhq.close()

    subprocess.Popen(['sed', '-i', '-r', 's@(SA:[^\\t]+)\\t(LB:[^\\t]+)\t(RG:[^\\t]+)@\\3\\t\\1\\t\\2@', f2sam], stdout=subprocess.PIPE).stdout.read()

    subprocess.Popen(['sed', '-i', '-r', 's@\\tFI:i:[0-9]+@@', f1sam], stdout=subprocess.PIPE).stdout.read()
    subprocess.Popen(['sed', '-i', '-r', 's@\\tFI:i:[0-9]+@@', f2sam], stdout=subprocess.PIPE).stdout.read()

    # one time only
    # subprocess.Popen(['sed', '-i' , '-r', 's@\\tSA:Z:[^\\t]+@@', f1sam], stdout=subprocess.PIPE).stdout.read()
    # subprocess.Popen(['sed', '-i' , '-r', 's@\\tSA:Z:[^\\t]+@@', f2sam], stdout=subprocess.PIPE).stdout.read()

    return filecmp.cmp(f1sam, f2sam), f1sam, f2sam


def get_diff(test_file, output_file):
	out = "diff '" + test_file + "' '" + output_file + "':\n"
	p = subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE)
	s = str(p.communicate()[0])
	return out + s

def main():
    unittest.main()
