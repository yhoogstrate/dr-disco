#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
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

import drdisco
from drdisco.IntronDecomposition import BAMExtract

import unittest
import logging
import sys
import filecmp
import pysam
import os

logging.basicConfig(level=logging.DEBUG, format=drdisco.__log_format__, stream=sys.stdout)

TEST_DIR = "tests/bam-extract/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01a(self):
        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"
        test_file = TEST_DIR + "test_terg_02.filtered.sam"

        c = BAMExtract(input_file)
        c.extract("chr21:39000000-40000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))

    def test_01b(self):
        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"
        test_file = TEST_DIR + "test_terg_02.filtered.sam"

        c = BAMExtract(input_file)
        c.extract("chr5:1-2", "chr21:39000000-40000000", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))

    def test_02a(self):
        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"
        test_file = TEST_DIR + "test_terg_02.filtered.sam"

        c = BAMExtract(input_file)
        c.extract("chr7:151000000-153000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))

    def test_02b(self):
        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"
        test_file = TEST_DIR + "test_terg_02.filtered.sam"

        c = BAMExtract(input_file)
        c.extract("chr5:1-2", "chr7:151000000-153000000", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))

    def test_03(self):
        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"

        c = BAMExtract(input_file)
        c.extract("chr12:151000000-153000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        with open(output_file_s, "r") as fh:
            self.assertEqual(fh.read(), "")  # empty file check


def main():
    unittest.main()

if __name__ == '__main__':
    main()
