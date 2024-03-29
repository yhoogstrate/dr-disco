#!/usr/bin/env python
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

from drdisco.IntronDecomposition import BAMExtract

import unittest
import filecmp
import pysam
import os

from tests.utils import main, get_diff


TEST_DIR = "tests/bam-extract/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01_a(self):
        # Tests a file that has not (yet) been fixed with `dr-disco fix`

        input_file = TEST_DIR + "test_01_terg.bam"
        output_file = T_TEST_DIR + "test_01_terg.filtered.bam"
        output_file_s = T_TEST_DIR + "test_01_terg.filtered.sam"
        test_file = TEST_DIR + "test_02_terg.filtered.sam"

        c = BAMExtract(input_file, False)
        c.extract("chr21:39000000-40000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_02_a(self):
        input_file = TEST_DIR + "test_02_terg.bam"
        output_file = T_TEST_DIR + "test_02_terg.filtered.bam"
        output_file_s = T_TEST_DIR + "test_02_terg.filtered.sam"
        test_file = TEST_DIR + "test_02_terg.filtered.sam"

        c = BAMExtract(input_file, False)
        c.extract("chr21:39000000-40000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_02_b(self):
        input_file = TEST_DIR + "test_02_terg.bam"
        output_file = T_TEST_DIR + "test_02_terg.filtered.bam"
        output_file_s = T_TEST_DIR + "test_02_terg.filtered.sam"
        test_file = TEST_DIR + "test_02_terg.filtered.sam"

        c = BAMExtract(input_file, False)
        c.extract("chr5:1-2", "chr21:39000000-40000000", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_02_c(self):
        input_file = TEST_DIR + "test_02_terg.bam"
        output_file = T_TEST_DIR + "test_02_terg.filtered.bam"
        output_file_s = T_TEST_DIR + "test_02_terg.filtered.sam"
        test_file = TEST_DIR + "test_02_terg.filtered.sam"

        c = BAMExtract(input_file, False)
        c.extract("chr7:151000000-153000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_02_d(self):
        input_file = TEST_DIR + "test_02_terg.bam"
        output_file = T_TEST_DIR + "test_02_terg.filtered.bam"
        output_file_s = T_TEST_DIR + "test_02_terg.filtered.sam"
        test_file = TEST_DIR + "test_02_terg.filtered.sam"

        c = BAMExtract(input_file, False)
        c.extract("chr5:1-2", "chr7:151000000-153000000", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_02_e(self):
        input_file = TEST_DIR + "test_02_terg.bam"
        output_file = T_TEST_DIR + "test_02_terg.filtered.bam"
        output_file_s = T_TEST_DIR + "test_02_terg.filtered.sam"

        c = BAMExtract(input_file, False)
        c.extract("chr12:151000000-153000000", "chr5:1-2", output_file)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        file_created = False
        with open(output_file_s, "r") as fh:
            file_created = True
            self.assertEqual(fh.read(), "")  # empty file check
        self.assertEqual(file_created, True)  # check presence of fily anyway

if __name__ == '__main__':
    main()
