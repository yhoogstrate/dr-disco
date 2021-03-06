#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Dr. Disco - testing fix-chimeric

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


from drdisco.ChimericAlignment import ChimericAlignment

import unittest
import subprocess
import filecmp
import pysam
import os
from tests.utils import main, get_diff


subprocess.call(["bash", "tests/rm_bai_files.sh"])


TEST_DIR = "tests/fix-chimeric/"
T_TEST_DIR = "tmp/" + TEST_DIR

if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestChimericAlignment(unittest.TestCase):
    def test_01(self):
        input_file = TEST_DIR + "test_terg_01.filtered.bam"
        output_file = T_TEST_DIR + "test_terg_01.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_terg_01.filtered.fixed.sam"

        test_file = TEST_DIR + "test_terg_01.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        #self.assertTrue(filecmp.cmp(test_file, output_file_s), msg="diff '" + test_file + "' '" + output_file_s + "':\n" + subprocess.Popen(['diff', test_file, output_file_s], stdout=subprocess.PIPE).stdout.read())
        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_02(self):
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = TEST_DIR + "test_terg_02.filtered.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.fixed.sam"

        test_file = TEST_DIR + "test_terg_02.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_03(self):
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = TEST_DIR + "test_terg_03.filtered.bam"
        test_file = TEST_DIR + "test_terg_03.filtered.fixed.sam"

        output_file = T_TEST_DIR + "test_terg_03.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_terg_03.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        #self.assertTrue(filecmp.cmp(test_file, output_file_s), msg="diff '" + test_file + "' '" + output_file_s + "':\n" + subprocess.Popen(['diff', test_file, output_file_s], stdout=subprocess.PIPE).stdout.read())
        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_04(self):
        # It used to be problematic if 2 mates have exactly the same SA tag (chr, pos, cigar)
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = TEST_DIR + "test_terg_04.filtered.bam"
        test_file = TEST_DIR + "test_terg_04.filtered.fixed.sam"

        output_file = T_TEST_DIR + "test_terg_04.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_terg_04.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        #self.assertTrue(filecmp.cmp(test_file, output_file_s), msg="diff '" + test_file + "' '" + output_file_s + "':\n" + subprocess.Popen(['diff', test_file, output_file_s], stdout=subprocess.PIPE).stdout.read())
        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))

    def test_05(self):
        # It used to be problematic if 2 mates have exactly the same SA tag (chr, pos, cigar)
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = TEST_DIR + "test_05.sam"

        output_file = T_TEST_DIR + "test_05.fixed.bam"
        alignment_handle = ChimericAlignment(input_file)

        # alignment_handle.convert(output_file, "tmp")
        self.assertRaises(Exception, alignment_handle.convert, output_file, "tmp")  # triggers exception because of empty file


if __name__ == '__main__':
    main()
