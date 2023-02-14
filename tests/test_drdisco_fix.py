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
import pysam


subprocess.call(["bash", "tests/rm_bai_files.sh"])


TEST_DIR = "tests/fix-chimeric/"
T_TEST_DIR = "tmp/" + TEST_DIR

if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestChimericAlignment(unittest.TestCase):
    def test_01(self):
        input_file = TEST_DIR + "test_01_terg.filtered.bam"
        output_file = T_TEST_DIR + "test_01_terg.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_01_terg.filtered.fixed.sam"

        test_file = TEST_DIR + "test_01_terg.filtered.fixed.sam"

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

        input_file = TEST_DIR + "test_02_terg.filtered.bam"
        output_file = T_TEST_DIR + "test_02_terg.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_02_terg.filtered.fixed.sam"

        test_file = TEST_DIR + "test_02_terg.filtered.fixed.sam"

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

        input_file = TEST_DIR + "test_03_terg.filtered.bam"
        test_file = TEST_DIR + "test_03_terg.filtered.fixed.sam"

        output_file = T_TEST_DIR + "test_03_terg.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_03_terg.filtered.fixed.sam"

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

        input_file = TEST_DIR + "test_04_terg.filtered.bam"
        test_file = TEST_DIR + "test_04_terg.filtered.fixed.sam"

        output_file = T_TEST_DIR + "test_04_terg.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_04_terg.filtered.fixed.sam"

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


    def test_06__dealing_with_incorrect_hi_tags_in_within_bam_alignments(self):
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file_chim = TEST_DIR + "test_06__Chimeric.out.sam"
        output_file_chim = T_TEST_DIR + "test_06__Chimeric.out.fixed.bam"
        alignment_handle = ChimericAlignment(input_file_chim)
        alignment_handle.convert(output_file_chim, "tmp")

        input_file_within_bam = TEST_DIR + "test_06__Aligned.sortedByCoord.out.sam"
        output_file_within_bam = T_TEST_DIR + "test_06__Aligned.sortedByCoord.out.fixed.bam"
        alignment_handle = ChimericAlignment(input_file_within_bam)
        alignment_handle.convert(output_file_within_bam, "tmp")
        
        fixed_chim = pysam.AlignmentFile(output_file_chim, "rb")
        fixed_within_bam = pysam.AlignmentFile(output_file_within_bam, "rb")
        
        for read1, read2 in zip(fixed_chim.fetch(), fixed_within_bam.fetch()):
            self.assertEqual(read1.query_name, read2.query_name)
            self.assertEqual(read1.pos, read2.pos)
            if read1.get_tag('RG').find("discordant_mates") == -1 and read1.get_tag('RG').find("silent_mate"):
                self.assertEqual(read1.get_tag('HI'), read2.get_tag('HI'))
            self.assertEqual(read1.get_tag('RG'), read2.get_tag('RG'))


if __name__ == '__main__':
    main()
