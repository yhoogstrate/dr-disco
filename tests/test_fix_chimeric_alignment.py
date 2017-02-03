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


import drdisco
from drdisco.ChimericAlignment import ChimericAlignment

import unittest
import logging
import sys
import subprocess
import filecmp
import pysam
import os

logging.basicConfig(level=logging.DEBUG, format=drdisco.__log_format__, stream=sys.stdout)

subprocess.call(["bash", "tests/rm_bai_files.sh"])


class TestChimericAlignment(unittest.TestCase):
    def test_01(self):
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = "tests/fix-chimeric/test_terg_01.filtered.bam"
        output_file = "tmp/test_terg_01.filtered.fixed.bam"
        output_file_s = "tmp/test_terg_01.filtered.fixed.sam"

        test_file = "tests/fix-chimeric/test_terg_01.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg="diff '" + test_file + "' '" + output_file_s + "':\n" + subprocess.Popen(['diff', test_file, output_file_s], stdout=subprocess.PIPE).stdout.read())

    def test_02(self):
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = "tests/fix-chimeric/test_terg_02.filtered.bam"
        output_file = "tmp/test_terg_02.filtered.fixed.bam"
        output_file_s = "tmp/test_terg_02.filtered.fixed.sam"

        test_file = "tests/fix-chimeric/test_terg_02.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))

    def test_03(self):
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = "tests/fix-chimeric/test_terg_03.filtered.bam"
        test_file = "tests/fix-chimeric/test_terg_03.filtered.fixed.sam"

        output_file = "tmp/test_terg_03.filtered.fixed.bam"
        output_file_s = "tmp/test_terg_03.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg="diff '" + test_file + "' '" + output_file_s + "':\n" + subprocess.Popen(['diff', test_file, output_file_s], stdout=subprocess.PIPE).stdout.read())

    def test_04(self):
        # It used to be problematic if 2 mates have exactly the same SA tag (chr, pos, cigar)
        if not os.path.exists("tmp"):
            os.mkdir("tmp")

        input_file = "tests/fix-chimeric/test_terg_04.filtered.bam"
        test_file = "tests/fix-chimeric/test_terg_04.filtered.fixed.sam"

        output_file = "tmp/test_terg_04.filtered.fixed.bam"
        output_file_s = "tmp/test_terg_04.filtered.fixed.sam"

        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file, "tmp")

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg="diff '" + test_file + "' '" + output_file_s + "':\n" + subprocess.Popen(['diff', test_file, output_file_s], stdout=subprocess.PIPE).stdout.read())


def main():
    unittest.main()


if __name__ == '__main__':
    main()
