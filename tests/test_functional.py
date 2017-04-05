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

import unittest
import filecmp
import pysam
import os
import subprocess

from test_intronic_break_detection import sam_to_fixed_bam

# Nosetests doesn't use main()

subprocess.call(["bash", "tests/rm_bai_files.sh"])


class TestFunctional_bam_extract(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/bam-extract/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_bam_extract_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"
        test_file = TEST_DIR + "test_terg_02.filtered.sam"

        # c = BAMExtract(input_file)
        # c.extract("chr21:39000000-40000000", "chr5:1-2", output_file)
        command = ["bin/dr-disco",
                   "bam-extract",
                   "chr21:39000000-40000000",
                   "chr5:1-2",
                   output_file,
                   input_file]

        self.assertEqual(subprocess.call(command), 0)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))


class TestFunctional_detect(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/detect-intronic/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_detect_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()
        input_file_s = TEST_DIR + "test_01.sam"
        input_file_a = T_TEST_DIR + "test_01.bam"

        sam_to_fixed_bam(input_file_s, input_file_a)
        test_file = TEST_DIR + "test_01.out.dbed"
        output_file = T_TEST_DIR + "test_01.out.dbed"

        # ic = IntronDecomposition(input_file_a)
        # ic.decompose(0)
        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   output_file,
                   input_file_a]

        self.assertEqual(subprocess.call(command), 0)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_detect_02(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()
        input_file_s = TEST_DIR + "test_02.sam"
        input_file_a = T_TEST_DIR + "test_02.bam"

        sam_to_fixed_bam(input_file_s, input_file_a)

        test_file = TEST_DIR + "test_02.out.dbed"
        output_file = T_TEST_DIR + "test_02.out.dbed"

        # ic = IntronDecomposition(input_file_a)
        # ic.decompose(0)
        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   output_file,
                   input_file_a]

        self.assertEqual(subprocess.call(command), 0)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


def main():
    unittest.main()


if __name__ == '__main__':
    main()
