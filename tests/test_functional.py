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
from utils import *

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
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        if not filecmp.cmp(output_file_s, test_file):
            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(output_file_s, test_file))


class TestFunctional_fix_chimeric(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/fix-chimeric/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_fix_chimeric_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        input_file = TEST_DIR + "test_terg_01.filtered.bam"
        output_file = T_TEST_DIR + "test_terg_01.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_terg_01.filtered.fixed.sam"
        test_file = TEST_DIR + "test_terg_01.filtered.fixed.sam"

        command = ["bin/dr-disco",
                   "fix",
                   input_file,
                   output_file]

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

        test_id = '01'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   fixed_bam,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_detect_02(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = '02'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   fixed_bam,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


class TestFunctional_classify(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/classify/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_classify_16(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()
        test_id = '16'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        command = ["bin/dr-disco",
                   "classify",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_classify_16__only_valid(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()
        test_id = '16'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        command = ["bin/dr-disco",
                   "classify",
                   "--only-valid",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


class TestFunctional_integrate(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/integrate/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_02_s041(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'terg_s041'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        gtf_file = TEST_DIR + "test_" + test_id + ".in.gtf"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        command = ["bin/dr-disco",
                   "integrate",
                   "--gtf", gtf_file,
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_02_s041_no_gtf(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'terg_s041_b'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        command = ["bin/dr-disco",
                   "integrate",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


if __name__ == '__main__':
    main()
