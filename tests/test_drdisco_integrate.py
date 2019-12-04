#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Dr. Disco - testing integrate

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
import subprocess
from drdisco.DetectOutput import DetectOutput
from tests.utils import main, get_diff


TEST_DIR = "tests/integrate/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_s041(self):
        test_id = 'terg_s041'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        gtf_files = [TEST_DIR + "test_" + test_id + ".in.gtf", TEST_DIR + "test_" + test_id + ".in.no_chr_prefix.gtf"]
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        for gtf_file in gtf_files:
            cl = DetectOutput(input_file)
            cl.integrate(output_file, gtf_file, None)

            self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))
            #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_s041_nocrash(self):
        test_id = 'terg_s041_b'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        gtf_file = TEST_DIR + "example_refseq.gff"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        cl = DetectOutput(input_file)
        cl.integrate(output_file, gtf_file, None)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_s041_no_gtf(self):
        test_id = 'terg_s041_b'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        gtf_file = None
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        cl = DetectOutput(input_file)
        cl.integrate(output_file, gtf_file, None)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_in_frame_non_hybrid_protein(self):
        test_id = 'in_frame_non_hybrid_protein'
        # Transcript ID's necessary:
        # - TMPRSS2: ENST00000424093
        # - ERG: ENST00000398910

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        gtf_files = [TEST_DIR + "test_" + test_id + ".gtf"]
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        for gtf_file in gtf_files:
            cl = DetectOutput(input_file)
            cl.integrate(output_file, gtf_file, None)

            #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())
            self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))


if __name__ == '__main__':
    main()
