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


import unittest
import filecmp
import os
import subprocess

from drdisco.Classify import Blacklist
from drdisco.DetectOutput import DetectOutput

from utils import *


D_TEST_DIR = "tests/detect-intronic/"
TEST_DIR = "tests/classify/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIDetectOutputCalssification(unittest.TestCase):
    def test_01(self):
        test_id = '01'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_01__only_valid(self):
        test_id = '01'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_02(self):
        test_id = '02'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_02__only_valid(self):
        test_id = '02'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_03(self):
        test_id = '03'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_03__only_valid(self):
        test_id = '03'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_04(self):
        test_id = '04'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_04__only_valid(self):
        test_id = '04'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_05(self):
        test_id = '05'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_05__only_valid(self):
        test_id = '05'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_08(self):
        test_id = '08'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_08__only_valid(self):
        test_id = '08'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_09(self):
        test_id = '09'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_09__only_valid(self):
        test_id = '09'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_10(self):
        test_id = '10'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_10__only_valid(self):
        test_id = '10'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_11(self):
        test_id = '11'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_11__only_valid(self):
        test_id = '11'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_12(self):
        test_id = '12'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_12__only_valid(self):
        test_id = '12'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_13(self):
        test_id = '13'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_13__only_valid(self):
        test_id = '13'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_14(self):
        test_id = '14'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_14__only_valid(self):
        test_id = '14'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_15(self):
        test_id = '15'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_15__only_valid(self):
        test_id = '15'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_16(self):
        test_id = '16'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_16__only_valid(self):
        test_id = '16'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_17(self):
        test_id = '17'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_17__only_valid(self):
        test_id = '17'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_18(self):
        test_id = '18'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_18__only_valid(self):
        test_id = '18'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_19(self):
        test_id = '19'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_19__only_valid(self):
        test_id = '19'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_21(self):
        test_id = '21'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_21__only_valid(self):
        test_id = '21'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_22(self):
        test_id = '22'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_22__only_valid(self):
        test_id = '22'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_23(self):
        test_id = '23'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_23__only_valid(self):
        test_id = '23'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_24(self):
        test_id = '24'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_24__only_valid(self):
        test_id = '24'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_27(self):
        test_id = '27'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, False, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_27__only_valid(self):
        test_id = '27'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        cl = DetectOutput(input_file)
        cl.classify(output_file, True, Blacklist())

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_blacklists(self):  # only test if they don't crash - do not test actual output
        test_id = '01'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        blacklists = [Blacklist(), Blacklist(), Blacklist(), Blacklist()]
        blacklists[0].add_junctions_from_file('share/blacklist-junctions.hg19.txt')
        blacklists[1].add_junctions_from_file('share/blacklist-junctions.hg38.txt')
        blacklists[2].add_regions_from_bed('share/blacklist-regions.hg19.bed')
        blacklists[3].add_regions_from_bed('share/blacklist-regions.hg38.bed')

        for blacklist in blacklists:
            cl = DetectOutput(input_file)
            cl.classify(output_file, False, blacklist)


if __name__ == '__main__':
    main()
