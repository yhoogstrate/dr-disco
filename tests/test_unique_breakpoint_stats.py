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


import os
import unittest

from drdisco.Classify import Blacklist
from drdisco.DetectOutput import DetectOutput
from drdisco.IntronDecomposition import IntronDecomposition

from tests.utils import main, sam_to_fixed_bam


TEST_DIR = "tests/unique-breakpoint-stats/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIDetectOutputCalssification(unittest.TestCase):
    def test_01(self):
        test_id = 'vcap_err_01'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        detect_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.classified.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        fh = open(detect_file, "w")
        ic.export(fh)
        fh.close()

        cl = DetectOutput(detect_file)
        cl.classify(output_file, False, Blacklist(), 1, True)

        n_valid = 0
        with open(output_file) as fh:
            for line in fh:
                if line.find('valid') > -1:
                    n_valid += 1

        self.assertEqual(n_valid, 0)

    def test_02(self):
        test_id = 'vcap_err_02'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        detect_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.classified.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        fh = open(detect_file, "w")
        ic.export(fh)
        fh.close()

        cl = DetectOutput(detect_file)
        cl.classify(output_file, False, Blacklist(), 1, True)

        n_valid = 0
        with open(output_file) as fh:
            for line in fh:
                if line.find('valid') > -1:
                    n_valid += 1

        self.assertEqual(n_valid, 0)


if __name__ == '__main__':
    main()
