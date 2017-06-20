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


from drdisco.DetectOutput import DetectOutput


TEST_DIR = "tests/integrate/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_02_s041(self):
        test_id = 'terg_s041'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        gtf_file = None

        cl = DetectOutput(input_file)
        cl.integrate(gtf_file, output_file)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


def main():
    unittest.main()


if __name__ == '__main__':
    main()
