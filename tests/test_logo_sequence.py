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

import drdisco
import logging
import sys
import os
import unittest
import subprocess
import filecmp

logging.basicConfig(level=logging.DEBUG, format=drdisco.__log_format__, stream=sys.stdout)

subprocess.call(["bash", "tests/rm_bai_files.sh"])

TEST_DIR = "tests/logo-sequence/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestLogoSequence(unittest.TestCase):
    def test_01(self):
        input_file = TEST_DIR + "test_01.ref.fa"
        output_file_n = T_TEST_DIR + "test_01.logo-n.fa"
        output_file_p = T_TEST_DIR + "test_01.logo-p.fa"
        test_file_n = TEST_DIR + "test_01.logo-n.fa"
        test_file_p = TEST_DIR + "test_01.logo-p.fa"

        command = ['dr-disco', 'logo-sequence', 'chr1:25', input_file, '-n', '5', '-p', '5', output_file_n, output_file_p]
        print " ".join(command)
        self.assertEqual(subprocess.call(command), 0)  # Ensure error code is 0 - no exceptions have been thrown

        if not filecmp.cmp(output_file_n, test_file_n):
            print 'diff \'' + output_file_n + '\' \'' + test_file_n + '\''

        if not filecmp.cmp(output_file_p, test_file_p):
            print 'diff \'' + output_file_p + '\' \'' + test_file_p + '\''

        self.assertTrue(filecmp.cmp(output_file_n, test_file_n))
        self.assertTrue(filecmp.cmp(output_file_p, test_file_p))

    def test_02(self):
        input_file = TEST_DIR + "test_01.ref.fa"
        output_file_n = T_TEST_DIR + "test_02.logo-n.fa"
        output_file_p = T_TEST_DIR + "test_02.logo-p.fa"
        test_file_n = TEST_DIR + "test_02.logo-n.fa"
        test_file_p = TEST_DIR + "test_02.logo-p.fa"

        command = ['dr-disco', 'logo-sequence', 'chr2:3', input_file, '-n', '4', '-p', '4', output_file_n, output_file_p]
        print " ".join(command)
        self.assertEqual(subprocess.call(command), 0)  # Ensure error code is 0 - no exceptions have been thrown

        if not filecmp.cmp(output_file_n, test_file_n):
            print 'diff \'' + output_file_n + '\' \'' + test_file_n + '\''

        if not filecmp.cmp(output_file_p, test_file_p):
            print 'diff \'' + output_file_p + '\' \'' + test_file_p + '\''

        self.assertTrue(filecmp.cmp(output_file_n, test_file_n))
        self.assertTrue(filecmp.cmp(output_file_p, test_file_p))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
