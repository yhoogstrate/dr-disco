#!/usr/bin/env python2
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
from drdisco.IntronDecomposition import IntronDecomposition

import unittest
import logging
import sys
import subprocess
import filecmp
import os

logging.basicConfig(level=logging.DEBUG, format=drdisco.__log_format__, stream=sys.stdout)

TEST_DIR = "tests/detect-intronic/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01(self):
        input_file_a = TEST_DIR + "test_01.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_01.out.dbed"
        output_file = T_TEST_DIR + "test_01.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        fh = open(output_file, "w")
        ic.export(fh)
        fh.close()

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_02(self):
        input_file_a = TEST_DIR + "test_02.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_02.out.dbed"
        output_file = T_TEST_DIR + "test_02.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_03(self):
        input_file_a = TEST_DIR + "test_03.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_03.out.dbed"
        output_file = T_TEST_DIR + "test_03.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_04(self):
        input_file_a = TEST_DIR + "test_04.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_04.out.dbed"
        output_file = T_TEST_DIR + "test_04.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_05(self):
        input_file_a = TEST_DIR + "test_05.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_05.out.dbed"
        output_file = T_TEST_DIR + "test_05.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_06(self):
        input_file_a = TEST_DIR + "test_06.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        output_file = T_TEST_DIR + "test_06.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Not sure what 'true' here is exactly - probably reporting err
        # by STAR? at least, throw a warning and don't terminate

    def test_07(self):
        input_file_a = TEST_DIR + "test_07.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        output_file = T_TEST_DIR + "test_07.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Not sure what 'true' here is exactly, as long as it does
        # not throw an exception

    def test_08_test_inclusion_of_disco_reads(self):
        input_file_a = TEST_DIR + "test_08.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_08.out.dbed"
        output_file = T_TEST_DIR + "test_08.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_09(self):
        input_file_a = TEST_DIR + "test_09.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_09.out.dbed"
        output_file = T_TEST_DIR + "test_09.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_10(self):
        input_file_a = TEST_DIR + "test_10.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_10.out.dbed"
        output_file = T_TEST_DIR + "test_10.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_11(self):
        input_file_a = TEST_DIR + "test_11.bam"
        # input_file_f = TEST_DIR + "test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file = TEST_DIR + "test_11.out.dbed"
        output_file = T_TEST_DIR + "test_11.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_12_extract_subnetworks(self):
        input_file_a = TEST_DIR + "test_12.bam"
        test_file = TEST_DIR + "test_12.out.dbed"
        output_file = T_TEST_DIR + "test_12.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_13_merge_overlapping_subnetworks(self):
        input_file_a = TEST_DIR + "test_13.bam"
        test_file = TEST_DIR + "test_13.out.dbed"
        output_file = T_TEST_DIR + "test_13.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_14_test_inserting_spanning_paired_12_s(self):
        input_file_a = TEST_DIR + "test_14.bam"
        test_file = TEST_DIR + "test_14.out.dbed"
        output_file = T_TEST_DIR + "test_14.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_15_pruning_offset(self):
        input_file_a = TEST_DIR + "test_15.bam"
        test_file = TEST_DIR + "test_15.out.dbed"
        output_file = T_TEST_DIR + "test_15.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_16_final(self):
        input_file_a = TEST_DIR + "test_16.bam"
        test_file = TEST_DIR + "test_16.out.dbed"
        output_file = T_TEST_DIR + "test_16.out.dbed"

        ic = IntronDecomposition(input_file_a)
        # ic.annotate_genes(gobj)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_17(self):
        input_file_a = TEST_DIR + "test_17.bam"
        test_file = TEST_DIR + "test_17.out.dbed"
        output_file = T_TEST_DIR + "test_17.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_18(self):
        input_file_a = TEST_DIR + "test_18.bam"
        test_file = TEST_DIR + "test_18.out.dbed"
        output_file = T_TEST_DIR + "test_18.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_19_tests_parsing_of_inversed_TERG_from_s55(self):
        input_file_a = TEST_DIR + "test_19.bam"
        test_file = TEST_DIR + "test_19.out.dbed"
        output_file = T_TEST_DIR + "test_19.out.dbed"

        ic = IntronDecomposition(input_file_a)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


def main():
    unittest.main()


if __name__ == '__main__':
    main()
