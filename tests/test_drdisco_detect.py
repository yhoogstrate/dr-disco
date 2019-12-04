#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Dr. Disco - testing `dr-disco detect ...`. This includes running `dr-disco fix ...` first to ensure the whole pipeline is in working.

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

from drdisco.IntronDecomposition import IntronDecomposition

import unittest
import subprocess
import filecmp
import os
from utils import main, sam_to_fixed_bam, get_diff


subprocess.call(["bash", "tests/rm_bai_files.sh"])

TEST_DIR = "tests/detect-intronic/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01(self):
        test_id = '01'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        fh = open(output_file, "w")
        ic.export(fh)
        fh.close()

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_02(self):
        test_id = '02'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_03(self):
        test_id = '03'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_04(self):
        test_id = '04'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_05(self):
        test_id = '05'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_06(self):
        test_id = '06'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Not sure what 'true' here is exactly - probably reporting err
        # by STAR? at least, throw a warning and don't terminate

    def test_07(self):
        test_id = '07'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Not sure what 'true' here is exactly, as long as it does
        # not throw an exception

    def test_08_test_inclusion_of_disco_reads(self):
        test_id = '08'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_09_insert_size_filter(self):
        test_id = '09'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_10(self):
        test_id = '10'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_11(self):
        test_id = '11'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_12_extract_subnetworks(self):
        test_id = '12'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_13_merge_overlapping_subnetworks_s54(self):
        test_id = '13'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_14_test_inserting_spanning_paired_12_s(self):
        test_id = '14'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_15_pruning_offset(self):
        test_id = '15'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_16_final(self):
        test_id = '16'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_17(self):
        test_id = '17'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_18_s27(self):
        test_id = '18'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))
        
        # order is different, but content matching is still valid:
        self.assertTrue( sorted(open(test_file,'r').read().split("\n"))  == sorted(open(output_file,'r').read().split("\n")) )

    def test_19_tests_parsing_of_inversed_TERG_s55(self):
        test_id = '19'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_20_tests_trigger_error_on_non_fixed_file(self):
        input_file_a = TEST_DIR + "test_20.bam"

        ic = IntronDecomposition(input_file_a)
        self.assertRaises(Exception, ic.decompose, 0)  # ic.decompose(0) triggers exception

    def test_21_tests_extracting_subnetworks_in_ideal_optimization_usecase(self):
        test_id = '21'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_22_interchromosomal_tp_sample_130(self):
        test_id = '22'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_23_circRNA_s14(self):
        test_id = '23'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_24_spanning_singleton_x_r(self):
        test_id = '24'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_27_bp_entropy_s54(self):
        test_id = '27'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_28_bp_entropy_s_b6729(self):
        test_id = '28'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_29_mismatch_ratio_and_stddev_b11383(self):
        # sum matches needs to be: 5675 (2919 + 2756 = 5675)
        # sum mismatches needs to be: 149 (73 + 76 = 149)
        test_id = '29'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_30_many_muts_b11745(self):
        test_id = '30'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_31_many_muts_on_test_12_for_testing_manymuts_merging_in_subgraphs(self):
        """
            mm needs to be 4*1=4
            m needs to be 4*50=200
        """
        test_id = '31'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(output_file, "w") as fh:
            ic.export(fh)

        # Test data not checked, should just not throw an exception
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))


if __name__ == '__main__':
    main()
