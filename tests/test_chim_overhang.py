#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Dr. Disco - testing whether the appropriate output is included in `dr-disco detect` and whether `dr-disco classify` indeed classifies based on imbalanced chim_overhang

[License: GNU General Public License v3 (GPLv3)]

 This file is part of Dr. Disco.

 FuMa is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.leker

 Dr. Disco is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
"""


import unittest
import subprocess
import filecmp
import os

from drdisco.IntronDecomposition import IntronDecomposition
from drdisco.DetectOutput import DetectOutput
from drdisco.Classify import Blacklist
from tests.utils import main, sam_to_fixed_bam, get_diff


TEST_DIR = "tests/chim_overhang/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestChimOverhang(unittest.TestCase):
    """super softs:

    test1
    chr2	17281785	17281796	+	chr11	5566699	5566710	-	recurrent in normals & typical 'super soft' (94S21M11S)	024,035
    [9464:521]	chr2:17281790->HBG2,OR52T1P	chr2	17281790	+	0	22	chr11	5566704	-	22	0	inf	valid	linear	intronic	33	22	11	01530	10	1	1	1	0	0	0.7372	0.7372	0.0000	0.0000	0.3818	18.0000	0.8367	0.0013	0.0833	2.0455	38.7727	0.9652	0.0000	0.1847	0.0000	0.3333	2.0000				chr2:17281790/17281791(+)->chr11:5566704/5566705(-):(spanning_paired_1:11,spanning_paired_2:11)

    test2
    chr5	105166407	105166418	+	chr7	12213271	12213282	+	recurrent in normals & supersoft	329,141,151
    [18017:431]	AC091987.1,AC099520.1->TMEM106B	chr5	105166412	+	0	22	chr7	12213276	+	22	0	inf	valid	linear	intronic	33	22	11	0	1528	14	1	1	1	0	0	0.7372	0.7372	0.0000	0.0000	0.1273	15.0000	0.8367	0.0013	0.0278	0.6091	53.9545	0.9419	0.0000	0.0724	0.0000	0.3333	2.0000		chr5:105166412/105166413(+)->chr7:12213276/12213277(+):(spanning_paired_1:11,spanning_paired_2:11)

    test3
    chr8	86472663	86472674	-	chr17	40685756	40685767	+	recurrent in normals & supersoft	036,004
    [11745:251]	RMDN1,WWP1->chr17:40685761	chr8	86472668	-	6	8	chr17	40685761	+	8	6	inf	valid	linear	intronic	21	8	70	998	6	1	1	1	0	0	0.8277	0.8277	0.6547	0.0000	0.3214	53.7500	0.9186	0.0035	0.0619	0.2143	16.9286	0.8660	0.0117	0.0553	0.0000	0.1905	2.0000			chr8:86472668/86472669(-)->chr17:40685761/40685762(+):(spanning_paired_1:3,spanning_paired_1_t:4,spanning_paired_2:3,spanning_paired_2_t:4)
    """

    def test_01(self):
        test_id = '01'

        unfixed_sam = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"

        drdisco_detect = T_TEST_DIR + "test_" + test_id + "_detect.out.txt"
        drdisco_detect_test = TEST_DIR + "test_" + test_id + "_detect.out.txt"

        drdisco_classify = T_TEST_DIR + "test_" + test_id + "_classify.out.txt"
        drdisco_classify_test = TEST_DIR + "test_" + test_id + "_classify.out.txt"

        drdisco_integrate = T_TEST_DIR + "test_" + test_id + "_integrate.out.txt"
        drdisco_integrate_test = TEST_DIR + "test_" + test_id + "_integrate.out.txt"

        # Step 01: dr-disco fix (don't check please)
        sam_to_fixed_bam(unfixed_sam, fixed_bam, T_TEST_DIR)

        # Step 02: dr-disco detect (check appropriate values and columns)
        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(drdisco_detect, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(drdisco_detect_test, drdisco_detect), msg=get_diff( drdisco_detect_test , drdisco_detect ))

        # Step 03: dr-disco classify
        cl = DetectOutput(drdisco_detect)
        cl.classify(drdisco_classify, False, Blacklist(), 25, True)

        self.assertTrue(filecmp.cmp(drdisco_classify_test, drdisco_classify), msg=get_diff( drdisco_classify_test , drdisco_classify ))

        # Step 04: dr-disco integrate
        cl = DetectOutput(drdisco_classify)
        cl.integrate(drdisco_integrate, None, None)

        self.assertTrue(filecmp.cmp(drdisco_integrate_test, drdisco_integrate), msg=get_diff( drdisco_integrate_test , drdisco_integrate ))

    def test_02(self):
        test_id = '02'

        unfixed_sam = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"

        drdisco_detect = T_TEST_DIR + "test_" + test_id + "_detect.out.txt"
        drdisco_detect_test = TEST_DIR + "test_" + test_id + "_detect.out.txt"

        drdisco_classify = T_TEST_DIR + "test_" + test_id + "_classify.out.txt"
        drdisco_classify_test = TEST_DIR + "test_" + test_id + "_classify.out.txt"

        drdisco_integrate = T_TEST_DIR + "test_" + test_id + "_integrate.out.txt"
        drdisco_integrate_test = TEST_DIR + "test_" + test_id + "_integrate.out.txt"

        # Step 01: dr-disco fix (don't check please)
        sam_to_fixed_bam(unfixed_sam, fixed_bam, T_TEST_DIR)

        # Step 02: dr-disco detect (check appropriate values and columns)
        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(drdisco_detect, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(drdisco_detect_test, drdisco_detect), msg=get_diff( drdisco_detect_test , drdisco_detect ))

        # Step 03: dr-disco classify
        cl = DetectOutput(drdisco_detect)
        cl.classify(drdisco_classify, False, Blacklist(), 25, True)

        self.assertTrue(filecmp.cmp(drdisco_classify_test, drdisco_classify), msg=get_diff( drdisco_classify_test , drdisco_classify))

        # Step 04: dr-disco integrate
        cl = DetectOutput(drdisco_classify)
        cl.integrate(drdisco_integrate, None, None)

        self.assertTrue(filecmp.cmp(drdisco_integrate_test, drdisco_integrate), msg=get_diff( drdisco_integrate_test , drdisco_integrate))

    def test_03(self):
        test_id = '03'

        unfixed_sam = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"

        drdisco_detect = T_TEST_DIR + "test_" + test_id + "_detect.out.txt"
        drdisco_detect_test = TEST_DIR + "test_" + test_id + "_detect.out.txt"

        drdisco_classify = T_TEST_DIR + "test_" + test_id + "_classify.out.txt"
        drdisco_classify_test = TEST_DIR + "test_" + test_id + "_classify.out.txt"

        drdisco_integrate = T_TEST_DIR + "test_" + test_id + "_integrate.out.txt"
        drdisco_integrate_test = TEST_DIR + "test_" + test_id + "_integrate.out.txt"

        # Step 01: dr-disco fix (don't check please)
        sam_to_fixed_bam(unfixed_sam, fixed_bam, T_TEST_DIR)

        # Step 02: dr-disco detect (check appropriate values and columns)
        ic = IntronDecomposition(fixed_bam)
        ic.decompose(0)

        with open(drdisco_detect, "w") as fh:
            ic.export(fh)

        self.assertTrue(filecmp.cmp(drdisco_detect_test, drdisco_detect), msg=get_diff( drdisco_detect_test , drdisco_detect ))

        # Step 03: dr-disco classify
        cl = DetectOutput(drdisco_detect)
        cl.classify(drdisco_classify, False, Blacklist(), 25, True)

        self.assertTrue(filecmp.cmp(drdisco_classify_test, drdisco_classify), msg=get_diff( drdisco_classify_test, drdisco_classify ))

        # Step 04: dr-disco integrate
        cl = DetectOutput(drdisco_classify)
        cl.integrate(drdisco_integrate, None, None)

        self.assertTrue(filecmp.cmp(drdisco_integrate_test, drdisco_integrate), msg=get_diff( drdisco_integrate_test , drdisco_integrate ))


if __name__ == '__main__':
    main()
