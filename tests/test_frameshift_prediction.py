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
 (at your option) any later version.leker

 Dr. Disco is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
"""


import unittest
import os

from drdisco.DetectFrameShifts import DetectFrameShifts
from utils import main


TEST_DIR = "tests/integrate/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestFrameShiftPrediction(unittest.TestCase):
    def test_01(self):  # example of in-frame fusion - strands are RNA strand
        fusions = [(['chr1', 1035203, '+'], ['chr1', 999610, '-']), (['1', 1035203, '+'], ['1', 999610, '-'])]  # strands are at RNA level, and gene order is DONOR, ACCEPTOR
        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']

        for fusion in fusions:
            for gtf_file in gtf_files:
                dfs = DetectFrameShifts(gtf_file)
                exons_from, exons_to, frameshift_annotation = dfs.evaluate(fusion[0], fusion[1], 0)
                self.assertEqual(str(frameshift_annotation[0]), "[(('AGRN(ENST00000620552.4)-ensembl', 0), ('HES4(ENST00000304952.10)-ensembl_havana', 0))]")
                self.assertEqual(len(frameshift_annotation[1]), 0)
                self.assertEqual(len(frameshift_annotation[2]), 0)

                self.assertEqual(",".join(exons_from), "AGRN(ENST00000620552.4)-ensembl")
                self.assertEqual(",".join(exons_to), "HES4(ENST00000304952.10)-ensembl_havana")

    def test_02(self):  # 0, +2
        fusions = [(['chr1', 1035203, '+'], ['chr1', 999020, '-']), (['1', 1035203, '+'], ['1', 999020, '-'])]  # (from), (to)  and strands are at RNA level!
        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']

        for fusion in fusions:
            for gtf_file in gtf_files:
                dfs = DetectFrameShifts(gtf_file)
                exons_from, exons_to, frameshift_annotation = dfs.evaluate(fusion[0], fusion[1], 0)
                self.assertEqual(len(frameshift_annotation[0]), 0)
                self.assertEqual(len(frameshift_annotation[1]), 0)
                self.assertEqual(str(frameshift_annotation[2]), "[(('AGRN(ENST00000620552.4)-ensembl', 0), ('HES4(ENST00000304952.10)-ensembl_havana', 2))]")

                self.assertEqual(",".join(exons_from), "AGRN(ENST00000620552.4)-ensembl")
                self.assertEqual(",".join(exons_to), "HES4(ENST00000304952.10)-ensembl_havana")

    def test_03(self):  # +1, +2 -> 0
        fusions = [(['chr1', 1040604, '+'], ['chr1', 999020, '-']), (['1', 1040604, '+'], ['1', 999020, '-'])]
        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']

        for fusion in fusions:
            for gtf_file in gtf_files:
                dfs = DetectFrameShifts(gtf_file)
                exons_from, exons_to, frameshift_annotation = dfs.evaluate(fusion[0], fusion[1], 0)
                self.assertEqual(str(frameshift_annotation[0]), "[(('AGRN(ENST00000620552.4)-ensembl', 1), ('HES4(ENST00000304952.10)-ensembl_havana', 2))]")
                self.assertEqual(len(frameshift_annotation[1]), 0)
                self.assertEqual(len(frameshift_annotation[2]), 0)

                self.assertEqual(",".join(exons_from), "AGRN(ENST00000620552.4)-ensembl")
                self.assertEqual(",".join(exons_to), "HES4(ENST00000304952.10)-ensembl_havana")


if __name__ == '__main__':
    main()
