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


from drdisco.DetectFrameShifts import DetectFrameShifts


TEST_DIR = "tests/integrate/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestFrameShiftPrediction(unittest.TestCase):
    def test_01(self):# example of in-frame fusion - strands are RNA strand
        fusion = ('chr1',1035203,'+'),('chr1',999610,'-')# (from), (to)  and strands are at RNA level!
        gtf_file = TEST_DIR + 'frameshift_example.gtf'
        
        dfs = DetectFrameShifts(gtf_file)
        frameshift_annotation = dfs.evaluate(fusion[0],fusion[1])
        self.assertEqual(str(frameshift_annotation[0]), "[(('AGRN(ENST00000620552.4)-ensembl', 0), ('HES4(ENST00000304952.10)-ensembl_havana', 0))]")
        self.assertEqual(len(frameshift_annotation[1]), 0)
        self.assertEqual(len(frameshift_annotation[2]), 0)

    def test_02(self):# 0, +2
        fusion = ('chr1',1035203,'+'),('chr1',999020,'-')# (from), (to)  and strands are at RNA level!
        gtf_file = TEST_DIR + 'frameshift_example.gtf'
        
        dfs = DetectFrameShifts(gtf_file)
        frameshift_annotation = dfs.evaluate(fusion[0],fusion[1])
        self.assertEqual(len(frameshift_annotation[0]), 0)
        self.assertEqual(len(frameshift_annotation[1]), 0)
        self.assertEqual(str(frameshift_annotation[2]), "[(('AGRN(ENST00000620552.4)-ensembl', 0), ('HES4(ENST00000304952.10)-ensembl_havana', 2))]")

    def test_03(self):# +1, +2 -> 0
        fusion = ('chr1',1040604,'+'),('chr1',999020,'-')
        gtf_file = TEST_DIR + 'frameshift_example.gtf'
        
        dfs = DetectFrameShifts(gtf_file)
        frameshift_annotation = dfs.evaluate(fusion[0],fusion[1])
        self.assertEqual(str(frameshift_annotation[0]), "[(('AGRN(ENST00000620552.4)-ensembl', 1), ('HES4(ENST00000304952.10)-ensembl_havana', 2))]")
        self.assertEqual(len(frameshift_annotation[1]), 0)
        self.assertEqual(len(frameshift_annotation[2]), 0)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
