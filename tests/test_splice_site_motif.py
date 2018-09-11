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


from drdisco.IntronDecomposition import IntronDecomposition

import unittest
import filecmp
import os
import subprocess
from drdisco.DetectOutput import DetectOutput
from utils import main, sam_to_fixed_bam


TEST_DIR = "tests/splice_site_motif/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestSpliceJunctions(unittest.TestCase):
    def test_sj_01(self):
        test_id = 'splice_site_motif_01'

        input_sam = TEST_DIR + "test_" + test_id + ".in.sam"
        input_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        input_file = T_TEST_DIR + "test_" + test_id + ".dbed"

        gtf_file = None
        fasta_file = TEST_DIR + "test_" + test_id + ".in.fa"

        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"

        # sam -> fixed bam
        sam_to_fixed_bam(input_sam, input_bam, T_TEST_DIR)

        # fixed bam -> dr-disco detect
        ic = IntronDecomposition(input_bam)
        ic.decompose(0)
        fh = open(input_file, "w")
        ic.export(fh)
        fh.close()

        # dr-disco-detect (skip classify) -> dr-disco integrate
        cl = DetectOutput(input_file)
        cl.integrate(output_file, gtf_file, fasta_file)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_sj_02(self):
        test_id = 'splice_site_motif_02'

        input_sam = TEST_DIR + "test_" + test_id + ".in.sam"
        input_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        input_file = T_TEST_DIR + "test_" + test_id + ".dbed"

        gtf_file = None
        fasta_file = TEST_DIR + "test_" + test_id + ".in.fa"

        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"

        # sam -> fixed bam
        sam_to_fixed_bam(input_sam, input_bam, T_TEST_DIR)

        # fixed bam -> dr-disco detect
        ic = IntronDecomposition(input_bam)
        ic.decompose(0)
        fh = open(input_file, "w")
        ic.export(fh)
        fh.close()

        # dr-disco-detect (skip classify) -> dr-disco integrate
        cl = DetectOutput(input_file)
        cl.integrate(output_file, gtf_file, fasta_file)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_sj_03__go_out_of_bound_in_the_fasta_file(self):
        test_id = 'splice_site_motif_03'

        input_sam = TEST_DIR + "test_" + test_id + ".in.sam"
        input_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        input_file = T_TEST_DIR + "test_" + test_id + ".dbed"

        gtf_file = None
        fasta_file = TEST_DIR + "test_" + test_id + ".in.fa"

        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"

        # sam -> fixed bam
        sam_to_fixed_bam(input_sam, input_bam, T_TEST_DIR)

        # fixed bam -> dr-disco detect
        ic = IntronDecomposition(input_bam)
        ic.decompose(0)
        fh = open(input_file, "w")
        ic.export(fh)
        fh.close()

        # dr-disco-detect (skip classify) -> dr-disco integrate
        cl = DetectOutput(input_file)

        # originally, this triggered an exception, now we just log an error
        # cl.integrate(output_file, gtf_file, fasta_file)
        # self.assertRaises(Exception, cl.integrate, output_file, gtf_file, fasta_file)

        cl.integrate(output_file, gtf_file, fasta_file)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_sj_04__CA_GT_d4(self):
        test_id = 'splice_site_motif_04'

        input_sam = TEST_DIR + "test_" + test_id + ".in.sam"
        input_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        input_file = T_TEST_DIR + "test_" + test_id + ".dbed"

        gtf_file = None
        fasta_file = TEST_DIR + "test_" + test_id + ".in.fa"

        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"

        # sam -> fixed bam
        sam_to_fixed_bam(input_sam, input_bam, T_TEST_DIR)

        # fixed bam -> dr-disco detect
        ic = IntronDecomposition(input_bam)
        ic.decompose(0)
        fh = open(input_file, "w")
        ic.export(fh)
        fh.close()

        # dr-disco-detect (skip classify) -> dr-disco integrate
        cl = DetectOutput(input_file)

        # originally, this triggered an exception, now we just log an error
        # cl.integrate(output_file, gtf_file, fasta_file)
        # self.assertRaises(Exception, cl.integrate, output_file, gtf_file, fasta_file)

        cl.integrate(output_file, gtf_file, fasta_file)

        self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


if __name__ == '__main__':
    main()
