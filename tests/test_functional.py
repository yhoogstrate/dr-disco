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

import unittest
import filecmp
import pysam
import os
import subprocess
from tests.utils import main, sam_to_fixed_bam, get_diff

# Nosetests doesn't use main()

subprocess.call(["bash", "tests/rm_bai_files.sh"])


class TestFunctional_bam_extract(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/bam-extract/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_bam_extract_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        input_file = TEST_DIR + "test_terg_02.bam"
        output_file = T_TEST_DIR + "test_terg_02.filtered.bam"
        output_file_s = T_TEST_DIR + "test_terg_02.filtered.sam"
        test_file = TEST_DIR + "test_terg_02.filtered.sam"

        # c = BAMExtract(input_file)
        # c.extract("chr21:39000000-40000000", "chr5:1-2", output_file)
        command = ["bin/dr-disco",
                   "bam-extract",
                   "chr21:39000000-40000000",
                   "chr5:1-2",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

        #if not filecmp.cmp(output_file_s, test_file):
        #    print 'diff \'' + output_file_s + '\' \'' + test_file + '\''
        #self.assertTrue(filecmp.cmp(test_file, output_file_s))
        
        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))


class TestFunctional_fix_chimeric(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/fix-chimeric/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_fix_chimeric_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        input_file = TEST_DIR + "test_terg_01.filtered.bam"
        output_file = T_TEST_DIR + "test_terg_01.filtered.fixed.bam"
        output_file_s = T_TEST_DIR + "test_terg_01.filtered.fixed.sam"
        test_file = TEST_DIR + "test_terg_01.filtered.fixed.sam"

        command = ["bin/dr-disco",
                   "fix",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)

        # Bam2Sam
        fhq = open(output_file_s, "w")
        fhq.write(pysam.view(output_file))
        fhq.close()

#        if not filecmp.cmp(output_file_s, test_file):
#            print 'diff \'' + output_file_s + '\' \'' + test_file + '\''

        self.assertTrue(filecmp.cmp(test_file, output_file_s), msg=get_diff(test_file, output_file_s))


class TestFunctional_detect(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/detect-intronic/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_detect_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = '01'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   fixed_bam,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_detect_02(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = '02'

        input_file_a = TEST_DIR + "test_" + test_id + ".sam"
        fixed_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        sam_to_fixed_bam(input_file_a, fixed_bam, T_TEST_DIR)

        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   fixed_bam,
                   output_file]

        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))


class TestFunctional_classify(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/classify/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_classify_16(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()
        test_id = '16'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"

        command = ["bin/dr-disco",
                   "classify",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_classify_16__only_valid(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()
        test_id = '16'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.only-valid.dbed"

        command = ["bin/dr-disco",
                   "classify",
                   "--only-valid",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))


class TestFunctional_integrate(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/integrate/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_02_s041(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'terg_s041'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        gtf_file = TEST_DIR + "test_" + test_id + ".in.gtf"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        command = ["bin/dr-disco",
                   "integrate",
                   "--gtf", gtf_file,
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))
        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_02_s041_no_gtf(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'terg_s041_b'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        command = ["bin/dr-disco",
                   "integrate",
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))


class TestFrameShiftPrediction(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/integrate/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_01(self):  # example of in-frame fusion - strands are RNA strand
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'frameshift-prediction_01'

        # both do have their DNA strand at minus!! :
        #
        #             <=(-)=| acceptor in negative strand at RNA
        #     =====(+)=====>| donor in positive strand at RNA

        #           donor                   acceptor
        # fusions = chr1', 1035203, '+'], ['chr1', 999610, '-'])
        #           1', 1035203, '+'], ['1', 999610, '-'])]  # strands are at RNA level, and gene order is DONOR, ACCEPTOR

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']
        for gtf_file in gtf_files:
            command = ["bin/dr-disco",
                       "integrate",
                       "--gtf", gtf_file,
                       input_file,
                       output_file]

            self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

            self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

            # must statisfy:
            # self.assertEqual(str(frameshift_annotation[0]), "[(('AGRN(ENST00000620552.4)-ensembl', 0), ('HES4(ENST00000304952.10)-ensembl_havana', 0))]")
            # self.assertEqual(len(frameshift_annotation[1]), 0)
            # self.assertEqual(len(frameshift_annotation[2]), 0)

    def test_01_complementary(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'frameshift-prediction_01-complementary'

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']
        for gtf_file in gtf_files:
            command = ["bin/dr-disco",
                       "integrate",
                       "--gtf", gtf_file,
                       input_file,
                       output_file]

            self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

            self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

    def test_02(self):  # 0, +2
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'frameshift-prediction_02'

        # fusions = [(['chr1', 1035203, '+'], ['chr1', 999020, '-']), (['1', 1035203, '+'], ['1', 999020, '-'])]  # (from), (to)  and strands are at RNA level!
        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']
        for gtf_file in gtf_files:
            command = ["bin/dr-disco",
                       "integrate",
                       "--gtf", gtf_file,
                       input_file,
                       output_file]

            self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

            self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

            # must statisfy:
            # self.assertEqual(len(frameshift_annotation[0]), 0)
            # self.assertEqual(len(frameshift_annotation[1]), 0)
            # self.assertEqual(str(frameshift_annotation[2]), "[(('AGRN(ENST00000620552.4)-ensembl', 0), ('HES4(ENST00000304952.10)-ensembl_havana', 2))]")

    def test_03(self):  # +1, +2 -> 0
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'frameshift-prediction_03'

        # fusions = [(['chr1', 1040604, '+'], ['chr1', 999020, '-']), (['1', 1040604, '+'], ['1', 999020, '-'])]

        input_file = TEST_DIR + "test_" + test_id + ".in.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.txt"
        output_file = T_TEST_DIR + "test_" + test_id + ".out.txt"

        gtf_files = [TEST_DIR + 'frameshift_example.gtf', TEST_DIR + 'frameshift_example.no_chr_prefix.gtf']
        for gtf_file in gtf_files:
            command = ["bin/dr-disco",
                       "integrate",
                       "--gtf", gtf_file,
                       input_file,
                       output_file]

            self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

            self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))

            # must statisfy:
            # self.assertEqual(str(frameshift_annotation[0]), "[(('AGRN(ENST00000620552.4)-ensembl', 1), ('HES4(ENST00000304952.10)-ensembl_havana', 2))]")
            # self.assertEqual(len(frameshift_annotation[1]), 0)
            # self.assertEqual(len(frameshift_annotation[2]), 0)


class TestFunctional_integrate_splice_site_motif(unittest.TestCase):
    def __get_temp_dirs(self):
        TEST_DIR = "tests/splice_site_motif/"
        T_TEST_DIR = "tmp/" + TEST_DIR

        if not os.path.exists(T_TEST_DIR):
            os.makedirs(T_TEST_DIR)

        return TEST_DIR, T_TEST_DIR

    def test_sj_01(self):
        TEST_DIR, T_TEST_DIR = self.__get_temp_dirs()

        test_id = 'splice_site_motif_01'

        input_sam = TEST_DIR + "test_" + test_id + ".in.sam"
        input_bam = T_TEST_DIR + "test_" + test_id + ".fixed.bam"
        input_file = T_TEST_DIR + "test_" + test_id + ".dbed"

        # gtf_file = None
        fasta_file = TEST_DIR + "test_" + test_id + ".in.fa"

        output_file = T_TEST_DIR + "test_" + test_id + ".out.dbed"
        test_file = TEST_DIR + "test_" + test_id + ".out.dbed"

        # sam -> fixed bam
        command = ["bin/dr-disco",
                   "fix",
                   input_sam,
                   input_bam]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

        # fixed bam -> dr-disco detect
        command = ["bin/dr-disco",
                   "detect",
                   "-m", "0",
                   input_bam,
                   input_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

        # dr-disco-detect (skip classify) -> dr-disco integrate
        command = ["bin/dr-disco",
                   "integrate",
                   "--fasta", fasta_file,
                   input_file,
                   output_file]

        self.assertEqual(subprocess.call(command), 0, msg=" ".join([str(x) for x in command]))

        self.assertTrue(filecmp.cmp(test_file, output_file), msg=get_diff(test_file, output_file))


if __name__ == '__main__':
    main()
