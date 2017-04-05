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
import subprocess
import filecmp
import os
import pysam
from drdisco.ChimericAlignment import ChimericAlignment

subprocess.call(["bash", "tests/rm_bai_files.sh"])

TEST_DIR = "tests/detect-intronic/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


def sam_to_fixed_bam(sam, fixed_bam):
    basename, ext = os.path.splitext(os.path.basename(sam))
    bam_file = T_TEST_DIR+'/'+basename+'.bam'
    
    fhq = open(bam_file, "wb")
    fhq.write(pysam.view('-b',sam))
    fhq.close()

    alignment_handle = ChimericAlignment(bam_file)
    return alignment_handle.convert(fixed_bam, T_TEST_DIR)

def bam_diff(f1, f2):
    basename, ext = os.path.splitext(os.path.basename(f1))
    
    f1sorted = T_TEST_DIR + basename + '.f1.sorted.bam'
    f2sorted = T_TEST_DIR + basename + '.f2.sorted.bam'
    
    pysam.sort(f1,'-n', '-o', f1sorted)
    pysam.sort(f2,'-n', '-o', f2sorted)
    
    f1sam = T_TEST_DIR + basename + '.f1.sam'
    f2sam = T_TEST_DIR + basename + '.f2.sam'
    
    fhq = open(f1sam, "w")
    fhq.write(pysam.view('-h',f1sorted))
    fhq.close()
    
    fhq = open(f2sam, "w")
    fhq.write(pysam.view('-h',f2sorted))
    fhq.close()
    
    subprocess.Popen(['sed', '-i' , '-r', 's@(SA:[^\\t]+)\\t(LB:[^\\t]+)\t(RG:[^\\t]+)@\\3\\t\\1\\t\\2@', f2sam], stdout=subprocess.PIPE).stdout.read()
    
    
    return filecmp.cmp(f1sam, f2sam),f1sam, f2sam

class TestIntronicBreakDetection(unittest.TestCase):
    #def test_01(self):
        #input_file_a = TEST_DIR + "test_01.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_01.sam"
        #fixed_bam = T_TEST_DIR + "test_01.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ###
        
        #test_file = TEST_DIR + "test_01.out.dbed"
        #output_file = T_TEST_DIR + "test_01.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #fh = open(output_file, "w")
        #ic.export(fh)
        #fh.close()

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_02(self):
        #input_file_a = TEST_DIR + "test_02.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_02.sam"
        #fixed_bam = T_TEST_DIR + "test_02.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ###
        

        #test_file = TEST_DIR + "test_02.out.dbed"
        #output_file = T_TEST_DIR + "test_02.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_03(self):
        #input_file_a = TEST_DIR + "test_03.bam"
        
        ## Dev stuff // insert
        #sam = TEST_DIR + "test_03.sam"
        #fixed_bam = T_TEST_DIR + "test_03.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ###
        
        #test_file = TEST_DIR + "test_03.out.dbed"
        #output_file = T_TEST_DIR + "test_03.out.dbed"


        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_04(self):
        #input_file_a = TEST_DIR + "test_04.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_04.sam"
        #fixed_bam = T_TEST_DIR + "test_04.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##
        

        #test_file = TEST_DIR + "test_04.out.dbed"
        #output_file = T_TEST_DIR + "test_04.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_05(self):
        #input_file_a = TEST_DIR + "test_05.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_05.sam"
        #fixed_bam = T_TEST_DIR + "test_05.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##
        
        #test_file = TEST_DIR + "test_05.out.dbed"
        #output_file = T_TEST_DIR + "test_05.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_06(self):
        #input_file_a = TEST_DIR + "test_06.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_06.sam"
        #fixed_bam = T_TEST_DIR + "test_06.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #output_file = T_TEST_DIR + "test_06.out.dbed"
        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Not sure what 'true' here is exactly - probably reporting err
        ## by STAR? at least, throw a warning and don't terminate

    #def test_07(self):
        #input_file_a = TEST_DIR + "test_07.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_07.sam"
        #fixed_bam = T_TEST_DIR + "test_07.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #output_file = T_TEST_DIR + "test_07.out.dbed"
        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Not sure what 'true' here is exactly, as long as it does
        ## not throw an exception

    #def test_08_test_inclusion_of_disco_reads(self):
        #input_file_a = TEST_DIR + "test_08.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_08.sam"
        #fixed_bam = T_TEST_DIR + "test_08.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_08.out.dbed"
        #output_file = T_TEST_DIR + "test_08.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_09_insert_size_filter(self):
        #input_file_a = TEST_DIR + "test_09.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_09.sam"
        #fixed_bam = T_TEST_DIR + "test_09.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_09.out.dbed"
        #output_file = T_TEST_DIR + "test_09.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_10(self):
        #input_file_a = TEST_DIR + "test_10.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_10.sam"
        #fixed_bam = T_TEST_DIR + "test_10.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_10.out.dbed"
        #output_file = T_TEST_DIR + "test_10.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_11(self):
        #input_file_a = TEST_DIR + "test_11.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_11.sam"
        #fixed_bam = T_TEST_DIR + "test_11.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_11.out.dbed"
        #output_file = T_TEST_DIR + "test_11.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_12_extract_subnetworks(self):
        #input_file_a = TEST_DIR + "test_12.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_12.sam"
        #fixed_bam = T_TEST_DIR + "test_12.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_12.out.dbed"
        #output_file = T_TEST_DIR + "test_12.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_13_merge_overlapping_subnetworks(self):
        #input_file_a = TEST_DIR + "test_13.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_13.sam"
        #fixed_bam = T_TEST_DIR + "test_13.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_13.out.dbed"
        #output_file = T_TEST_DIR + "test_13.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_14_test_inserting_spanning_paired_12_s(self):
        #input_file_a = TEST_DIR + "test_14.bam"
        
        ## Dev stuff // insert
        #sam = TEST_DIR + "test_14.sam"
        #fixed_bam = T_TEST_DIR + "test_14.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##
        
        #test_file = TEST_DIR + "test_14.out.dbed"
        #output_file = T_TEST_DIR + "test_14.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_15_pruning_offset(self):
        #input_file_a = TEST_DIR + "test_15.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_15.sam"
        #fixed_bam = T_TEST_DIR + "test_15.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_15.out.dbed"
        #output_file = T_TEST_DIR + "test_15.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_16_final(self):
        #input_file_a = TEST_DIR + "test_16.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_16.sam"
        #fixed_bam = T_TEST_DIR + "test_16.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_16.out.dbed"
        #output_file = T_TEST_DIR + "test_16.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_17(self):
        #input_file_a = TEST_DIR + "test_17.bam"

        ## Dev stuff // insert
        #sam = TEST_DIR + "test_17.sam"
        #fixed_bam = T_TEST_DIR + "test_17.fixed.bam"
        
        #sam_to_fixed_bam(sam, fixed_bam)
        
        #d = bam_diff(fixed_bam, input_file_a)
        #print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        #input_file_a = fixed_bam
        ##

        #test_file = TEST_DIR + "test_17.out.dbed"
        #output_file = T_TEST_DIR + "test_17.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    def test_18_s27(self):
        input_file_a = TEST_DIR + "test_18.bam"

        # Dev stuff // insert
        sam = TEST_DIR + "test_18.sam"
        fixed_bam = T_TEST_DIR + "test_18.fixed.bam"
        
        sam_to_fixed_bam(sam, fixed_bam)
        
        d = bam_diff(fixed_bam, input_file_a)
        print subprocess.Popen(['diff', d[1], d[2]], stdout=subprocess.PIPE).stdout.read()
        
        input_file_a = fixed_bam
        #

        #test_file = TEST_DIR + "test_18.out.dbed"
        #output_file = T_TEST_DIR + "test_18.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_19_tests_parsing_of_inversed_TERG_from_s55(self):
        #input_file_a = TEST_DIR + "test_19.bam"
        #test_file = TEST_DIR + "test_19.out.dbed"
        #output_file = T_TEST_DIR + "test_19.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_20_tests_trigger_error_on_non_fixed_file(self):
        #input_file_a = TEST_DIR + "test_20.bam"

        #ic = IntronDecomposition(input_file_a)
        #self.assertRaises(Exception, ic.decompose, 0)  # ic.decompose(0) triggers exception

    #def test_21_tests_extracting_subnetworks_in_ideal_optimization_usecase(self):
        #input_file_a = TEST_DIR + "test_21.bam"
        #test_file = TEST_DIR + "test_21.out.dbed"
        #output_file = T_TEST_DIR + "test_21.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_22_interchromosomal_tp_sample_130(self):
        #input_file_a = TEST_DIR + "test_22.bam"
        #test_file = TEST_DIR + "test_22.out.dbed"
        #output_file = T_TEST_DIR + "test_22.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_23_circRNA(self):
        #input_file_a = TEST_DIR + "test_23.bam"
        #test_file = TEST_DIR + "test_23.out.dbed"
        #output_file = T_TEST_DIR + "test_23.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())

    #def test_24_spanning_singleton_x_r(self):
        ## Positions are 100% correct, strands I don't know
        #input_file_a = TEST_DIR + "test_24.bam"
        #test_file = TEST_DIR + "test_24.out.dbed"
        #output_file = T_TEST_DIR + "test_24.out.dbed"

        #ic = IntronDecomposition(input_file_a)
        #ic.decompose(0)

        #with open(output_file, "w") as fh:
            #ic.export(fh)

        ## Test data not checked, should just not throw an exception
        #self.assertTrue(filecmp.cmp(test_file, output_file), msg="diff '" + test_file + "' '" + output_file + "':\n" + subprocess.Popen(['diff', test_file, output_file], stdout=subprocess.PIPE).stdout.read())


def main():
    unittest.main()


if __name__ == '__main__':
    main()
