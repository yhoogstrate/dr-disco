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

import unittest,logging,sys,subprocess,filecmp,pysam,os
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.Readers import ReadFusionCatcherFinalList as FusionCatcher
from drdisco.IntronDecomposition import IntronDecomposition


TEST_DIR = "tests/detect-intronic/"
T_TEST_DIR = "tmp/"+TEST_DIR


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01(self):
        #print("\n")
        
        input_file_a =   TEST_DIR+"test_terg_01.sub_01.filtered.fixed.bam"
        input_file_f =   TEST_DIR+"test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file  =     TEST_DIR+"test_01.out.dbed"
        output_file  = T_TEST_DIR+"test_01.out.dbed"
        
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        n_candidates = ic.decompose(input_file_a)
        
        with open(output_file, "w") as fh:
            ic.export(fh)
        
        self.assertEqual(n_candidates, 1)
        self.assertTrue(filecmp.cmp(test_file, output_file))

    def test_02(self):
        #print("\n")
        
        input_file_a =   TEST_DIR+"test_terg_01.sub_02.filtered.fixed.bam"
        input_file_f =   TEST_DIR+"test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file  =     TEST_DIR+"test_02.out.dbed"
        output_file  = T_TEST_DIR+"test_02.out.dbed"
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        print "vvv"
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        n_candidates = ic.decompose(input_file_a)
        
        with open(output_file, "w") as fh:
            ic.export(fh)
        
        print "^^^"
        
        self.assertEqual(n_candidates, 1)
        self.assertTrue(filecmp.cmp(test_file, output_file))


    def test_03(self):
        #print("\n")
        
        input_file_a =   TEST_DIR+"test_terg_01.sub_03.filtered.fixed.bam"
        input_file_f =   TEST_DIR+"test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file  =     TEST_DIR+"test_03.out.dbed"
        output_file  = T_TEST_DIR+"test_03.out.dbed"
        
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        n_candidates = ic.decompose(input_file_a)
        
        with open(output_file, "w") as fh:
            ic.export(fh)
        
        self.assertEqual(n_candidates, 1)
        self.assertTrue(filecmp.cmp(test_file, output_file))


    def test_04(self):
        #print("\n")
        
        input_file_a =    TEST_DIR+"test_terg_01.sub_04.filtered.fixed.bam"
        input_file_f =    TEST_DIR+"test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file  =      TEST_DIR+"test_04.out.dbed"
        output_file  =  T_TEST_DIR+"test_04.out.dbed"
        
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        n_candidates = ic.decompose(input_file_a)
        
        with open(output_file, "w") as fh:
            ic.export(fh)
        
        self.assertEqual(n_candidates, 2)
        self.assertTrue(filecmp.cmp(test_file, output_file))


    def test_05(self):
        #print("\n")
        
        input_file_a =    TEST_DIR+"test_terg_01.sub_05.filtered.fixed.bam"
        input_file_f =    TEST_DIR+"test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        test_file  =      TEST_DIR+"test_05.out.dbed"
        output_file  =  T_TEST_DIR+"test_05.out.dbed"
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        n_candidates = ic.decompose(input_file_a)
        
        with open(output_file, "w") as fh:
            ic.export(fh)
        
        self.assertEqual(n_candidates, 2)
        self.assertTrue(filecmp.cmp(test_file, output_file))



def main():
    if not os.path.exists(T_TEST_DIR):
        os.makedirs(T_TEST_DIR)
    
    unittest.main()

if __name__ == '__main__':
    main()
