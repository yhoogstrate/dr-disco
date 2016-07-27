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

import unittest,logging,sys,subprocess,filecmp,pysam
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.Readers import ReadFusionCatcherFinalList as FusionCatcher
from drdisco.IntronDecomposition import IntronDecomposition


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01(self):
        print("\n")
        
        input_file_a =    "tests/detect-intronic/test_terg_01.sub_01.filtered.fixed.bam"
        input_file_f =    "tests/detect-intronic/test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        candidates = ic.decompose(input_file_a)
        
        self.assertEqual(str(candidates[0][0]), "chr21:39877811/39877812(+)->chr21:42873374/42873375(-):(spanning_paired_1:3)")
        self.assertEqual(candidates[0][1], 1.0)


    def test_02(self):
        print("\n")
        
        input_file_a =    "tests/detect-intronic/test_terg_01.sub_02.filtered.fixed.bam"
        input_file_f =    "tests/detect-intronic/test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        candidates = ic.decompose(input_file_a)
        
        self.assertEqual(str(candidates[0][0]), "chr21:39877811/39877812(+)->chr21:42873374/42873375(-):(discordant_mates:3,spanning_paired_1:4)")
        self.assertEqual(candidates[0][1], 1.0)


    def test_03(self):
        print("\n")
        
        input_file_a =    "tests/detect-intronic/test_terg_01.sub_03.filtered.fixed.bam"
        input_file_f =    "tests/detect-intronic/test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        candidates = ic.decompose(input_file_a)
        self.assertEqual(str(candidates[0][0]), 'chr21:39817544/39817545(-)->chr21:42880007/42880008(+):(spanning_paired_1:5)')


    def test_04(self):
        print("\n")
        
        input_file_a =    "tests/detect-intronic/test_terg_01.sub_04.filtered.fixed.bam"
        input_file_f =    "tests/detect-intronic/test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        
        bps = FusionCatcher(input_file_f,"")
        bps_i = bps.__iter__()
        bp = bps_i.next()
        
        ic = IntronDecomposition(bp)
        #ic.annotate_genes(gobj)
        candidates = ic.decompose(input_file_a)
        #self.assertEqual(str(candidates[0][0]), '')
            




def main():
    unittest.main()

if __name__ == '__main__':
    main()
