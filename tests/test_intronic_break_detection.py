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
        
        input_file_a =    "tests/detect-intronic/test_terg_01.sub_02.filtered.fixed.bam"
        input_file_f =    "tests/detect-intronic/test_terg_01_final-list_candidate-fusion-genes.GRCh37.txt"
        
        bps = FusionCatcher(input_file_f,"")
        for bp in bps:
            ic = IntronDecomposition(bp)
            #ic.annotate_genes(gobj)
            ic.decompose(input_file_a)
            
            #decomposed_bp = bp.decompose(input_file_a)
            #print(decomposed_bp)
            return True


def main():
    unittest.main()

if __name__ == '__main__':
    main()
