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
from drdisco.CigarAlignment import *


class TestIntronicBreakDetection(unittest.TestCase):
    def test_01(self):
        cig1 = "15S15M"
        cig2 = "15M15S"
        
        ans1 = [(4, 15), (0, 15)]
        ans2 = [(0, 15), (4, 15)]
        
        ca = CigarAlignment(cigar_to_cigartuple(cig1),cigar_to_cigartuple(cig2))
        ca.fill_matrix()
        aligned_cigars = ca.traceback_matrix()
        
        #print aligned_cigars[0]
        #print aligned_cigars[1]
        
        self.assertEqual(aligned_cigars[0],ans1)
        self.assertEqual(aligned_cigars[1],ans2)
        
    def test_02(self):
        cig1 = "15M15S"
        cig2 = "15S15M"
        
        ans1 = [(0, 15), (4, 15)]
        ans2 = [(4, 15), (0, 15)]
        
        ca = CigarAlignment(cigar_to_cigartuple(cig1),cigar_to_cigartuple(cig2))
        ca.fill_matrix()
        aligned_cigars = ca.traceback_matrix()
        
        #print aligned_cigars[0]
        #print aligned_cigars[1]
        
        self.assertEqual(aligned_cigars[0],ans1)
        self.assertEqual(aligned_cigars[1],ans2)

    def test_03(self):
        cig1 = "16M28431N69M41S"
        cig2 = "85S41M"
        
        ans1 = [(0, 85), (4, 41)]
        ans2 = [(4, 85), (0, 41)]
        
        ca = CigarAlignment(cigar_to_cigartuple(cig1),cigar_to_cigartuple(cig2))
        ca.fill_matrix()
        aligned_cigars = ca.traceback_matrix()
        
        #print aligned_cigars[0]
        #print aligned_cigars[1]
        
        self.assertEqual(aligned_cigars[0],ans1)
        self.assertEqual(aligned_cigars[1],ans2)
    
    def test_04(self):
        cig1 = "85S41M"
        cig2 = "16M28431N69M41S"
        
        ans1 = [(4, 85), (0, 41)]
        ans2 = [(0, 85), (4, 41)]
        
        ca = CigarAlignment(cigar_to_cigartuple(cig1),cigar_to_cigartuple(cig2))
        ca.fill_matrix()
        aligned_cigars = ca.traceback_matrix()
        
        #print aligned_cigars[0]
        #print aligned_cigars[1]
        
        self.assertEqual(aligned_cigars[0],ans1)
        self.assertEqual(aligned_cigars[1],ans2)
        
        
        
        
        
    


def main():
    unittest.main()

if __name__ == '__main__':
    main()