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

from drdisco.ChimericAlignment import ChimericAlignment


class TestChimericAlignment(unittest.TestCase):
    def test_01(self):
        print("\n")
        
        if not os.path.exists("tmp"):
            os.mkdir("tmp")
        
        input_file =    "tests/fix-chimeric/test_terg_01.filtered.bam"
        output_file =   "tmp/test_terg_01.filtered.fixed.bam"
        output_file_s = "tmp/test_terg_01.filtered.fixed.sam"
        
        test_file = "tests/fix-chimeric/test_terg_01.filtered.fixed.sam"
        
        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file,"tmp")
        
        
        # Bam2Sam
        fhq = open(output_file_s,"w")
        fhq.write(pysam.view(output_file))
        fhq.close()
        
        
        self.assertTrue(filecmp.cmp(output_file_s, test_file))

    def test_02(self):
        print("\n")
        
        if not os.path.exists("tmp"):
            os.mkdir("tmp")
        
        input_file =    "tests/fix-chimeric/test_terg_02.filtered.bam"
        output_file =   "tmp/test_terg_02.filtered.fixed.bam"
        output_file_s = "tmp/test_terg_02.filtered.fixed.sam"
        
        test_file = "tests/fix-chimeric/test_terg_02.filtered.fixed.sam"
        
        alignment_handle = ChimericAlignment(input_file)
        alignment_handle.convert(output_file,"tmp")
        
        
        # Bam2Sam
        fhq = open(output_file_s,"w")
        fhq.write(pysam.view(output_file))
        fhq.close()
        
        
        self.assertTrue(filecmp.cmp(output_file_s, test_file))

def main():
    unittest.main()

if __name__ == '__main__':
    main()
