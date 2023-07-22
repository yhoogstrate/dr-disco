#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]

    Dr. Disco: fusion gene detection in random hexamer RNA-seq data
    Copyright (C) 2017  Youri Hoogstrate

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    You can contact me via the github repository at the following url:
    <https://github.com/yhoogstrate/dr-disco>

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
    gmail dot com
"""


import pysam
from tqdm import tqdm


from drdisco import __version__
from drdisco import log


class alignment():
    def __init__(self, input_alignment_file):
        self.input_alignment_file = input_alignment_file

    def extract_reads_by_names(self, read_names: set, output_file: str) -> int:
        found = 0
        log.info("Accessing alignment file, trying to extract "+str(len(read_names)) + " unique reads")
        
        with pysam.AlignmentFile(self.input_alignment_file, "rb") as fh_in:
            fh_out = pysam.AlignmentFile(output_file, "wb", header = fh_in.header)

            for read in tqdm(fh_in, total=fh_in.mapped):
                if read.qname in read_names:
                    fh_out.write(read)
                    found += 1

        log.info("Extracted " + str(found) + "/" + str(len(read_names)) + " unique reads")

        return found


    def get_star_version(self) -> tuple:
        try:
            fh = pysam.AlignmentFile(self.input_alignment_file, "rb")
        except:
            fh = pysam.AlignmentFile(self.input_alignment_file, "r")

        header = fh.header()

        fh.close()
        
        
        print(header)
