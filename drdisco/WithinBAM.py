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

#from __future__ import annotations
#from abc import ABC, abstractmethod


import os
import shutil
import hashlib
import random
import string

import pysam

from tqdm import tqdm


from drdisco import __version__
from drdisco import log
from drdisco.alignment import alignment
from drdisco.utils import str_to_bytearray


class WithinBAM(alignment):
    #def __init__(self, input_alignment_file):
    #    self.input_alignment_file = input_alignment_file

    def obtain_chimeric_read_names(self) -> set:
        names = set([])
        
        log.info('Searching for chimeric reads (read names) in bam-file')
        
        with pysam.AlignmentFile(self.input_alignment_file, "rb") as samfile_in:
            with tqdm(total = samfile_in.mapped) as pbar:

                for read in samfile_in.fetch():
                    if read.has_tag('SA') or read.is_supplementary or read.flag in [65,81,97,113,129,145,161,177]: # last one are not properly mapped + not suppl
                        names.add(read.qname)
                    pbar.update(1)

        log.info('Found '+str(len(names))+' chimeric reads (read names) in bam-file')

        return names

    def extract_chimeric_reads(self, output_alignment_file: str, temp_dir):
        
        names = self.obtain_chimeric_read_names()

        log.info('Extracting '+str(len(names))+' chimeric reads from bam file')

        with pysam.AlignmentFile(self.input_alignment_file, "rb") as samfile_in:
            with pysam.AlignmentFile(output_alignment_file, "wb", header=samfile_in.header) as samfile_out:
                with tqdm(total = samfile_in.mapped) as pbar:

                    for read in samfile_in.fetch():
                        if read.qname in names:
                            samfile_out.write(read)
                        pbar.update(1)



