#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import HTSeq

from drdisco import log


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


class Blacklist:
    """
    They came... with souls... from the blacklist :)
    """

    def __init__(self):
        self.idx_regions = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.idx_junctions = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.r = 0
        self.j = 0

    def add_regions_from_bed(self, regions_bed_file):
        log.info("Parsing regions blacklist file: " + str(regions_bed_file))

        header = True
        with open(regions_bed_file, 'r') as fh:
            for line in fh:
                if not header:
                    params = line.strip("\t\n ").split("\t")
                    if len(params) > 1:
                        for i in [1, 2]:
                            params[i] = int(params[i])

                        d = params[2] - params[1]

                        if d < 1:
                            raise ValueError("Too small region (starts are 0-based, ends are 1-based, like BED):\n" + line)

                        if len(params) >= 5:
                            self.add_region(params[0], params[1], params[2], params[3], params[4])
                        else:
                            self.add_region(params[0], params[1], params[2], params[3], None)
                else:
                    header = False

        log.info("Added " + str(self.r) + " regions to the blacklist")

    def add_region(self, _chr, _pos_s, _pos_e, _strand, id_suffix):
        uid = 'r-' + str(self.r)
        if id_suffix:
            uid += '-' + str(id_suffix)

        try:
            self.idx_regions[HTSeq.GenomicInterval(_chr, _pos_s, _pos_e, _strand)] += uid
        except:
            raise ValueError("Invalid entry in region file:", uid, ":", _chr, _pos_s, _pos_e, _strand, id_suffix)

        self.r += 1

    def add_junctions_from_file(self, junction_file):
        log.info("Parsing junction blacklist file: " + str(junction_file))
        header = True
        with open(junction_file, 'r') as fh:
            for line in fh:
                if not header:
                    params = line.strip().split("\t")
                    if len(params) > 1:
                        for i in [1, 2, 5, 6]:
                            params[i] = int(params[i])

                        d1 = params[2] - params[1]
                        d2 = params[6] - params[5]

                        if d1 < 1 or d2 < 1:
                            raise ValueError("Too small junction (starts are 0-based, ends are 1-based, like BED):\n" + line)

                        if (params[4] < params[0]) or (params[0] == params[4] and params[5] < params[1]):
                            reg1 = (params[4], params[5], params[6], params[7])
                            reg2 = (params[0], params[1], params[2], params[3])
                        else:
                            reg1 = (params[0], params[1], params[2], params[3])
                            reg2 = (params[4], params[5], params[6], params[7])

                        if len(params) >= 9:
                            self.add_junction(reg1, reg2, params[8])
                        else:
                            self.add_junction(reg1, reg2, None)
                else:
                    header = False

        log.info("Added " + str(self.j) + " junctions to the blacklist")

    def add_junction(self, reg1, reg2, id_suffix):
        uid = 'j-' + str(self.j)
        if id_suffix:
            uid += '-' + str(id_suffix)

        try:
            self.idx_junctions[HTSeq.GenomicInterval(reg1[0], reg1[1], reg1[2], reg1[3])] += uid
            self.idx_junctions[HTSeq.GenomicInterval(reg2[0], reg2[1], reg2[2], reg2[3])] += uid
        except:
            raise ValueError("Invalid entry in junction file:", uid, ":", reg1, reg2, id_suffix)

        self.j += 1

    def is_blacklisted_by_junctions(self, pos1, pos2):
        if (pos2[0] < pos1[0]) or (pos1[0] == pos2[0] and pos2[1] < pos1[1]):
            pos3 = pos1
            pos1 = pos2
            pos2 = pos3

        ids1 = set()
        position = self.idx_junctions[HTSeq.GenomicPosition(pos1[0], pos1[1], pos1[2])]
        for step in position:
            ids1.add(step)

        ids2 = set()
        position = self.idx_junctions[HTSeq.GenomicPosition(pos2[0], pos2[1], pos2[2])]
        for step in position:
            ids2.add(step)

        return ids1.intersection(ids2)

    def is_blacklisted_by_regions(self, pos1, pos2):
        if (pos2[0] < pos1[0]) or (pos1[0] == pos2[0] and pos2[1] < pos1[1]):
            pos3 = pos1
            pos1 = pos2
            pos2 = pos3

        ids = set()
        position = self.idx_regions[HTSeq.GenomicPosition(pos1[0], pos1[1], pos1[2])]
        for step in position:
            ids.add(step)

        position = self.idx_regions[HTSeq.GenomicPosition(pos2[0], pos2[1], pos2[2])]
        for step in position:
            ids.add(step)

        return ids

    def is_blacklisted(self, pos1, pos2):
        return self.is_blacklisted_by_junctions(pos1, pos2).union(self.is_blacklisted_by_regions(pos1, pos2))
