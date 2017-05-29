#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

import math

from __init__ import MIN_DISCO_PER_SUBNET_PER_NODE

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

    You can e-mail me via 'yhoogstrate' at the following webmail domain:
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

        self.idx_regions[HTSeq.GenomicInterval(_chr, _pos_s, _pos_e, _strand)] += uid

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
                            raise ValueError("Too small region (starts are 0-based, ends are 1-based, like BED):\n" + line)

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

        self.idx_junctions[HTSeq.GenomicInterval(reg1[0], reg1[1], reg1[2], reg1[3])] += uid
        self.idx_junctions[HTSeq.GenomicInterval(reg2[0], reg2[1], reg2[2], reg2[3])] += uid

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


class Entry:
    def __init__(self, line_in_results_file):
        self.line = line_in_results_file.strip().split("\t")
        self.parse()

    def parse(self):
        """
            0. chr-A
            1. pos-A
            2. direction-A
            3. chr-B
            4. pos-B
            5. direction-B
            6. genomic-distance
            7. filter-status
            8. circRNA
            9. intronic/exonic
            10. score
            11. soft+hardclips
            12. n-split-reads
            13. n-discordant-reads
            14. n-edges
            15. n-nodes-A
            16. n-nodes-B
            17. n-splice-junc-A
            18. n-splice-junc-B
            19. entropy-bp-edge
            20. entropy-all-edges
            21. bp-pos-stddev
            22. entropy-disco-bps
            23. data-structure
        """
        self.chrA = self.line[0]
        self.posA = int(self.line[1])
        self.strandA = self.line[2]

        self.chrB = self.line[3]
        self.posB = int(self.line[4])
        self.strandB = self.line[5]

        self.dist = self.line[6]
        self.status = self.line[7]
        self.score = int(self.line[10])
        self.clips = int(self.line[11])
        self.n_split_reads = int(self.line[12])
        self.n_discordant_reads = int(self.line[13])
        self.n_supporting_reads = self.n_split_reads + self.n_discordant_reads
        self.n_edges = self.line[14]
        self.n_nodes_A = int(self.line[15])
        self.n_nodes_B = int(self.line[16])
        self.n_nodes = self.n_nodes_A + self.n_nodes_B
        self.n_splice_junc_A = self.line[17]
        self.n_splice_junc_B = self.line[18]
        self.entropy_bp_edge = float(self.line[19])
        self.entropy_all_edges = float(self.line[20])
        self.bp_pos_stddev = float(self.line[21])
        self.entropy_disco_bps = self.line[22]
        self.lr_A_slope = self.line[23]
        self.lr_A_intercept = self.line[24]
        self.lr_A_rvalue = self.line[25]
        self.lr_A_pvalue = self.line[26]
        self.lr_A_stderr = self.line[27]
        self.lr_B_slope = self.line[28]
        self.lr_B_intercept = self.line[29]
        self.lr_B_rvalue = self.line[30]
        self.lr_B_pvalue = self.line[31]
        self.lr_B_stderr = self.line[32]
        self.disco_split = self.line[33]
        self.clips_score = self.line[34]
        self.nodes_edge = float(self.line[35])

    def __str__(self):
        line = self.line
        line[7] = self.status
        return "\t".join(line) + "\n"


class Classify:
    def __init__(self, input_results_file):
        self.input_alignment_file = input_results_file

    def classify(self, output_file, only_valid, blacklist):
        log.info("Loading " + output_file + "[only_valid=" + {True: 'true', False: 'false'}[only_valid] + "]")
        n = 0
        k = 0
        with open(output_file, 'w') as fh:
            header = True
            with open(self.input_alignment_file, 'r') as fh_in:
                for line in fh_in:
                    if not header:
                        n += 1
                        status = []
                        e = Entry(line)

                        all_entropy_min = -1.0 * (e.score) / (1.2 + e.score) + (1.0 + 0.74)
                        all_entropy_max = -1.0 * (max(e.score, 171) - 175.0) / (5.0 + max(e.score, 171) - 175.0) + (1.0 + 0.965)
                        if e.entropy_all_edges < all_entropy_min:
                            status.append("entropy=" + str(e.entropy_bp_edge) + '<' + str(round(all_entropy_min, 4)))
                        if e.entropy_all_edges > all_entropy_max:
                            status.append("entropy=" + str(e.entropy_bp_edge) + '>' + str(round(all_entropy_max, 4)))

                        # @todo subfunc
                        n_disco_min = MIN_DISCO_PER_SUBNET_PER_NODE * int(round(math.sqrt(e.n_nodes)))
                        if e.n_discordant_reads < n_disco_min:
                            status.append("n_discordant_reads=" + str(e.n_discordant_reads) + "<" + str(n_disco_min))

                        # @todo subfunc
                        n_support_min = (0.12 * pow(max(0, e.n_nodes), 1.7)) + 6.5
                        n_support_min = int(round(n_support_min))

                        if e.n_supporting_reads < n_support_min:
                            status.append("n_support=" + str(e.n_supporting_reads) + "<" + str(n_support_min))

                        # @todo subfunc
                        # n_disco_max = int(round(35 + (0.55 * e.n_split_reads)))
                        n_disco_max = int(round(math.pow(22 * e.n_split_reads, 0.9) + 13))
                        n_disco_min = int(round(math.pow(0.0195 * e.n_split_reads, 1.95)))
                        if e.n_discordant_reads > n_disco_max:
                            status.append("n_disco=" + str(e.n_discordant_reads) + ">" + str(n_disco_max))
                        if e.n_discordant_reads < n_disco_min:
                            status.append("n_disco=" + str(e.n_discordant_reads) + "<" + str(n_disco_min))

                        # @todo subfunc
                        n_split_min = int(round((0.32 * e.n_supporting_reads) - pow((0.1 * e.n_supporting_reads), 1.15) - 4))
                        n_split_max = int(round((0.985 * e.n_supporting_reads) - pow(0.014 * e.n_supporting_reads, 2.20 - ((1 / 15000) * e.n_supporting_reads))))
                        if e.n_split_reads < n_split_min:
                            status.append("n_split=" + str(e.n_split_reads) + "<" + str(n_split_min))
                        if e.n_split_reads > n_split_max:
                            status.append("n_split=" + str(e.n_split_reads) + ">" + str(n_split_max))

                        # @todo subfunc
                        slope = 51
                        bp_pos_stddev_max = -(slope * e.nodes_edge) + 15 + (2 * slope)
                        if e.bp_pos_stddev > bp_pos_stddev_max:
                            status.append("bp_pos_stddev=" + str(e.bp_pos_stddev) + ">" + str(bp_pos_stddev_max))

                        # @todo subfunc
                        clips_min = (0.19 * e.score) - 25
                        clips_max = (0.78 * e.score) + 97
                        if e.clips < clips_min:
                            status.append("clips=" + str(e.clips) + "<" + str(clips_min))
                        if e.clips > clips_max:
                            status.append("clips=" + str(e.clips) + ">" + str(clips_max))

                        # @todo subfunc
                        blacklisted = blacklist.is_blacklisted((e.chrA, e.posA, e.strandA), (e.chrB, e.posB, e.strandB))
                        if len(blacklisted) > 0:
                            status.append("blacklist=" + '&'.join(blacklisted))

                        if len(status) == 0:
                            e.status = 'valid'
                            fh.write(str(e))
                            k += 1
                        elif not only_valid:
                            e.status = ','.join(status)
                            fh.write(str(e))

                    else:
                        fh.write(line)
                        header = False

        log.info("Classified " + str(k) + "/" + str(n) + " as valid")
