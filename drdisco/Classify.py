#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

import math

from __init__ import MIN_SUBNET_ENTROPY, MIN_DISCO_PER_SUBNET_PER_NODE, MIN_SUPPORTING_READS_PER_SUBNET_PER_NODE


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
        self.score = self.line[10]
        self.clips = self.line[11]
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
        self.bp_pos_stddev = self.line[21]
        self.entropy_disco_bps = self.line[22]

    def __str__(self):
        line = self.line
        line[7] = self.status
        return "\t".join(line) + "\n"


class Classify:
    def __init__(self, input_results_file):
        self.input_alignment_file = input_results_file

    def classify(self, output_file, only_valid):
        with open(output_file, 'w') as fh:
            header = True
            with open(self.input_alignment_file, 'r') as fh_in:
                for line in fh_in:
                    if not header:
                        status = []
                        e = Entry(line)

                        # @todo subfunc
                        if e.entropy_bp_edge < MIN_SUBNET_ENTROPY:
                            status.append("entropy=" + str(e.entropy_bp_edge) + '<' + str(MIN_SUBNET_ENTROPY))

                        # @todo subfunc
                        n_disco_min = MIN_DISCO_PER_SUBNET_PER_NODE * int(round(math.sqrt(e.n_nodes)))
                        if e.n_discordant_reads < n_disco_min:
                            status.append("n_discordant_reads=" + str(e.n_discordant_reads) + "<" + str(n_disco_min))

                        # @todo subfunc
                        n_support_min = (MIN_SUPPORTING_READS_PER_SUBNET_PER_NODE * e.n_nodes)
                        n_support_min_new = int(round(pow(1.2 * n_support_min, 0.913)))
                        if e.n_supporting_reads < n_support_min_new:
                            status.append("n_support=" + str(e.n_supporting_reads) + "<" + str(n_support_min))  # err msg is wrong, should be 'new'!

                        # @todo subfunc
                        n_disco_max = int(round(35 + (0.55 * e.n_split_reads)))
                        if e.n_discordant_reads > n_disco_max:
                            status.append("n_disco" + str(e.n_discordant_reads) + ">" + str(n_disco_max))

                        # @todo subfunc
                        n_split_min = int(round((0.52 * e.n_supporting_reads) - pow((0.1 * e.n_supporting_reads), 1.2) - 2))
                        if e.n_split_reads < n_split_min:
                            status.append("n_split" + str(e.n_split_reads) + "<" + str(n_split_min))

                        if len(status) == 0:
                            e.status = 'valid'
                            fh.write(str(e))
                        elif not only_valid:
                            e.status = ','.join(status)
                            fh.write(str(e))

                    else:
                        fh.write(line)
                        header = False
