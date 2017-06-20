#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

import math

from __init__ import MIN_DISCO_PER_SUBNET_PER_NODE

from drdisco import log
import HTSeq


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


class DetectOutputEntry:
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
        self.circ_lin = self.line[8]
        self.x_onic = self.line[9]
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
        self.lr_A_slope = float(self.line[23])
        self.lr_A_intercept = self.line[24]
        self.lr_A_rvalue = float(self.line[25])
        self.lr_A_pvalue = self.line[26]
        self.lr_A_stderr = self.line[27]
        self.lr_B_slope = float(self.line[28])
        self.lr_B_intercept = self.line[29]
        self.lr_B_rvalue = float(self.line[30])
        self.lr_B_pvalue = self.line[31]
        self.lr_B_stderr = self.line[32]
        self.disco_split = self.line[33]
        self.clips_score = self.line[34]
        self.nodes_edge = float(self.line[35])
        self.structure = self.line[36]

    def get_donors_acceptors(self, gene_tree):
        def structure_to_acceptor_donor_order():
            idx = {}
            for a in self.structure.split('&'):
                for b in a.split(':',3)[3].strip('()').split(','):
                    c = b.split(':')
                    c[0] = c[0].replace('_1','_[12]').replace('_2','_[12]')
                    if c[0] != 'discordant_mates':
                        if c[0] not in idx:
                            idx[c[0]] = 0
                        
                        idx[c[0]] += int(c[1])
            
            print idx

            """
                sp[12] ratio: 184 / (184 + 5)
                sp[12]_t ratio: 5 / (184 + 5)

                resulting fusion:

                <<2 + <<1
            """
            
            return (), ()

        acceptor_pos, donor_pos = structure_to_acceptor_donor_order()

        #if ratio_highest == spanning_paired_12:
        #   acceptor = (self.chrB, self.posB, self.strandB)
        #   donor = (self.chrA, self.posA, self.strandA)

        return ['TMPRSS2'], ['ERG']


    def __str__(self):
        line = self.line
        line[7] = self.status
        return "\t".join(line) + "\n"


class DetectOutput:
    def __init__(self, input_results_file):
        self.input_alignment_file = input_results_file
        self.header = self.get_header()

    def get_header(self):
        with open(self.input_alignment_file, 'r') as fh_in:
            for line in fh_in:
                return line

    def __iter__(self):
        header = True
        with open(self.input_alignment_file, 'r') as fh_in:
            for line in fh_in:
                if not header:
                    e = DetectOutputEntry(line)
                    yield e
                else:
                    header = False

    def classify(self, output_file, only_valid, blacklist):
        log.info("Loading " + output_file + "[only_valid=" + {True: 'true', False: 'false'}[only_valid] + "]")
        n = 0
        k = 0

        with open(output_file, 'w') as fh:
            fh.write(self.get_header())

            for e in self:
                if isinstance(e, basestring):
                    fh.write(e)
                else:
                    status = []
                    n += 1

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

                    # @todo subfunc
                    log_ratio_slope_max = (3.6 / 2)
                    log_ratio_rvalue_max = (0.8 / 2)
                    log_ratio_slope = abs(math.log((e.lr_A_slope + 0.0001) / (e.lr_B_slope + 0.0001)))
                    log_ratio_rvalue = abs(math.log((e.lr_A_rvalue + 0.0001) / (e.lr_B_rvalue + 0.0001)))
                    if log_ratio_slope > log_ratio_slope_max:
                        status.append("log_ratio_slope=" + str(round(log_ratio_slope, 2)) + ">" + str(round(log_ratio_slope_max, 2)))
                    if log_ratio_rvalue > log_ratio_rvalue_max:
                        status.append("log_ratio_rvalue=" + str(round(log_ratio_rvalue, 2)) + ">" + str(round(log_ratio_rvalue_max, 2)))

                    if len(status) == 0:
                        e.status = 'valid'
                        fh.write(str(e))
                        k += 1
                    elif not only_valid:
                        e.status = ','.join(status)
                        fh.write(str(e))

        log.info("Classified " + str(k) + "/" + str(n) + " as valid")

    def integrate(self, output_table, gtf_file):
        def insert_in_index(index, entries, score):
            if score not in index:
                index[score] = {}

            key = entries[0].chrA + ':' + str(entries[0].posA) + '(' + entries[0].strandA + ')-' + entries[0].chrB + ':' + str(entries[0].posB) + '(' + entries[0].strandB + ')'
            index[score][key] = entries

        with open(output_table, 'w') as fh_out:
            fh_out.write("shared-id\tfusion\t" + self.header)
            self.idx = HTSeq.GenomicArrayOfSets("auto", stranded=True)

            intronic_linear = []
            remainder = []

            # Find 'duplicates' or fusions that belong to each other
            for e in self:
                if e.x_onic == 'intronic' and e.circ_lin == 'linear':
                    intronic_linear.append(e)
                else:
                    remainder.append(e)

                def insert(pos, e):
                    # position_accession = HTSeq.GenomicPosition(pos[0], pos[1], pos[2])
                    position_accession = HTSeq.GenomicInterval(pos[0], pos[1], pos[1] + 1, pos[2])
                    position = self.idx[position_accession]
                    position += e

                insert((e.chrA, e.posA, e.strandA), e)
                insert((e.chrB, e.posB, e.strandB), e)

            # Reorder
            idx2 = {}

            for e in intronic_linear:
                results = {}
                positions = [(e.chrA, e.posA, e.strandA), (e.chrB, e.posB, e.strandB)]

                for pos in positions:
                    if pos[2] == '-':
                        pos1 = pos[1] - 200000
                        pos2 = pos[1]
                    else:
                        pos1 = pos[1]
                        pos2 = pos[1] + 200000

                    for step in self.idx[HTSeq.GenomicInterval(pos[0], pos1, pos2, pos[2])].steps():
                        for e2 in step[1]:
                            if e != e2:
                                if e2 not in results:
                                    results[e2] = 0

                                results[e2] += 1

                top_result = (None, 9999999999999)
                for r in results:
                    if results[r] >= 2:
                        d1 = (r.posA - e.posA)
                        d2 = (r.posB - e.posB)
                        sq_d = math.sqrt(pow(d1, 2) + pow(d2, 2))

                        shared_score = max(e.score, r.score) - min(e.score, r.score)
                        penalty = 1.0 * sq_d / shared_score
                        if penalty < top_result[1]:
                            top_result = (r, penalty)

                if top_result[0]:
                    insert_in_index(idx2, [e, top_result[0]], e.score + top_result[0].score)
                else:
                    insert_in_index(idx2, [e], e.score)

            for e in remainder:
                insert_in_index(idx2, [e], e.score)

            # Generate output
            i = 1
            exported = set([])
            for score in sorted(idx2.keys(), reverse=True):
                for key in sorted(idx2[score].keys()):
                    added = 0
                    for entry in idx2[score][key]:
                        if entry not in exported:
                            donors, acceptors = entry.get_donors_acceptors(None)

                            fh_out.write(str(i) + "\t" + ','.join(donors)+'->'+','.join(acceptors) + "\t" + str(entry))
                            exported.add(entry)
                            added += 1

                    if added > 0:
                        i += 1
