#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

import math

from drdisco import log
from drdisco.DetectFrameShifts import DetectFrameShifts
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

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
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
            3. n-acceptor A
            4. n-donor A
            5. chr-B
            6. n-acceptor B
            7. n-donor B
            8. pos-B
            9. direction-B
            10. genomic-distance
            11. filter-status
            12. circRNA
            13. intronic/exonic
            14. score
        """
        self.chrA = self.line[0]
        self.posA = int(self.line[1])
        self.strandA = self.line[2]
        
        self.acceptorA = int(self.line[3])
        self.donorA = int(self.line[4])
        
        self.chrB = self.line[5]
        self.posB = int(self.line[6])
        self.strandB = self.line[7]
        
        self.acceptorB = int(self.line[8])
        self.donorB = int(self.line[9])
        
        self.dist = self.line[10]
        self.status = self.line[11]
        self.circ_lin = self.line[12]
        self.x_onic = self.line[13]
        self.score = int(self.line[14])
        self.clips = int(self.line[15])
        self.n_split_reads = int(self.line[16])
        self.n_discordant_reads = int(self.line[17])
        self.n_supporting_reads = self.n_split_reads + self.n_discordant_reads
        self.n_edges = self.line[18]
        self.n_nodes_A = int(self.line[19])
        self.n_nodes_B = int(self.line[20])
        self.n_nodes = self.n_nodes_A + self.n_nodes_B
        self.n_splice_junc_A = self.line[21]
        self.n_splice_junc_B = self.line[22]
        self.entropy_bp_edge = float(self.line[23])
        self.entropy_all_edges = float(self.line[24])
        self.bp_pos_stddev = float(self.line[25])
        self.entropy_disco_bps = self.line[26]
        self.lr_A_slope = float(self.line[27])
        self.lr_A_intercept = self.line[28]
        self.lr_A_rvalue = float(self.line[29])
        self.lr_A_pvalue = self.line[30]
        self.lr_A_stderr = self.line[31]
        self.lr_B_slope = float(self.line[32])
        self.lr_B_intercept = self.line[33]
        self.lr_B_rvalue = float(self.line[34])
        self.lr_B_pvalue = self.line[35]
        self.lr_B_stderr = self.line[36]
        self.disco_split = self.line[37]
        self.clips_score = self.line[38]
        self.nodes_edge = float(self.line[39])
        self.structure = self.line[40]
        
        inv = {'-': '+', '+': '-'}
        if self.acceptorA > self.donorA:
            # TMPRSS2 ERG DNA: - + 
            # TMPRSS2 ERG RNA: - -
            self.RNAstrandA = self.strandA
            self.RNAstrandB = inv[self.strandB]
        elif self.donorA < self.acceptorA:
            self.RNAstrandA = inv[self.strandA]
            self.RNAstrandB = self.strandB
        else:
            self.RNAstrandA = '.'
            self.RNAstrandB = '.'

        
        self.frameshift_0 = ''
        self.frameshift_1 = ''
        self.frameshift_2 = ''

    def get_donors_acceptors(self, gtf_file):
        idx = {}
        for a in self.structure.split('&'):
            for b in a.split(':', 3)[3].strip('()').split(','):
                c = b.split(':')
                c[0] = c[0].replace('_1', '_[12]').replace('_2', '_[12]')
                if c[0] != 'discordant_mates':
                    if c[0] not in idx:
                        idx[c[0]] = 0

                    idx[c[0]] += int(c[1])

        def pos_to_gene_str(pos_chr, pos_pos):
            pos = HTSeq.GenomicInterval(pos_chr, pos_pos, pos_pos + 1, ".")
            genes = set([])

            for step in gtf_file[pos]:
                genes = genes.union(step)

            genes_str = ','.join(sorted(list(genes)))
            if not genes:
                return pos_chr + ':' + str(pos_pos)
            else:
                return genes_str

        genesA = pos_to_gene_str(self.chrA, self.posA)
        genesB = pos_to_gene_str(self.chrB, self.posB)

        if self.donorA < self.acceptorA:
            return genesA + '->' + genesB
        elif self.donorA > self.acceptorA:
            return genesB + '->' + genesA
        else:
            return genesB + '<->' + genesA

    def __str__(self):
        line = self.line
        line[11] = self.status
        return "\t".join(line) + "\n"


class DetectOutput:
    def __init__(self, input_results_file):
        self.input_alignment_file = input_results_file
        self.header = self.get_header()

    def get_header(self):
        with open(self.input_alignment_file, 'r') as fh_in:
            for line in fh_in:
                return line
        raise Exception("Invalid file: " + str(self.input_alignment_file))

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

                    all_entropy_min = 0.705 + (math.atan((e.score - 150) * 0.005) * 0.035)
                    all_entropy_max = -1.0 * (max(e.score, 171) - 175.0) / (5.0 + max(e.score, 171) - 175.0) + (1.0 + 0.965)
                    if e.entropy_all_edges < all_entropy_min:
                        status.append("entropy=" + str(e.entropy_bp_edge) + '<' + str(round(all_entropy_min, 4)))
                    if e.entropy_all_edges > all_entropy_max:
                        status.append("entropy=" + str(e.entropy_bp_edge) + '>' + str(round(all_entropy_max, 4)))

                    # @todo subfunc
                    n_disco_min = int(round(pow(((e.n_nodes - 2) * 0.22), 1.7)))
                    if e.n_discordant_reads < n_disco_min:
                        status.append("n_discordant_reads=" + str(e.n_discordant_reads) + "<" + str(n_disco_min))

                    # @todo subfunc
                    n_support_min = (0.215 * pow(max(0, e.n_nodes) - 1.0, 1.59)) + 6.5
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
            gene_annotation = HTSeq.GenomicArrayOfSets("auto", stranded=False)
            dfs = None
            
            if gtf_file:
                dfs = DetectFrameShifts(gtf_file)
                gtf_file = HTSeq.GFF_Reader(gtf_file, end_included=True)
                for feature in gtf_file:
                    if feature.type == "gene":
                        if 'gene_name' in feature.attr:
                            name = feature.attr['gene_name']
                        elif 'Name' in feature.attr:
                            name = feature.attr['Name']
                        elif 'gene' in feature.attr:
                            name = feature.attr['gene']
                        else:
                            name = feature.name
                        gene_annotation[feature.iv] += name

            intronic_linear = []
            remainder = []

            # Find 'duplicates' or fusions that belong to each other
            for e in self:
                #print e
                #print (e.chrA,e.posA,e.RNAstrandA),(e.chrB,e.posB,e.RNAstrandB)
                
                #if dfs and e.RNAstrandA != '.' and e.RNAstrandB:
                    #print 'dna', (e.chrA,e.posA,e.strandA),(e.chrB,e.posB,e.strandB)
                    #print 'rna', (e.chrA,e.posA,e.RNAstrandA),(e.chrB,e.posB,e.RNAstrandB)
                    #print dfs.evaluate((e.chrA,e.posA,e.RNAstrandA), (e.chrB,e.posB,e.RNAstrandB))
                    #print
                
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

                    for step in self.idx[HTSeq.GenomicInterval(pos[0], max(0, pos1), pos2, pos[2])].steps():
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

                        shared_score = math.sqrt((pow(e.score, 2) + pow(r.score, 2)) * 0.5)
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
                            acceptors_donors = entry.get_donors_acceptors(gene_annotation)

                            fh_out.write(str(i) + "\t" + acceptors_donors + "\t" + str(entry))
                            exported.add(entry)
                            added += 1

                    if added > 0:
                        i += 1
