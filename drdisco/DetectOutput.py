#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

import math

from drdisco import log
from drdisco.DetectFrameShifts import DetectFrameShifts
from drdisco.utils import reverse_complement, is_gzip
import gzip
import HTSeq
from tqdm import tqdm
from pyfaidx import Fasta


# python2 devision compatibility
import decimal
decimal.getcontext().rounding = decimal.ROUND_HALF_UP


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
        self.line = line_in_results_file.strip("\r\n").split("\t")
        self.parse()

    def __lt__(self, other_doe):
        return str(self) < str(other_doe)

    def __gt__(self, other_doe):
        return str(self) > str(other_doe)

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

        self.alignment_score = int(self.line[18])
        self.mismatches = int(self.line[19])

        self.n_edges = int(self.line[20])
        self.n_nodes_A = int(self.line[21])
        self.n_nodes_B = int(self.line[22])
        self.n_nodes = self.n_nodes_A + self.n_nodes_B
        self.n_splice_junc_A = self.line[23]
        self.n_splice_junc_B = self.line[24]
        self.entropy_bp_edge = float(self.line[25])
        self.entropy_all_edges = float(self.line[26])
        self.bp_pos_stddev = float(self.line[27])
        self.entropy_disco_bps = self.line[28]

        self.stddev_disco_bps_B = float(self.line[29])
        self.entropy_disco_bps_A = float(self.line[30])
        self.stddev_disco_bps_B = float(self.line[31])
        self.entropy_disco_bps_B = float(self.line[32])

        self.lr_A_slope = float(self.line[29 + 4])
        self.lr_A_intercept = float(self.line[30 + 4])
        self.lr_A_rvalue = float(self.line[31 + 4])
        self.lr_A_pvalue = float(self.line[32 + 4])
        self.lr_A_stderr = self.line[33 + 4]

        self.lr_B_slope = float(self.line[34 + 4])
        self.lr_B_intercept = float(self.line[35 + 4])
        self.lr_B_rvalue = float(self.line[36 + 4])
        self.lr_B_pvalue = float(self.line[37 + 4])
        self.lr_B_stderr = self.line[38 + 4]

        self.disco_split = self.line[39 + 4]  # ratio1
        self.clips_score = self.line[40 + 4]  # ratio2
        self.nodes_edge = float(self.line[41 + 4])  # ratio3 -> 1.0 * (nodes_a + nodes_b) / len(self.edges)

        self.break_A_median_AS = int(self.line[42 + 4])
        self.break_B_median_AS = int(self.line[43 + 4])
        self.break_A_max_AS = int(self.line[44 + 4])
        self.break_B_max_AS = int(self.line[45 + 4])

        self.edit_dist_to_splice_motif = ""

        self.exons_from = []
        self.exons_to = []

        self.structure = self.line[46 + 4]

        inv = {'-': '+', '+': '-'}
        if self.acceptorA > self.donorA:
            self.RNAstrandA = self.strandA
            self.RNAstrandB = inv[self.strandB]
        elif self.acceptorA < self.donorA:
            self.RNAstrandA = inv[self.strandA]
            self.RNAstrandB = self.strandB
        else:
            self.RNAstrandA = '.'
            self.RNAstrandB = '.'

        self.fgd = ''
        self.frameshift_0 = ''
        self.frameshift_1 = ''
        self.frameshift_2 = ''

    def get_donors_acceptors(self, gtf_file):
        idx = {}
        for a in self.structure.split('&'):
            if a != '':
                for b in a.split(':', 3)[3].strip('()').split(','):
                    c = b.split(':')
                    c[0] = c[0].replace('_1', '_[12]').replace('_2', '_[12]')
                    if c[0] != 'discordant_mates':
                        if c[0] not in idx:
                            idx[c[0]] = 0

                        idx[c[0]] += int(c[1])

        def pos_to_gene_str(pos_chr, pos_pos):
            if pos_chr[0:3] == 'chr':
                chrom = pos_chr[3:]
            else:
                chrom = pos_chr

            pos = HTSeq.GenomicInterval(chrom, pos_pos, pos_pos + 1, ".")
            genes = set([])

            for step in gtf_file.idx[pos]:
                genes = genes.union(step)

            genes_str = ','.join(sorted(list(genes)))
            if not genes:
                return pos_chr + ':' + str(pos_pos)
            else:
                return genes_str

        genesA = pos_to_gene_str(self.chrA, self.posA)
        genesB = pos_to_gene_str(self.chrB, self.posB)

        if self.donorA > self.donorB:
            return genesA + '->' + genesB
        elif self.donorB > self.donorA:
            return genesB + '->' + genesA
        else:
            return genesB + '<->' + genesA

    # def is_on_splice_junction_motif_long(self, fasta_fh):
        # """

# motif:

# 5' exon:

# [ ...{AC}{A}{G} ] {G}{T}{AG}{A}{G}{T} . . . {C}{A}{G} [ {G}... ]

        # """

        # pos5_in_exon_length = 3
        # pos5_post_exon_length = 6

        # pos3_pre_exon_length = 3
        # pos3_in_exon_length = 1

        # if self.donorA > self.donorB:
            # pos5p = [self.chrA, self.posA, self.strandA]
            # pos3p = [self.chrB, self.posB, self.strandB]
        # elif self.donorA < self.donorB:
            # pos5p = [self.chrB, self.posB, self.strandB]
            # pos3p = [self.chrA, self.posA, self.strandA]
        # else:
            # pos5p = None

        # if pos5p:
            # if pos5p[2] == '-':
                # try:
                    # seq_in_5p_exon = str(fasta_fh[pos5p[0]][pos5p[1] - pos5_in_exon_length:pos5p[1]]).upper()
                    # seq_post_5p_exon = str(fasta_fh[pos5p[0]][pos5p[1]:pos5p[1] + pos5_post_exon_length]).upper()
                # except ValueError:
                    # log.error("Accessing region that mismatches FASTA file: " + pos5p[0] + ":" + str(pos5p[1] - pos5_in_exon_length) + "-" + str(pos5p[1] + pos5_post_exon_length) + ". Skipping this for splice junction site estimation.")
                    # return ""
            # else:
                # try:
                    # seq_in_5p_exon = reverse_complement(str(fasta_fh[pos5p[0]][pos5p[1]:pos5p[1] + pos5_in_exon_length]))
                    # seq_post_5p_exon = reverse_complement(str(fasta_fh[pos5p[0]][pos5p[1] - pos5_post_exon_length:pos5p[1]]))
                # except ValueError:
                    # log.error("Accessing region that mismatches FASTA file: " + pos5p[0] + ":" + str(pos5p[1] - pos5_post_exon_length) + "-" + str(pos5p[1] + pos5_in_exon_length) + ". Skipping this for splice junction site estimation.")
                    # return ""

            # if pos3p[2] == '+':
                # try:
                    # seq_pre_3p_exon = str(fasta_fh[pos3p[0]][pos3p[1] - pos3_pre_exon_length:pos3p[1]]).upper()
                    # seq_in_3p_exon = str(fasta_fh[pos3p[0]][pos3p[1]:pos3p[1] + pos3_in_exon_length]).upper()
                # except ValueError:
                    # log.error("Accessing region that mismatches FASTA file: " + pos3p[0] + ":" + str(pos3p[1] - pos3_pre_exon_length) + "-" + str(pos3p[1] + pos3_in_exon_length) + ". Skipping this for splice junction site estimation.")
                    # return ""
            # else:
                # try:
                    # seq_in_3p_exon = reverse_complement(str(fasta_fh[pos3p[0]][pos3p[1] - pos3_in_exon_length:pos3p[1]]))
                    # seq_pre_3p_exon = reverse_complement(str(fasta_fh[pos3p[0]][pos3p[1]:pos3p[1] + pos3_pre_exon_length]))
                # except ValueError:
                    # log.error("Accessing region that mismatches FASTA file: " + pos3p[0] + ":" + str(pos3p[1] - pos3_in_exon_length) + "-" + str(pos3p[1] + pos3_pre_exon_length) + ". Skipping this for splice junction site estimation.")
                    # return ""

            # def calc_dist(pat, subseq):
                # d = 0

                # if len(pat) != len(subseq):  # may happen at beginning and end of chromosomes (suchs as chrM)
                    # return max(len(pat), len(subseq))
                # else:
                    # for i in range(len(pat)):
                        # if subseq[i] not in pat[i]:
                            # d += 1

                # return d

            # dist = calc_dist(["AC", "A", "G"], seq_in_5p_exon) + calc_dist(["G", "T", "AG", "A", "G", "T"], seq_post_5p_exon) + calc_dist(["C", "A", "G"], seq_pre_3p_exon) + calc_dist(["G"], seq_in_3p_exon)
            # # print "[ ... " + seq_in_5p_exon + " ] " + seq_post_5p_exon + " ... ... " + seq_pre_3p_exon + " [ " + seq_in_3p_exon + " ... ] ---> " + str(dist)
            # self.edit_dist_to_splice_motif = str(dist)

            # return dist
        # else:
            # return ""

    def is_on_splice_junction_motif(self, fasta_fh):
        """
motif:

5' exon: -------------------------------->> 3'
[ . . . ] {G}{T} . . . . . {A}{G} [ . . . ]
  . . .    G  T  . . . . .  A  G    . . .   DNA

inversed, like TMPRSS2-ERG:

3' <<----------------------------------- 5'
[ . . . ] {C}{T} . . . . . {A}{C} [ . . . ] RNA
  . . .    G  A  . . . . .  T  G    . . .   DNA
        """

        pos5_in_exon_length = 0
        pos5_post_exon_length = 2

        pos3_pre_exon_length = 2
        pos3_in_exon_length = 0

        if self.donorA > self.donorB:
            pos5p = [self.chrA, self.posA, self.strandA]
            pos3p = [self.chrB, self.posB, self.strandB]
        elif self.donorA < self.donorB:
            pos5p = [self.chrB, self.posB, self.strandB]
            pos3p = [self.chrA, self.posA, self.strandA]
        else:
            pos5p = None

        if pos5p:
            if pos5p[2] == '-':
                try:
                    # seq_in_5p_exon = ''  str(fasta_fh[pos5p[0]][pos5p[1] - pos5_in_exon_length:pos5p[1]]).upper()
                    seq_post_5p_exon = str(fasta_fh[pos5p[0]][pos5p[1]:pos5p[1] + pos5_post_exon_length]).upper()

                    if len(seq_post_5p_exon) != pos5_post_exon_length:
                        log.error("[1] Accessing region that mismatches FASTA file: " + pos5p[0] + ":" + str(pos5p[1] - pos5_in_exon_length) + "-" + str(pos5p[1] + pos5_post_exon_length) + ". Skipping this for splice junction site estimation.")
                        return ""

                except (ValueError, KeyError):
                    log.error("[1] Accessing region that mismatches FASTA file: " + pos5p[0] + ":" + str(pos5p[1] - pos5_in_exon_length) + "-" + str(pos5p[1] + pos5_post_exon_length) + ". Skipping this for splice junction site estimation.")
                    return ""
            else:
                try:
                    # seq_in_5p_exon = '' reverse_complement(str(fasta_fh[pos5p[0]][pos5p[1]:pos5p[1] + pos5_in_exon_length]))
                    seq_post_5p_exon = reverse_complement(str(fasta_fh[pos5p[0]][pos5p[1] - pos5_post_exon_length:pos5p[1]]))

                    # Later versions of pyfaidx return empty array without value error if out of bound
                    if len(seq_post_5p_exon) != pos5_post_exon_length:
                        log.error("[2] Accessing region that mismatches FASTA file: " + pos5p[0] + ":" + str(pos5p[1] - pos5_post_exon_length) + "-" + str(pos5p[1]) + ". Skipping this for splice junction site estimation.")
                        return ""

                except (ValueError, KeyError):
                    log.error("[2] Accessing region that mismatches FASTA file: " + pos5p[0] + ":" + str(pos5p[1] - pos5_post_exon_length) + "-" + str(pos5p[1]) + ". Skipping this for splice junction site estimation.")
                    return ""

            if pos3p[2] == '+':
                try:
                    seq_pre_3p_exon = str(fasta_fh[pos3p[0]][pos3p[1] - pos3_pre_exon_length:pos3p[1]]).upper()
                    # seq_in_3p_exon = ''  str(fasta_fh[pos3p[0]][pos3p[1]:pos3p[1] + pos3_in_exon_length]).upper()

                    if len(seq_pre_3p_exon) != pos3_pre_exon_length:
                        log.error("[3] Accessing region that mismatches FASTA file: " + pos3p[0] + ":" + str(pos3p[1] - pos3_pre_exon_length) + "-" + str(pos3p[1] + pos3_in_exon_length) + ". Skipping this for splice junction site estimation.")
                        return ""

                except (ValueError, KeyError):
                    log.error("[3] Accessing region that mismatches FASTA file: " + pos3p[0] + ":" + str(pos3p[1] - pos3_pre_exon_length) + "-" + str(pos3p[1] + pos3_in_exon_length) + ". Skipping this for splice junction site estimation.")
                    return ""
            else:
                try:
                    # seq_in_3p_exon = ''  reverse_complement(str(fasta_fh[pos3p[0]][pos3p[1] - pos3_in_exon_length:pos3p[1]]))
                    seq_pre_3p_exon = reverse_complement(str(fasta_fh[pos3p[0]][pos3p[1]:pos3p[1] + pos3_pre_exon_length]))

                    if len(seq_pre_3p_exon) != pos3_pre_exon_length:
                        log.error("[4] Accessing region that mismatches FASTA file: " + pos3p[0] + ":" + str(pos3p[1] - pos3_in_exon_length) + "-" + str(pos3p[1] + pos3_pre_exon_length) + ". Skipping this for splice junction site estimation.")
                        return ""

                except (ValueError, KeyError):
                    log.error("[4] Accessing region that mismatches FASTA file: " + pos3p[0] + ":" + str(pos3p[1] - pos3_in_exon_length) + "-" + str(pos3p[1] + pos3_pre_exon_length) + ". Skipping this for splice junction site estimation.")
                    return ""

            def calc_dist(pat, subseq):
                d = 0

                if len(pat) != len(subseq):  # may happen at beginning and end of chromosomes (suchs as chrM)
                    return max(len(pat), len(subseq))
                else:
                    for i in range(len(pat)):
                        if subseq[i] not in pat[i]:
                            d += 1

                return d

            dist = 0 + calc_dist(["G", "T"], seq_post_5p_exon) + calc_dist(["A", "G"], seq_pre_3p_exon) + 0
            self.edit_dist_to_splice_motif = str(dist) + "|" + seq_post_5p_exon + "->" + seq_pre_3p_exon

            return self.edit_dist_to_splice_motif
        else:
            return ""

    def __str__(self):
        line = self.line
        line[11] = self.status
        line[13] = self.x_onic
        return "\t".join(line) + "\n"


class GeneAnnotation:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.idx = self.parse()

    def parse(self):
        idx = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        if self.gtf_file:  # could be None of no gtf file is provided
            log.info("Loading " + self.gtf_file)

            gtf_file = HTSeq.GFF_Reader(self.gtf_file, end_included=True)
            n = 0

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

                    if feature.iv.chrom[0:3] == 'chr':
                        feature.iv.chrom = feature.iv.chrom[3:]

                    idx[feature.iv] += name
                    n += 1

            log.info("Loaded " + str(n) + " features")
        return idx


class DetectOutput:
    def __init__(self, input_results_file):
        self.input_alignment_file = input_results_file
        self.is_gzip = is_gzip(input_results_file)
        self.header = self.get_header()

    def get_header(self):
        if self.is_gzip:
            with gzip.open(self.input_alignment_file, 'rb') as fh_in:
                for line in fh_in:
                    return line.decode("utf-8")
        else:
            with open(self.input_alignment_file, 'r') as fh_in:
                for line in fh_in:
                    return line
        raise Exception("Invalid file: " + str(self.input_alignment_file))

    def __iter__(self):
        header = True

        if self.is_gzip:
            with gzip.open(self.input_alignment_file, 'r') as fh_in:
                for line in fh_in:
                    if not header:
                        e = DetectOutputEntry(line.decode("utf-8"))
                        yield e
                    else:
                        header = False
        else:
            with open(self.input_alignment_file, 'r') as fh_in:
                for line in fh_in:
                    if not header:
                        e = DetectOutputEntry(line)
                        yield e
                    else:
                        header = False

    def classify(self, output_file, only_valid, blacklist, min_chim_overhang, ffpe_mismatch_ratio):
        log.info("Loading " + output_file + "[only_valid=" + {True: 'true', False: 'false'}[only_valid] + "]")
        n = 0
        k = 0

        with open(output_file, 'w') as fh:
            fh.write(str(self.get_header()))

            for e in self:
                if isinstance(e, str):
                    fh.write(e)
                else:
                    def classify_intronic_exonic():
                        n_edges_max = int(round(0.00575 * e.score + 5.75, 0))

                        if e.n_edges >= n_edges_max:
                            e.x_onic = 'exonic'

                    classify_intronic_exonic()

                    status = []
                    n += 1

                    # all_entropy_min_f1 = 0.705 + (math.atan((e.score - 150.0) * 0.005) * 0.035)
                    # all_entropy_min_f2 = -0.1375 * math.tanh((e.score - 2250.0) / 900.0) + 0.6175
                    all_entropy_min = -0.26 * math.tanh((e.score - 25.0) / 20.0) + 0.6225
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
                    n_disco_max = int(round(math.pow(22.0 * e.n_split_reads, 0.9) + 13))
                    if e.n_split_reads < 100:
                        n_disco_min = int(round(math.pow(0.0195 * e.n_split_reads, 1.95)))
                    elif e.n_split_reads >= 100 and e.n_split_reads < 125:
                        n_disco_min = 4
                    elif e.n_split_reads >= 125 and e.n_split_reads < 325:
                        n_disco_min = int(round((0.135 * (e.n_split_reads - 200.0)) + 14.0))
                    else:
                        n_disco_min = int(round(30.875 + (e.n_split_reads - 325) * 0.024))
                    if e.n_discordant_reads > n_disco_max:
                        status.append("n_disco=" + str(e.n_discordant_reads) + ">" + str(n_disco_max))
                    if e.n_discordant_reads < n_disco_min:
                        status.append("n_disco=" + str(e.n_discordant_reads) + "<" + str(n_disco_min))

                    # @todo subfunc
                    n_split_min = int(round((0.32 * e.n_supporting_reads) - pow((0.1 * e.n_supporting_reads), 1.15) - 4.0))
                    if e.n_supporting_reads < 385:
                        n_split_max = int(round((0.986 * e.n_supporting_reads) - pow(0.00535 * e.n_supporting_reads, 3.99 - ((1.0 / 15000.0) * e.n_supporting_reads))))
                    else:
                        n_split_max = int(round(0.94 * e.n_supporting_reads))
                    if e.n_split_reads < n_split_min:
                        status.append("n_split=" + str(e.n_split_reads) + "<" + str(n_split_min))
                    if e.n_split_reads > n_split_max:
                        status.append("n_split=" + str(e.n_split_reads) + ">" + str(n_split_max))

                    # @todo subfunc
                    slope = 51.0
                    bp_pos_stddev_max = -(slope * e.nodes_edge) + 15 + (2.0 * slope)
                    if e.bp_pos_stddev > bp_pos_stddev_max:
                        status.append("bp_pos_stddev=" + str(e.bp_pos_stddev) + ">" + str(bp_pos_stddev_max))

                    # @todo subfunc
                    clips_min = (0.19 * e.score) - 25.0
                    clips_max = (0.84 * e.score) + 550.0
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

                    # @todo subfunc
                    # FFPE material seems to have a substantial higher amount of mismatches per base, though randomly distributed
                    # if we ever make a v2 of dr-disco that incorporates the concordant reads, this variable can be dertemined with some kind of calibration
                    # now we're only estimating the MM ratio without entropy per position
                    # [CGCGCTATATCTCGATCGCCCTTAGAGATCCTTTCGAGAGAGCTCTAGAGCG] SOME KIND OF REFERENCE SEQUENCE
                    #  CGCG*TATAT*TC                  TTTC*AGAGAGCT*TAG      The more randomly dispersed mismatches are more trustworthy (right side example)
                    #  CGCG*TATAT*TCGAT                TTCGAGAG*GCTCT
                    #   GCG*TATAT*TCG                  T*CGAGAGAG*TCTA
                    #   GCG*TATAT*TCGA                 TTCG*GAGAGCTCTA
                    #    CG*TATAT*TCGAT                TTCGAG*GAGCTCTAG
                    #     G*TATAT*TCG                   TCGA*AGA*CTCT
                    #
                    if ffpe_mismatch_ratio:
                        log_value_max = -6.45 - ((e.score + 6750.0) / (4400.0 - (e.score + 6750.0)))
                    else:
                        log_value_max = -4.7
                    log_value = math.log((float(e.mismatches) + 0.0000001) / float(e.alignment_score))
                    if log_value > log_value_max:
                        status.append("many_muts=" + str(round(log_value, 2)) + ">" + str(round(log_value_max, 2)))

                    # @todo subfunc
                    lr_a = e.lr_A_pvalue * e.lr_A_intercept
                    lr_b = e.lr_A_pvalue * e.lr_B_intercept
                    lr_symmetry_max = -e.score / (0.11 + (0.0246 * (e.score))) + 41
                    n_lr_symmetry = pow(pow(lr_a, 2) + pow(lr_b, 2), 0.5)
                    if n_lr_symmetry >= lr_symmetry_max:
                        status.append("n_lr_symmetry=" + str(round(n_lr_symmetry, 2)) + ">=" + str(round(lr_symmetry_max, 2)))

                    # @todo subfunc
                    chim_overhang = min(e.break_A_max_AS, e.break_B_max_AS)
                    if chim_overhang < min_chim_overhang:
                        status.append("chim_overhang=" + str(chim_overhang) + "<" + str(min_chim_overhang))

                    # @todo subfunc
                    if e.score <= 150:
                        lr_intercept_max = (-31.0 * ((e.score + 100.0) / (1800.0 + e.score + 100.0))) + 85.5
                    else:
                        lr_intercept_max = ((e.score - 150.0) * 0.0225) + 81.71951
                    if e.lr_A_intercept > lr_intercept_max:
                        status.append("lr_A_intercept=" + str(e.lr_A_intercept) + ">" + "{:.12g}".format(lr_intercept_max))
                    if e.lr_B_intercept > lr_intercept_max:
                        status.append("lr_B_intercept=" + str(e.lr_B_intercept) + ">" + "{:.12g}".format(lr_intercept_max))

                    # @todo subfunc
                    sqrt_entropy_bps_ab = pow(pow(e.entropy_disco_bps_A, 2) + pow(e.entropy_disco_bps_B, 2), 0.5)
                    if e.entropy_all_edges <= 0.85:
                        sqrt_entropy_bps_ab_max = 0.475
                    else:
                        sqrt_entropy_bps_ab_max = 3.4 * e.entropy_all_edges - 2.415
                    if sqrt_entropy_bps_ab > sqrt_entropy_bps_ab_max:
                        status.append("sqrt_entropy_bps_ab=" + "{:.12g}".format(sqrt_entropy_bps_ab) + ">" + str(round(sqrt_entropy_bps_ab_max,5)))

                    if len(status) == 0:
                        e.status = 'valid'
                        fh.write(str(e))
                        k += 1
                    elif not only_valid:
                        e.status = ','.join(status)
                        fh.write(str(e))

        log.info("Classified " + str(k) + "/" + str(n) + " as valid")

    def integrate(self, output_table, gtf_file, fasta_file):
        log.info("Integrating results")
        
        def insert_in_index(index, entries, score, i):
            if score not in index:
                index[score] = {}

            key = entries[0].chrA + ':' + str(entries[0].posA) + '(' + entries[0].strandA + ')-' + entries[0].chrB + ':' + str(entries[0].posB) + '(' + entries[0].strandB + ')|' + str(i)
            index[score][key] = entries

        with open(output_table, 'w') as fh_out:
            header = self.header.split("\t")
            header = "\t".join(header[:-5] + ['full-gene-dysregulation', 'frameshift=0', 'frameshift=+1', 'frameshift=+2', 'splice-motif-edit-distance', "exons from (5')", "exons to (3')"] + header[-5:])

            fh_out.write("shared-id\tfusion\t" + header)

            # index used to find duplicates
            self.idx = HTSeq.GenomicArrayOfSets("auto", stranded=True)

            # index used to annotate gene names: TMPRSS2->ERG
            gene_annotation = GeneAnnotation(gtf_file)
            dfs = DetectFrameShifts(gtf_file) if gtf_file else None

            ffs = Fasta(fasta_file) if fasta_file else None

            intronic_linear = []
            remainder = []

            # Find 'duplicates' or fusions that belong to each other
            log.info("Searching for intronic and exonic breaks that belong to the same event")
            for e in self:
                if dfs and e.RNAstrandA != '.' and e.RNAstrandB != '.':
                    done_breaks = set([])

                    if e.donorA > e.donorB:
                        exons_from, exons_to, frame_shifts = dfs.evaluate([e.chrA, e.posA, e.RNAstrandA], [e.chrB, e.posB, e.RNAstrandB], 2)
                    else:
                        exons_from, exons_to, frame_shifts = dfs.evaluate([e.chrB, e.posB, e.RNAstrandB], [e.chrA, e.posA, e.RNAstrandA], 2)

                    done_breaks.add(e.chrA + ':' + str(e.posA) + '/' + str(e.posA + 1) + '(' + e.strandA + ')->' + e.chrB + ':' + str(e.posB) + '/' + str(e.posB + 1) + '(' + e.strandB + ')')

                    fgd = [x[0] + '->' + x[1] for x in frame_shifts['fgd']]
                    frameshifts_0 = [x[0][0] + '->' + x[1][0] for x in frame_shifts[0]]
                    frameshifts_1 = [x[0][0] + '(+' + str(x[0][1]) + ')->' + x[1][0] + '(+' + str(x[1][1]) + ')' for x in frame_shifts[1]]
                    frameshifts_2 = [x[0][0] + '(+' + str(x[0][1]) + ')->' + x[1][0] + '(+' + str(x[1][1]) + ')' for x in frame_shifts[2]]

                    for additional_breaks in e.structure.split('&'):
                        if additional_breaks != '':
                            params = additional_breaks.split(':(')
                            n_split_reads = sum([int(x.split(':')[1]) for x in params[1].rstrip(')').split(',') if x.split(':')[0] != 'discordant_mates'])

                            posAB = params[0].split(':')
                            posA, posB = int(posAB[1].split('/')[0]), int(posAB[2].split('/')[0])

                            if params[0] not in done_breaks and n_split_reads > 0:
                                if e.donorA > e.donorB:  # nice, use same thing to swap if necessary
                                    exons_from_, exons_to_, frame_shifts = dfs.evaluate([e.chrA, posA, e.RNAstrandA], [e.chrB, posB, e.RNAstrandB], 2)
                                else:
                                    exons_from_, exons_to_, frame_shifts = dfs.evaluate([e.chrB, posB, e.RNAstrandB], [e.chrA, posA, e.RNAstrandA], 2)

                                exons_from += exons_from_
                                exons_to += exons_to_
                                del(exons_from_, exons_to_)

                                fgd += [x[0] + '->' + x[1] for x in frame_shifts['fgd']]
                                frameshifts_0 += [x[0][0] + '->' + x[1][0] for x in frame_shifts[0]]
                                frameshifts_1 += [x[0][0] + '(+' + str(x[0][1]) + ')->' + x[1][0] + '(+' + str(x[1][1]) + ')' for x in frame_shifts[1]]
                                frameshifts_2 += [x[0][0] + '(+' + str(x[0][1]) + ')->' + x[1][0] + '(+' + str(x[1][1]) + ')' for x in frame_shifts[2]]

                            done_breaks.add(params[0])

                    e.exons_from = sorted(list(set(exons_from)))
                    e.exons_to = sorted(list(set(exons_to)))
                    del(exons_from, exons_to)

                    e.fgd = ','.join(sorted(list(set(fgd))))
                    e.frameshift_0 = ','.join(sorted(list(set(frameshifts_0))))
                    e.frameshift_1 = ','.join(sorted(list(set(frameshifts_1))))
                    e.frameshift_2 = ','.join(sorted(list(set(frameshifts_2))))
                    del(fgd, frameshifts_0, frameshifts_1, frameshifts_2)

                if ffs:
                    e.is_on_splice_junction_motif(ffs)

                if e.x_onic == 'intronic' and e.circ_lin == 'linear':
                    intronic_linear.append(e)
                else:
                    remainder.append(e)

                def insert(pos, e):
                    if pos[0][0:3] == 'chr':
                        chrom = pos[0][3:]
                    else:
                        chrom = pos[0]

                    # position_accession = HTSeq.GenomicPosition(pos[0], pos[1], pos[2])
                    position_accession = HTSeq.GenomicInterval(chrom, pos[1], pos[1] + 1, pos[2])
                    position = self.idx[position_accession]
                    position += e

                insert((e.chrA, e.posA, e.strandA), e)
                insert((e.chrB, e.posB, e.strandB), e)

            if ffs != None:
                ffs.close()

            # Reorder
            log.info("Re-order and find matching entries")
            idx2 = {}
            q = 0
            for e in intronic_linear:
                results_split = [set([]), set([])]
                positions = [(e.chrA, e.posA, e.strandA), (e.chrB, e.posB, e.strandB)]

                for pos_i in [0,1]:
                    pos = positions[pos_i]

                    if pos[2] == '-':
                        pos1 = pos[1] - 200000
                        pos2 = pos[1]
                    else:
                        pos1 = pos[1]
                        pos2 = pos[1] + 200000

                    if pos[0][0:3] == 'chr':
                        chrom = pos[0][3:]
                    else:
                        chrom = pos[0]

                    for step in self.idx[HTSeq.GenomicInterval(chrom, max(0, pos1), pos2, pos[2])].steps():
                        for e2 in [_ for _ in step[1] if _ != e]:
                            if e2.strandA == e.strandA and e2.strandB == e.strandB:
                                results_split[pos_i].add(e2)

                results = results_split[0].intersection(results_split[1])
                top_result = (None, 9999999999999)
                for r in results:
                    d1 = (r.posA - e.posA)
                    d2 = (r.posB - e.posB)
                    sq_d = math.sqrt(pow(d1, 2) + pow(d2, 2))

                    shared_score = math.sqrt((pow(e.score, 2) + pow(r.score, 2)) * 0.5)
                    penalty = 1.0 * sq_d / shared_score

                    if penalty < top_result[1]:
                        top_result = (r, penalty)

                if top_result[0]:
                    insert_in_index(idx2, [e, top_result[0]], e.score + top_result[0].score, q)
                else:
                    insert_in_index(idx2, [e], e.score, q)
                q += 1

            for e in remainder:
                insert_in_index(idx2, [e], e.score, q)
                q += 1

            log.info("Determining fusion gene names and generate output")
            # Generate output
            i = 1
            exported = set([])
            for score in sorted(idx2.keys(), reverse=True):
                for key in sorted(idx2[score].keys()):
                    added = 0
                    for entry in idx2[score][key]:
                        if entry not in exported:
                            acceptors_donors = entry.get_donors_acceptors(gene_annotation)
                            line = entry.line[:-5] + [entry.fgd, entry.frameshift_0, entry.frameshift_1, entry.frameshift_2, entry.edit_dist_to_splice_motif, ",".join(entry.exons_from), ",".join(entry.exons_to)] + entry.line[-5:]

                            fh_out.write(str(i) + "\t" + acceptors_donors + "\t" + "\t".join(line) + "\n")
                            exported.add(entry)
                            added += 1

                    if added > 0:
                        i += 1
