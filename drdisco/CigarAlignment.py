#!/usr/bin/env python
#  *- coding: utf-8 -*-
#  vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
Aligns CIGAR strings by STAR to aligned softclipped and matches over
each other again. This allows to calculate back which read was the
first chunk.

e.g.:

cig1: 100S 25M
cig2: 100M 25S

are already oredered and aligned under each other, just a weight measure
of the M/(M+S) ratio in the first half will determine that cig2
corresponds to the first chunk of the read.

If the matching chunk however also is softclipping (for just a few
bases) we need an alignment to solve it:

cig1: 100S 25M
cig2: 2S 98M 25S

-> alignment ->

cig1: - 100S 25M
cig2: 2S 98M 25S
"""

# http://www.samformat.info/sam-format-flag

import re

from fuma.Fusion import STRAND_FORWARD, STRAND_REVERSE

pat_bam_parse_alignment_offset_using_cigar = re.compile("([0-9]+)([MIDNSHPX=])")


def cigar_to_cigartuple(cigar_str):
    """Converts a CIGAR string into a tuple compatible with Pysam.
    E.g. '10M5S' becomes: [(0, 10), (4, 5)]
    """
    tt = {'M': 0,  # BAM_CMATCH	0
          'I': 1,  # BAM_CINS	1
          'D': 2,  # BAM_CDEL	2
          'N': 3,  # BAM_CREF_SKIP	3
          'S': 4,  # BAM_CSOFT_CLIP	4
          'H': 5,  # BAM_CHARD_CLIP	5
          'P': 6,  # BAM_CPAD	6
          '=': 7,  # BAM_CEQUAL	7
          'X': 8}  # BAM_CDIFF	8

    cigartup = []

    for chunk in pat_bam_parse_alignment_offset_using_cigar.finditer(cigar_str):
        flag = chunk.group(2)
        length = chunk.group(1)
        cigartup.append((tt[flag], int(length)))

    return cigartup


class CigarAlignment:
    """
    STAR Fusion is not clear about which chunk is the one that started
    the read. We can calculate this back based on the CIGAR string.

    Usually you find a discordant pair as follows:

    part1: 46S80M
    part2: 80M46s

    Here the M and S are exectly complementary to each other. It does
    however occur that additional soft/hard clips are added. Hence we
    have to align them to each other and find the one with the least
    distance (in cases above r^2 = 0).

    We use matrices in a alignment context:

    tup1: 2S55M69S
    tup2: 57S69M

    Intiation:

                   2S     55M     69S
        ||======||=====||======||======||
        || 0    || 2^2 || 55^2 || 69^2 ||
        ||======||=====||======||======||
    57S || 57^2 ||      |       |       |
        ||======||------+-------+-------+
    69M || 69^2 ||      |       |       |
        ||======||------+-------+-------+

    fill function:
    i, j = min(
       (tup1[i] - tup2[j])   +   matrix[i-1, j-1]
       matrix[i-1]              , i >= 1
       matrix[j-1]              , j >= 1
    )

                   2S     55M     69S
        ||======||=====||======||======||
        || 0    || 2^2 || 55^2 || 69^2 ||
        ||======||=====||======||======||
    57S || 57^2 || i1   |  i2   |  i3   |
        ||======||------+-------+-------+
    69M || 69^2 || i2   |  i3   |  i4   |
        ||======||------+-------+-------+

    """
    def __init__(self, cigtup1, cigtup2):
        self.cigtup1 = self.cleanup_cigar(cigtup1, [1, 2, 3])
        self.cigtup2 = self.cleanup_cigar(cigtup2, [1, 2, 3])

        self.m = len(self.cigtup1)
        self.n = len(self.cigtup2)

        self.init_matrix()

    def cleanup_cigar(self, cigartup, invalid_chunks):
        clean = []

        # Removal of wrong chunks
        for chunk in cigartup:
            if chunk[0] not in invalid_chunks:
                clean.append(chunk)

        # Merging those that become adjacent to each other
        #   E.g. removal of 'N' in 6M10N6M -> 6M6M -> 12M
        # or:    [(4, 17), (0, 55), (1, 1), (0, 53)] -> [(4, 17), (0, 108)]

        last_type = -1
        concat = []

        for chunk in clean:
            if chunk[0] == last_type:

                # increase last insert
                concat[-1] = (chunk[0], concat[-1][1] + chunk[1])
            else:
                concat.append(chunk)
                last_type = chunk[0]

        return concat

    def init_matrix(self):
        # self.matrix = [[0] * (self.n + 1) ] * (self.m + 1) << copies references of lists >, <
        self.matrix = [[0 for i in xrange(self.n + 1)] for i in xrange(self.m + 1)]

        for i in xrange(1, self.m + 1):
            self.matrix[i][0] = pow(self.cigtup1[i - 1][1], 2)

        for j in xrange(1, self.n + 1):
            self.matrix[0][j] = pow(self.cigtup2[j - 1][1], 2)

        self.tb_matrix = [["?" for i in xrange(self.n + 1)] for i in xrange(self.m + 1)]
        self.tb_matrix[0][0] = 't'  # Terminate traceback; finished

        for i in xrange(1, self.m + 1):
            self.tb_matrix[i][0] = "|"

        for j in xrange(1, self.n + 1):
            self.tb_matrix[0][j] = "-"

    def str_tb_matrix(self):
        out = ''

        for j in xrange(self.n + 1):
            for i in xrange(self.m + 1):
                out += str(self.tb_matrix[i][j]) + "\t"
            out += "\n"

        return out

    def str_sc_matrix(self):
        out = ''

        for j in xrange(self.n + 1):
            for i in xrange(self.m + 1):
                out += str(self.matrix[i][j]) + "\t"
            out += "\n"

        return out

    def get_diagonal(self, diagonal):
            #  0, 0
            #  1, 0   0, 1
            #  2, 0   1, 1   0, 2

            i = diagonal
            j = 0

            for block in xrange(diagonal + 1):
                if i < self.m and j < self.n:
                    yield (i, j)

                i -= 1
                j += 1

    def cigar_diff(self, cig_chunk1, cig_chunk2):
        """Maybe sort on cigar type, to reduce if statements
        """
        if (cig_chunk1[0] == 0 and cig_chunk2[0] == 0) or (cig_chunk1[0] == 4 and cig_chunk2[0] == 4):
            #  M * M or S*S => sqaure
            return pow(cig_chunk1[1] + cig_chunk2[1], 2)

        elif (cig_chunk1[0] == 0 and cig_chunk2[0] == 4) or (cig_chunk1[0] == 4 and cig_chunk2[0] == 0):
            return pow((cig_chunk1[1] - cig_chunk2[1]), 2)

        else:
            raise Exception("Not yet implemented:", self.cigtup1, "\n", self.cigtup2, "\n\n", cig_chunk1[0], "\n", cig_chunk2[0])

    def calc_diff(self, i, j):
        """
                2S     55M     69S
     ||======||=====||======||======||
     || 0    || 2^2 || 55^2 || 69^2 ||
     ||======||=====||======||======||
 57S || 57^2 ||      |       |       |
     ||======||------+-------+-------+
 69M || 69^2 ||      |       |       |
     ||======||------+-------+-------+
"""
        # c_ins_1 = pow(self.cigtup1[i][1], 2)
        # c_ins_2 = pow(self.cigtup2[j][1], 2)
        c_ins_1 = self.matrix[i - 1][j] + pow(self.cigtup1[i][1], 2)  # Insertion vertical, in tup1
        c_ins_2 = self.matrix[i][j - 1] + pow(self.cigtup2[j][1], 2)  # Insertion horizontal, in tup2
        c_diag = self.matrix[i][j] + self.cigar_diff(self.cigtup1[i], self.cigtup2[j])

        if c_ins_1 < c_diag:
            _type = "|"
            _min = c_ins_1
        else:
            _type = "m"
            _min = c_diag

        if c_ins_2 < _min:
            _type = "-"
            _min = c_ins_2

        return (_min, _type)

    def fill_matrix(self):
        n_diagonals = (self.n + self.m) - 1
        for diagonal in xrange(n_diagonals):
            for cell in self.get_diagonal(diagonal):

                diff, _type = self.calc_diff(cell[0], cell[1])
                self.tb_matrix[cell[0] + 1][cell[1] + 1] = _type
                self.matrix[cell[0] + 1][cell[1] + 1] = diff

    def traceback_matrix(self):
        tup1_rev = []
        tup2_rev = []

        i = self.m + 1
        j = self.n + 1

        action = "s"
        while action != "t":
            if i < 0 and j < 0:
                raise Exception("Unpredicted out of bound")

            if action == "s":
                i -= 1
                j -= 1

            elif action == "m":
                tup1_rev.append(self.cigtup1[i - 1])
                tup2_rev.append(self.cigtup2[j - 1])

                i -= 1
                j -= 1

            elif action == "-":
                tup1_rev.append((-1, 0))
                tup2_rev.append(self.cigtup2[j - 1])

                j -= 1

            elif action == "|":
                tup1_rev.append(self.cigtup1[i - 1])
                tup2_rev.append((-1, 0))

                i -= 1

            action = self.tb_matrix[i][j]

        return tup1_rev[::-1], tup2_rev[::-1]

    def calculate_order(self, cigtup1_aligned, cigtup2_aligned):
        if len(cigtup1_aligned) != len(cigtup2_aligned):
            raise Exception("Lengths of aligned cigarstring tuples are different - error has taken place in the alignment", cigtup1_aligned, cigtup2_aligned)

        c1 = []
        c2 = []

        for i in xrange(len(cigtup1_aligned)):
            if cigtup1_aligned[i][0] != -1 and cigtup2_aligned[i][0] != -1:
                c1.append(cigtup1_aligned[i])
                c2.append(cigtup2_aligned[i])

        #  'M':0, #   	BAM_CMATCH	0
        #  'S':4, #   	BAM_CSOFT_CLIP	4

        offset = int(len(c1) / 2)
        c1_l = sum([x[1] for x in c1[0: offset] if x[0] == 0])
        c2_l = sum([x[1] for x in c2[0: offset] if x[0] == 0])

        c1_total = sum([x[1] for x in c1 if x[0] in [0, 4]])
        c2_total = sum([x[1] for x in c2 if x[0] in [0, 4]])

        if c1_total == 0 or c2_total == 0:
            # @todo figure out when this happens?
            return STRAND_FORWARD
        else:
            c1_r = 1.0 * c1_l / c1_total
            c2_r = 1.0 * c2_l / c2_total

            if c1_r >= c2_r:
                return STRAND_FORWARD
            else:
                return STRAND_REVERSE

    def get_order(self):
        #  1. align
        #    a. fill
        self.fill_matrix()
        #    b. trace back
        c1, c2 = self.traceback_matrix()

        #  2. calculate M / S only on those chunks that are aligned (M to S flags and vice versa)
        return self.calculate_order(c1, c2)
