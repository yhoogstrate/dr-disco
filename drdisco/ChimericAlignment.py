#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

import os

import pysam

from fuma.Fusion import STRAND_FORWARD
from .CigarAlignment import CigarAlignment

from drdisco import __version__
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

"""
 Sep 01: fix bam files:
 - SA tag at the end to link the multiple supplementary alignments IF present - by using samtools -n
 - Fix the mate of the other:
  [     A   >                [  B ]   [  C >
  A -> B
  B -> C
  C <- B
 - Update NH tags per mate
 - Mark as deletion if possible (N is for splicing, so do D/P)
 - Set read group to sample name
"""


class ChimericAlignment:
    def __init__(self, input_alignment_file):
        self.input_alignment_file = input_alignment_file
        self.test_pysam_version()

    def test_pysam_version(self):
        versions = pysam.__version__.split('.', 2)

        if int(versions[1]) < 9:
            raise Exception("Version of pysam needs to be at least 0.9 but is: " + pysam.__version__ + " instead")
        else:
            return True

    def set_read_group(self, all_reads_updated, group):
        for a in all_reads_updated:
            a.set_tag('RG', group)

    def fix_alignment_score(self, all_reads_updated):
        for a in all_reads_updated:
            a.is_paired = True
            a.is_read2 = True

    def set_qname_to_group(self, all_reads_updated):
        qnames = []
        for a in all_reads_updated:
            qnames.append(a.qname)

        qnames = list(set(qnames))

        if len(qnames) != 1:
            raise Exception("Not all reads belong to the same QNAME")
        else:
            qname = qnames[0]
            for i in xrange(len(all_reads_updated)):
                all_reads_updated[i].set_tag('LB', qname.replace(":", "."))

        return all_reads_updated

    def update_sa_tags(self, reads_updated, bam_file):
        sa_ids = {}

        for read in reads_updated:
            try:
                nm = read.get_tag('nM')
            except:
                nm = -1

            sa_id = [bam_file.get_reference_name(read.reference_id),
                     read.reference_start,
                     read.cigarstring,
                     read.mapping_quality, "-" if read.is_reverse else "+",
                     nm]

            # -- pysam 0.10.0 SPECIFIC --
            # Use of id() is some kind of hack because the objects seem
            # to get lost in memory after the set_tag?
            # i suppose set_tag does not change the object but generates
            # a new one
            sa_ids[id(read)] = ",".join([str(x) for x in sa_id])
        for read in reads_updated:
            reads = [k for k in reads_updated if k != read]

            sa_tags = [sa_ids[id(l)] for l in reads]
            sa_tag = ';'.join(sa_tags)
            read.set_tag('SA', sa_tag)

        return reads_updated

    def set_next_ref(self, aligned_segment, position):
        if aligned_segment.next_reference_id != position[0] or aligned_segment.next_reference_start != position[1]:
            a = pysam.AlignedSegment()
            a.query_name = aligned_segment.query_name
            a.query_sequence = aligned_segment.query_sequence
            a.flag = aligned_segment.flag
            a.reference_id = aligned_segment.reference_id
            a.reference_start = aligned_segment.reference_start
            a.mapping_quality = aligned_segment.mapping_quality
            a.cigartuples = aligned_segment.cigartuples
            a.template_length = aligned_segment.template_length
            a.query_qualities = aligned_segment.query_qualities
            a.set_tags(aligned_segment.get_tags())

            a.next_reference_id = position[0]
            a.next_reference_start = position[1]

            return a

        else:
            return aligned_segment

    def get_closest(self, location, alignments):
        d_chr = None
        d_pos = None
        closest = None

        for alignment in alignments:
            is_closest = False
            if (d_chr is None and d_pos is None):
                is_closest = True
            else:
                dd_chr = abs(d_chr - location[0])
                dd_pos = abs(d_pos - location[1])

                if (dd_chr < d_chr) or (dd_chr == d_chr and dd_pos < d_pos):
                    #  Ref to itself is ackward
                    # if dd_chr == 0 and dd_pos == 0:
                    #    pass
                    # else:
                    d_chr = dd_chr
                    d_pos = dd_pos

                    is_closest = True

            if is_closest:
                closest = alignment

        return closest

    def get_closest_by_hi(self, hi_closest, alignments):
        d_hi = None
        closest = None

        for alignment in alignments:
            is_closest = False

            cur_hi = alignment.get_tag('HI')
            dd_hi = abs(hi_closest - cur_hi)

            if d_hi is None:
                is_closest = True
            else:
                if dd_hi < d_hi:
                    is_closest = True

            if is_closest:
                d_hi = dd_hi
                closest = alignment

        return closest

    def fix_chain(self, alignments, bam_file, mates):
        """
        all segments in variable `aligments` must come from the same MATE (unless discordant pairs without suppl. alignments) and should not be marked as DUPLICATE or MULTIMAPPING
        """
        chains_from = {}
        chains_to = {}

        k = len(alignments)

        for alignment in alignments:
            chain_from = str(alignment.reference_id) + ":" + str(alignment.reference_start)
            chain_to = str(alignment.next_reference_id) + ":" + str(alignment.next_reference_start)

            if chain_from not in chains_from.keys():
                chains_from[chain_from] = []
            if chain_to not in chains_to.keys():
                chains_to[chain_to] = []

            chains_from[chain_from].append(alignment)
            chains_to[chain_to].append(alignment)

        chains_from = set(chains_from)
        chains_to = set(chains_to)

        # _from = chains_from.difference(chains_to)
        # _to = chains_to.difference(chains_from)
        _linked = list(chains_to.intersection(chains_from))

        # situaties:
        # 1. net zoveel chains_to als chains_from -> alleen unieke start posities
        #    a. alle chains_from en chains_to zijn identiek, op 1 in beide lijsten na
        #       [a, b, c, d, e, f, g, h]
        #         [b, c, d, e, f, g, h, i]
        #       In dit geval is de chain al perfect
        #    b. de overlap van chains heeft meer dan 1 mismatch per vector - er is iets
        #       mis met de mapping van de next reads. Dit klopt niet en spuug een error
        #       uit
        #        [a, b, c, c, d, e]
        #          [b, c, c, d, e, f]
        # 2. Er zijn chains met dubbele entries in de TO en dus ook minder die overeen
        #    komen met de chains from. Updaten moet door die to's aan te passen. De
        #    entrie die het dichtste bij de 'to' zit moet er naar blijven linken, de
        #    andere(n) naar elkaar - per chromosoom?

        new_alignments = []
        new_mates = []

        if len(mates) == 1 and len(alignments) >= 2:
            if mates[0].is_read1:
                hi_closest = -1
            else:
                hi_closest = len(alignments) + 1

            start = self.get_closest_by_hi(hi_closest, alignments)
            last_pos = [mates[0].reference_id, mates[0].reference_start]

            alignments = [a for a in alignments if a != start]
            # If the mate is not exactly matched but close, fix it:
            new_mates.append(self.set_next_ref(mates[0], [start.reference_id, start.reference_start]))

            hi_closest = start.get_tag('HI')
            i = 0
            while len(alignments) >= 1:
                closest = self.get_closest_by_hi(hi_closest, alignments)
                hi_closest = closest.get_tag('HI')

                next_pos = [closest.reference_id, closest.reference_start]

                s_fixed = self.set_next_ref(start, next_pos)
                s_fixed.set_tag('FI', i)
                new_alignments.append(s_fixed)
                alignments = [a for a in alignments if a != closest]

                start = closest
                i += 1

            # Map last one back to the mate again
            if len(alignments) == 0:
                start = self.set_next_ref(start, last_pos)
                start.set_tag('FI', i)
                new_alignments.append(start)

            # - Does it based on shortest genomic distance:
            # if str(mates[0].next_reference_id) + ":" + str(mates[0].next_reference_start) in chains_from:
                # next_pos = [mates[0].next_reference_id, mates[0].next_reference_start]
                # last_pos = [mates[0].reference_id, mates[0].reference_start]

                # start = self.get_closest(next_pos, alignments)

            # else:
                # logger.warn("Warning - mates do not correspond? - maybe empty (-1) as well?")

                # next_pos = [alignments[0].reference_id, alignments[0].reference_start]
                # last_pos = [mates[0].reference_id, mates[0].reference_start]

                # start = self.get_closest(next_pos, alignments)

            # alignments = [a for a in alignments if a != start]
            # # If the mate is not exactly matched but close, fix it:
            # new_mates.append(self.set_next_ref(mates[0],[start.reference_id, start.reference_start]))

            # i = 0
            # while len(alignments) >= 1:
                # closest = self.get_closest(next_pos, alignments)
                # next_pos = [closest.reference_id, closest.reference_start]

                # s_fixed = self.set_next_ref(start, next_pos)
                # s_fixed.set_tag('FI', i)
                # new_alignments.append(s_fixed)
                # alignments = [a for a in alignments if a != closest]

                # start = closest
                # i + = 1

            # # Map last one back to the mate again
            # if len(alignments) == 0:
                # start = self.set_next_ref(start, last_pos)
                # start.set_tag('FI', i)
                # new_alignments.append(start)

            # Now do it based on HI-tag

        elif len(mates) == 0:
            # Either 2 discordant mates
            # Or 2 discordant segments from one singleton

            if len(_linked) == len(alignments):
                # cross reffing each other - is already fine
                return alignments, new_mates

            seg_pos = None
            if len(_linked) > 0:
                seg_pos = [int(x) for x in list(_linked)[0].split(":")]

            if seg_pos is None or (seg_pos[0] == -1 and seg_pos[1] == -1):
                seg_pos = [alignments[0].reference_id, alignments[0].reference_start]

            closest = self.get_closest(seg_pos, alignments)
            alignments = [a for a in alignments if a != closest]
            new_alignments.append(self.set_next_ref(closest, seg_pos))

            while len(alignments) > 0:
                seg_pos = [closest.reference_id, closest.reference_start]
                closest = self.get_closest(seg_pos, alignments)

                alignments = [a for a in alignments if a != closest]
                new_alignments.append(self.set_next_ref(closest, seg_pos))

        else:
            raise Exception("Dunno how to handle junctions in both mates yet... not aware of STAR Fusion producing them either")

        if len(new_alignments) != k:
            raise Exception("Somewhere alignments got lost in this function")

        return new_alignments, new_mates

    def reconstruct_alignments(self, alignments, bam_file, fh_out):
        """
        input: only reads with the same qname
        """
        n = len(alignments)
        r1 = []
        r2 = []
        singletons = []

        for alignment in alignments:
            if alignment.is_read1:
                r1.append(alignment)
            elif alignment.is_read2:
                r2.append(alignment)
            else:
                singletons.append(alignment)

        n_r1 = len(r1)
        n_r2 = len(r2)
        n_s = len(singletons)

        all_reads_updated = []

        if n_r1 == 2:
            reads_updated, mates_updated = self.fix_chain(r1, bam_file, r2)

            ca = CigarAlignment(reads_updated[0].cigar, reads_updated[1].cigar)
            if ca.get_order() == STRAND_FORWARD:
                """These reads have the opposite strand because they are both read1
                """
                self.set_read_group([reads_updated[0]], 'spanning_paired_1_s')
                self.set_read_group([reads_updated[1]], 'spanning_paired_2_s')
            else:
                self.set_read_group([reads_updated[0]], 'spanning_paired_1_t')
                self.set_read_group([reads_updated[1]], 'spanning_paired_2_t')

            self.set_read_group(mates_updated, 'silent_mate')
            for a in reads_updated:
                all_reads_updated.append(a)

            for a in mates_updated:
                all_reads_updated.append(a)

        elif n_r2 == 2:
            reads_updated, mates_updated = self.fix_chain(r2, bam_file, r1)

            ca = CigarAlignment(reads_updated[0].cigar, reads_updated[1].cigar)
            if ca.get_order() == STRAND_FORWARD:
                self.set_read_group([reads_updated[0]], 'spanning_paired_1')
                self.set_read_group([reads_updated[1]], 'spanning_paired_2')
            else:
                # Tested in 'tests/fix-chimeric/test_terg_03.filtered.bam'
                self.set_read_group([reads_updated[0]], 'spanning_paired_2')
                self.set_read_group([reads_updated[1]], 'spanning_paired_1')

            self.set_read_group(mates_updated, 'silent_mate')
            for a in reads_updated:
                all_reads_updated.append(a)

            for a in mates_updated:
                all_reads_updated.append(a)

        elif n_s == 2:
            reads_updated, mates_updated = self.fix_chain(singletons, bam_file, [])

            """ deprecated
            ca = CigarAlignment(reads_updated[0].cigar, reads_updated[1].cigar)
            if ca.get_order() == STRAND_FORWARD:
                if reads_updated[0].get_tag('HI') == 2 and reads_updated[1].get_tag('HI') == 1:
                    self.set_read_group([reads_updated[0]],'spanning_singleton_1r')
                    self.set_read_group([reads_updated[1]],'spanning_singleton_2')
                elif reads_updated[0].get_tag('HI') == 1 and reads_updated[1].get_tag('HI') == 2:
                    self.set_read_group([reads_updated[0]],'spanning_singleton_2r')
                    self.set_read_group([reads_updated[1]],'spanning_singleton_1')
                else:
                    raise Exception("Unknown strand order for singletons: %s\n%s\n", reads_updated[0], reads_updated[1])
            """

            if reads_updated[0].get_tag('HI') == 2 and reads_updated[1].get_tag('HI') == 1:
                self.set_read_group([reads_updated[0]], 'spanning_singleton_1_r')
                self.set_read_group([reads_updated[1]], 'spanning_singleton_2_r')
            elif reads_updated[0].get_tag('HI') == 1 and reads_updated[1].get_tag('HI') == 2:
                self.set_read_group([reads_updated[0]], 'spanning_singleton_2')
                self.set_read_group([reads_updated[1]], 'spanning_singleton_1')
            else:
                raise Exception("Unknown strand order for singletons: %s (%i)\n%s (%i)\n",
                                reads_updated[0].query_name,
                                reads_updated[0].reference_start,
                                reads_updated[1].query_name,
                                reads_updated[1].reference_start)

            self.fix_alignment_score(reads_updated)

            for a in reads_updated:
                all_reads_updated.append(a)

            for a in mates_updated:
                all_reads_updated.append(a)

        elif n_r1 == 1 and n_r2 == 1 and n_s == 0:
            reads_updated, mates_updated = self.fix_chain(r1 + r2, bam_file, [])
            self.set_read_group(reads_updated, 'discordant_mates')

            for a in reads_updated:
                all_reads_updated.append(a)

            for a in mates_updated:
                all_reads_updated.append(a)

        else:
            if n == 1:
                log.warn("segments of mate are missing: " + alignments[0].query_name)
                all_reads_updated.append(alignments[0])
            else:
                raise Exception("what happens here?")

        all_reads_updated = self.update_sa_tags(all_reads_updated, bam_file)
        if len(all_reads_updated) != n:
            raise Exception("Error - reads have been lost")

        all_reads_updated = self.set_qname_to_group(all_reads_updated)

        for a in all_reads_updated:
            fh_out.write(a)

    def convert(self, bam_file_discordant_fixed, temp_dir):
        basename, ext = os.path.splitext(os.path.basename(self.input_alignment_file))
        basename = temp_dir.rstrip("/") + "/" + basename

        # @TODO / consider todo - start straight from sam
        # samtools view -bS samples/7046-004-041_discordant.Chimeric.out.sam > samples/7046-004-041_discordant.Chimeric.out.unsorted.bam

        log.info("Convert into a name-sorted bam file, to get all reads with the same name adjacent to each other")
        pysam.sort("-o", basename + ".name-sorted.bam", "-n", self.input_alignment_file)

        log.info("Fixing sam file")
        sam_file_discordant = pysam.AlignmentFile(basename + ".name-sorted.bam", "rb")
        header = sam_file_discordant.header
        header['RG'] = [
            {'ID': 'discordant_mates', 'DS': 'This read has discordant mate pair'},
            {'ID': 'silent_mate', 'DS': 'Reads of this type are not discordant while their mate is'},
            {'ID': 'spanning_paired_1', 'DS': 'This read was aligned to two locations and also has an aligned mate'},
            {'ID': 'spanning_paired_1_r', 'DS': 'This read was aligned to two locations and also has an aligned mate (strand type r)'},
            {'ID': 'spanning_paired_1_s', 'DS': 'This read was aligned to two locations and also has an aligned mate (strand type s)'},
            {'ID': 'spanning_paired_1_t', 'DS': 'This read was aligned to two locations and also has an aligned mate (strand type t)'},
            {'ID': 'spanning_paired_2', 'DS': 'This read was aligned to two locations and also has an aligned mate'},
            {'ID': 'spanning_paired_2_r', 'DS': 'This read was aligned to two locations and also has an aligned mate (strand type r)'},
            {'ID': 'spanning_paired_2_s', 'DS': 'This read was aligned to two locations and also has an aligned mate (strand type s)'},
            {'ID': 'spanning_paired_2_t', 'DS': 'This read was aligned to two locations and also has an aligned mate (strand type t)'},
            {'ID': 'spanning_singleton_1', 'DS': 'This read was aligned to two locations but no aligned mate'},
            {'ID': 'spanning_singleton_1_r', 'DS': 'This read was aligned to two locations but no aligned mate'},
            {'ID': 'spanning_singleton_2', 'DS': 'This read was aligned to two locations but no aligned mate'},
            {'ID': 'spanning_singleton_2_r', 'DS': 'This read was aligned to two locations but no aligned mate'}]

        header['PG'] = [
            {'ID': 'drdisco_fix_chimeric', 'PN': 'drdisco fix-chimeric', 'CL': '', 'VN': __version__}
        ]

        fh = pysam.AlignmentFile(basename + ".name-sorted.fixed.sam", "wb", header=header)
        last_read_name = False
        alignments = []
        for read in sam_file_discordant:
            if read.qname != last_read_name:
                if len(alignments) > 0:
                    self.reconstruct_alignments(alignments, sam_file_discordant, fh)
                alignments = []
                last_read_name = read.qname
            alignments.append(read)
        self.reconstruct_alignments(alignments, sam_file_discordant, fh)
        fh.close()

        log.info("Converting fixed file into BAM")
        fhq = open(basename + ".name-sorted.fixed.bam", "wb")
        fhq.write(pysam.view('-bS', basename + ".name-sorted.fixed.sam"))
        fhq.close()

        log.info("Sorting position based fixed file")
        pysam.sort("-o", basename + ".sorted.fixed.bam", basename + ".name-sorted.fixed.bam")

        log.info("Indexing the position sorted bam file")
        pysam.index(basename + ".sorted.fixed.bam")

        log.info("Cleaning up temp files")
        for fname in [basename + ".name-sorted.bam", basename + ".name-sorted.fixed.sam", basename + ".name-sorted.fixed.bam"]:
            log.debug("=> " + fname)
            os.remove(fname)

        log.info("Moving to final destination")
        os.rename(basename + ".sorted.fixed.bam", bam_file_discordant_fixed)
        os.rename(basename + ".sorted.fixed.bam" + ".bai", bam_file_discordant_fixed + ".bai")
