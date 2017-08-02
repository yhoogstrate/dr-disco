#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


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


class DetectFrameShifts:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file

        self.gene_annotation_from = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.gene_annotation_to = HTSeq.GenomicArrayOfSets("auto", stranded=True)

        self.index_gtf()

    def index_gtf(self):
        """
        GTF file must have:
        CDS entries for coding sequences
        each CDS entry must have:
            - source
            - gene_name attribute
            - transcript_id attribute
            - transcript_version attribute
            - exon_number attribute

        Such gtf files are provided by Ensembl
        """
        def index_gtf_transcripts(gtf_file_entries, transcript_id):
            cumulative_offset = 0
            i = 0
            nxt = False
            prev = False
            previous = False

            # @todo move into load function
            exon_index = {}
            for feature in gtf_file_entries:
                exon_number = int(feature.attr['exon_number'])
                if exon_number in exon_index:
                    raise Exception("Error in GTF file - same exon id multiple times included: " + transcript_id + " - exon number: " + str(exon_number))
                else:
                    exon_index[exon_number] = feature

            for exon_number in sorted(exon_index.keys()):
                feature = exon_index[exon_number]

                length = (feature.iv.end) - feature.iv.start + cumulative_offset
                off1 = length % 3
                off2 = -length % 3

                if feature.iv.chrom[0:3] == 'chr':
                    chrom = feature.iv.chrom[3:]
                else:
                    chrom = feature.iv.chrom

                if nxt:
                    if feature.iv.strand == '+':
                        itv_from = HTSeq.GenomicInterval(chrom, previous.iv.end, previous.iv.end + 1, previous.iv.strand)
                        itv_to = HTSeq.GenomicInterval(chrom, feature.iv.start, feature.iv.start + 1, feature.iv.strand)

                    elif feature.iv.strand == '-':
                        itv_from = HTSeq.GenomicInterval(chrom, previous.iv.start, previous.iv.start + 1, previous.iv.strand)
                        itv_to = HTSeq.GenomicInterval(chrom, feature.iv.end, feature.iv.end + 1, feature.iv.strand)

                    self.gene_annotation_from[itv_from] += (transcript_id, int(prev))
                    self.gene_annotation_to[itv_to] += (transcript_id, int(nxt))

                prev = str(off1)
                nxt = str(off2)
                previous = feature

                cumulative_offset = off1

                i += 1

        def load_gtf_per_transcript():
            transcript_idx = {}
            gtf_file = HTSeq.GFF_Reader(self.gtf_file, end_included=True)

            for feature in gtf_file:
                if feature.type == 'CDS':
                    transcript_id = feature.attr['gene_name'] + '(' + feature.attr['transcript_id'] + '.' + feature.attr['transcript_version'] + ')-' + feature.source

                    if transcript_id not in transcript_idx:
                        transcript_idx[transcript_id] = []
                    transcript_idx[transcript_id].append(feature)

            return transcript_idx

        transcript_idx = load_gtf_per_transcript()
        for transcript_uid in sorted(transcript_idx.keys()):
            index_gtf_transcripts(transcript_idx[transcript_uid], transcript_uid)

    def evaluate(self, _from, _to, offset):
        """
        Offset may be convenient because STAR sometimes has problems aligning/clipping the first 2 bases after an exon
        Values of 4 and larger do not make sense.
        """
        from_l = []
        to_l = []

        if _from[0][0:3] == 'chr':
            _from[0] = _from[0][3:]

        if _to[0][0:3] == 'chr':
            _to[0] = _to[0][3:]

        for step in self.gene_annotation_from[HTSeq.GenomicInterval(_from[0], max(0, _from[1] - offset), _from[1] + offset + 1, _from[2])].steps():
            for entry in step[1]:
                from_l.append(entry)

        for step in self.gene_annotation_to[HTSeq.GenomicInterval(_to[0], max(0, _to[1] - offset), _to[1] + offset + 1, _to[2])].steps():
            for entry in step[1]:
                to_l.append(entry)

        results = {0: [], 1: [], 2: []}

        for from_l_i in from_l:
            for to_l_i in to_l:
                frame_shift = ((from_l_i[1] + to_l_i[1]) % 3)
                results[frame_shift].append((from_l_i, to_l_i))

        return results
