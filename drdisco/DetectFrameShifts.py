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


class DetectFrameShifts:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file

        # fgd = full gene dysregulation (junction falls before coding sequences in both genes)
        # valid:
        # ===    ====        |     ====    ==[::]   CDS starts in middle of 2nd exon in acceptor gene: (from) pre-coding -> (to) pre-coding
        # ===    ====        |     ==[::]           CDS starts in middle of 1st exon in acceptor gene: (from) pre-coding -> (to) pre-coding
        # ===    ====        |     ====    [::::]   CDS starts at start of 2nd exon in acceptor gene:  (from) pre-coding -> (to) first coding
        # ===    ====        |     [::::]           CDS starts at start of 1st exon in acceptor gene:  (from) pre-coding -> (to) first coding

        self.gene_annotation_from_fgd = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        self.gene_annotation_to_fgd = HTSeq.GenomicArrayOfSets("auto", stranded=True)

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

        def load_gtf_per_transcript():
            transcript_idx = {}
            gtf_file = HTSeq.GFF_Reader(self.gtf_file, end_included=True)

            for feature in gtf_file:
                gtf_type = feature.type.lower()
                if gtf_type in ['cds', 'exon']:
                    try:
                        transcript_id = feature.attr['gene_name'] + '(' + feature.attr['transcript_id'] + '.' + feature.attr['transcript_version'] + ')-' + feature.source

                        if transcript_id not in transcript_idx:
                            transcript_idx[transcript_id] = {}

                        exon_number = int(feature.attr['exon_number'])
                        if exon_number not in transcript_idx[transcript_id]:
                            transcript_idx[transcript_id][exon_number] = {'exon': None, 'cds': None}

                        if gtf_type in ['exon', 'cds']:
                            transcript_idx[transcript_id][exon_number][gtf_type] = feature

                    except KeyError:
                        log.warn("Warning: GTF file misses certain attributes (gene_name, transcript_id or transcript_version) and is therefore skipping the frameshift detection.")
                        # there is no GFF_Reader.close() so break it dirty:
                        break

            return transcript_idx

        def insert_transcript_idx(transcript_idx):
            def clean_chrom(chrom):
                if chrom[0:3] == 'chr':
                    return chrom[3:]
                else:
                    return chrom

            def calc_from(feature):
                if feature.iv.strand == '+':
                    return HTSeq.GenomicInterval(clean_chrom(feature.iv.chrom), feature.iv.end, feature.iv.end + 1, feature.iv.strand)
                elif feature.iv.strand == '-':
                    return HTSeq.GenomicInterval(clean_chrom(feature.iv.chrom), feature.iv.start, feature.iv.start + 1, feature.iv.strand)

            def calc_to(feature):
                if feature.iv.strand == '+':
                    return HTSeq.GenomicInterval(clean_chrom(feature.iv.chrom), feature.iv.start, feature.iv.start + 1, feature.iv.strand)
                elif feature.iv.strand == '-':
                    return HTSeq.GenomicInterval(clean_chrom(feature.iv.chrom), feature.iv.end, feature.iv.end + 1, feature.iv.strand)

            """
            @todo change to:

            for exon in exons:
                if there is a CDS with similar exon-id:

                    if there is no stop codon or this exon id is the last exon id:
                        last_coding_exon = exon_id

                    if there is a start_codon with the same exon-id:
                        first_coding_exon = exon_id


                    - see if this is pre-coding
                        - add to pre-coding 'from' idx
                        - add to pre-coding 'to' list - als je naar dit exon fuseert moet transcriptie nog starten
                    - see if this is first coding exon
                        - add to normal 'from' list
                        - add to pre-coding 'to' list
                    - see if this is an inbetween coding exon
                        - add to normal 'from' list
                        - add to normal 'to' list
                    - see if this is the last coding exon
                        - add to normal 'to' list
                    """
            for transcript_id in transcript_idx:
                coding = "pre"

                cumulative_offset = 0
                exon_ids = sorted(transcript_idx[transcript_id].keys())

                for e in exon_ids:
                    exon = transcript_idx[transcript_id][e]

                    if coding == "pre":
                        if exon['cds'] is None:  # - pre coding
                            # distances are not relevant
                            self.gene_annotation_to_fgd[calc_to(exon['exon'])] += (transcript_id)
                            self.gene_annotation_from_fgd[calc_from(exon['exon'])] += (transcript_id)
                        else:  # - first coding
                            length = (exon['cds'].iv.end - exon['cds'].iv.start) + cumulative_offset

                            off1 = length % 3
                            off2 = -length % 3

                            self.gene_annotation_from[calc_from(exon['exon'])] += (transcript_id, off1)
                            self.gene_annotation_to_fgd[calc_to(exon['exon'])] += (transcript_id)

                            cumulative_offset = off1
                            coding = True
                    elif coding is True:
                        if e == exon_ids[-1] or transcript_idx[transcript_id][e + 1]['cds'] is None:  # - last coding
                            self.gene_annotation_to[calc_to(exon['exon'])] += (transcript_id, off2)

                            coding = "post"
                        else:  # - middle coding
                            self.gene_annotation_to[calc_to(exon['exon'])] += (transcript_id, off2)

                            length = (exon['cds'].iv.end - exon['cds'].iv.start) + cumulative_offset
                            off1 = length % 3
                            off2 = -length % 3

                            self.gene_annotation_from[calc_from(exon['exon'])] += (transcript_id, off1)

                            cumulative_offset = off1
                    #  else: # - post coding

        transcript_idx = load_gtf_per_transcript()
        insert_transcript_idx(transcript_idx)

    def evaluate(self, _from, _to, offset):
        """
        Offset may be convenient because STAR sometimes has problems aligning/clipping the first 2 bases after an exon
        Values of 4 and larger do not make sense.
        """
        from_l_fgd = []
        to_l_fgd = []

        from_l = []
        to_l = []

        if _from[0][0:3] == 'chr':
            _from[0] = _from[0][3:]

        if _to[0][0:3] == 'chr':
            _to[0] = _to[0][3:]

        # These are ends of the exons in which translation was not (yet?) started
        # A full gene disregulation is found when in the next gene a start_codon-containing exon is found
        # these still need to be implemented
        for step in self.gene_annotation_from_fgd[HTSeq.GenomicInterval(_from[0], max(0, _from[1] - offset), _from[1] + offset + 1, _from[2])].steps():
            for entry in step[1]:
                from_l_fgd.append(entry)

        for step in self.gene_annotation_to_fgd[HTSeq.GenomicInterval(_to[0], max(0, _to[1] - offset), _to[1] + offset + 1, _to[2])].steps():
            for entry in step[1]:
                to_l_fgd.append(entry)

        for step in self.gene_annotation_from[HTSeq.GenomicInterval(_from[0], max(0, _from[1] - offset), _from[1] + offset + 1, _from[2])].steps():
            for entry in step[1]:
                from_l.append(entry)

        for step in self.gene_annotation_to[HTSeq.GenomicInterval(_to[0], max(0, _to[1] - offset), _to[1] + offset + 1, _to[2])].steps():
            for entry in step[1]:
                to_l.append(entry)

        results = {0: [], 1: [], 2: [], 'fgd': []}

        for from_l_i in from_l_fgd:
            for to_l_i in to_l_fgd:
                results['fgd'].append((from_l_i, to_l_i))

        for from_l_i in from_l:
            for to_l_i in to_l:
                frame_shift = ((from_l_i[1] + to_l_i[1]) % 3)
                results[frame_shift].append((from_l_i, to_l_i))

        return results
