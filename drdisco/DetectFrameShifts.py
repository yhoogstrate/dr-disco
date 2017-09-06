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

        self.gene_annotation_from_pre_coding = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        
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
        def index_gtf_transcripts(gtf_file, transcript_id):
            #gtf_file_entries = gtf_file[transcript_id]['cds']
            
            cumulative_offset = 0
            i = 0
            nxt = False
            prev = False
            previous = False

            # @todo move into load function
            exon_index = {}
            for feature in gtf_file[transcript_id]['cds']:
                exon_number = int(feature.attr['exon_number'])
                if exon_number in exon_index:
                    raise Exception("Error in GTF file - same exon id multiple times included: " + transcript_id + " - exon number: " + str(exon_number))
                else:
                    exon_index[exon_number] = feature

            start_codon = None
            for exon_number in sorted(exon_index.keys()):
                feature = exon_index[exon_number]

                if start_codon == None:
                    # This is the first exon, drop it.
                    if feature.iv.strand == '+':
                        # |::========>
                        start_codon = feature.iv.start
                    elif feature.iv.strand == '-':
                        # <========::|
                        start_codon = feature.iv.end
                

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
            
            
            
            
            # --- 
            exon_index = {}
            for feature in gtf_file[transcript_id]['exon']:
                exon_number = int(feature.attr['exon_number'])
                if exon_number in exon_index:
                    raise Exception("Error in GTF file - same exon id multiple times included: " + transcript_id + " - exon number: " + str(exon_number))
                else:
                    exon_index[exon_number] = feature
            exon_ids = sorted(exon_index.keys())
            
            i = 0
            n = len(exon_ids)
            while i < n:
                exon = exon_index[exon_ids[i]]
                add = True
                
                if start_codon != None:
                    # strand minus
                    # 
                    #         ::|       start_codon
                    #  <=============]  exon
                    #  <:::::::::       cds
                    #  |                splice site which produces coding region
                    # 
                    
                    # strand plus:
                    #       |::         start_codon
                    #  [=============>  exon
                    #       :::::::::>  cds

                    if exon.iv.strand == '-' and exon.iv.start <= start_codon:
                        add = False
                    elif exon.iv.strand == '+' and exon.iv.end >= start_codon:
                        add = False
                
                if add:
                    if exon.iv.strand == '-':
                        itv_from = HTSeq.GenomicInterval(chrom, exon.iv.start, exon.iv.start + 1, exon.iv.strand)
                    elif exon.iv.strand == '+':
                        itv_from = HTSeq.GenomicInterval(chrom, exon.iv.end, exon.iv.end + 1, exon.iv.strand)
                    
                    self.gene_annotation_from_pre_coding[itv_from] += transcript_id
                    i += 1
                
                else:# once the start codon is within the exons, skip scanning because translation has started
                    i = n
            
            print

        def load_gtf_per_transcript():
            transcript_idx = {}
            gtf_file = HTSeq.GFF_Reader(self.gtf_file, end_included=True)

            for feature in gtf_file:
                gtf_type = feature.type.lower()
                if gtf_type in ['cds', 'exon']:
                    transcript_id = feature.attr['gene_name'] + '(' + feature.attr['transcript_id'] + '.' + feature.attr['transcript_version'] + ')-' + feature.source

                    if transcript_id not in transcript_idx:
                        transcript_idx[transcript_id] = {
                            'exon': [],
                            'cds': []
                        }
                    transcript_idx[transcript_id][gtf_type].append(feature)

            return transcript_idx

        transcript_idx = load_gtf_per_transcript()
        for transcript_uid in sorted(transcript_idx.keys()):
            index_gtf_transcripts(transcript_idx, transcript_uid)

    def evaluate(self, _from, _to, offset):
        """
        Offset may be convenient because STAR sometimes has problems aligning/clipping the first 2 bases after an exon
        Values of 4 and larger do not make sense.
        """
        from_l_pre_coding = []
        
        from_l = []
        to_l = []

        if _from[0][0:3] == 'chr':
            _from[0] = _from[0][3:]

        if _to[0][0:3] == 'chr':
            _to[0] = _to[0][3:]

        
        for step in self.gene_annotation_from_pre_coding[HTSeq.GenomicInterval(_from[0], max(0, _from[1] - offset), _from[1] + offset + 1, _from[2])].steps():
            for entry in step[1]:
                from_l_pre_coding.append(entry)
        
        for step in self.gene_annotation_from[HTSeq.GenomicInterval(_from[0], max(0, _from[1] - offset), _from[1] + offset + 1, _from[2])].steps():
            for entry in step[1]:
                from_l.append(entry)

        for step in self.gene_annotation_to[HTSeq.GenomicInterval(_to[0], max(0, _to[1] - offset), _to[1] + offset + 1, _to[2])].steps():
            for entry in step[1]:
                to_l.append(entry)


        results = {0: [], 1: [], 2: [], 'pre_coding': []}

        for from_l_i in from_l_pre_coding:
            for to_l_i in to_l:
                results['pre_coding'].append((from_l_i, to_l_i))

        for from_l_i in from_l:
            for to_l_i in to_l:
                frame_shift = ((from_l_i[1] + to_l_i[1]) % 3)
                results[frame_shift].append((from_l_i, to_l_i))

        return results
