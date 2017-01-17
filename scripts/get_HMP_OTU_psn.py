#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
get_HMP_OTU_psn?py

Description: Get OTU sequences represented in a given HMP PSN sample

  get_HMP_OTU_psn.py -i rep_set_v13.fna -s 700015979 -t otu_table_psn_v13.txt -o SRS011405.QIIME.OTU_seq.fa

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2017-01-17
Last Modified: 2017-01-17
Licence: GNU GPL 3.0

Copyright 2017 Pierre Pericard

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
"""


import sys
import os
import argparse
import string
import re

def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in fasta_file_handle:
        if line[0] == '>':
            # Yield the last read header and sequence
            if sequence_nb:
                yield (header, ''.join(seqlines))
                del seqlines[:]
            # Get header
            header = line[1:].rstrip()
            sequence_nb += 1
        else:
            # Concatenate sequence
            seqlines.append(line.strip())
    # Yield the input file last sequence
    yield (header, ''.join(seqlines))
    # Close input file
    fasta_file_handle.close()

def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get OTU sequences represented in a given HMP PSN sample')
    parser.add_argument('-i', '--input_fasta',
                        metavar='PATH',
                        action = 'store',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input OTU fasta file')
    parser.add_argument('-o', '--output_fasta',
                        metavar='PATH',
                        action = 'store',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output selected OTUs fasta file')
    parser.add_argument('-s', '--sample',
                        metavar='ID',
                        action = 'store',
                        type=str,
                        required = True,
                        help='PNS sample id')
    parser.add_argument('-t', '--otu_table',
                        metavar='PATH',
                        action = 'store',
                        type=argparse.FileType('r'),
                        required = True,
                        help='File with OTU-sample correspondance')
    args = parser.parse_args()


    #
    args.otu_table.readline() # Deal with the first line

    samples_tuple = tuple(args.otu_table.readline().strip().split('\t'))
    sys.stderr.write('INFO: {} samples ID are in the OTU table\n'.format(len(samples_tuple)))

    # Get sample index
    sample_index = samples_tuple.index(args.sample)
    sys.stderr.write('DEBUG: Sample ID was found at index #{}\n'.format(sample_index))

    # Read the rest of the otu table and get the OTU IDs present in the sample
    sample_otus_ids_list = list()
    for tab in (l.strip().split('\t') for l in args.otu_table if l.strip()):
        if tab[sample_index] == '1':
            sample_otus_ids_list.append(tab[0])

    sample_otus_ids_frozenset = frozenset(sample_otus_ids_list)

    # Read the input fasta file and filter
    for header, sequence in read_fasta_file_handle(args.input_fasta):
        seq_id = header.split()[0]
        if seq_id in sample_otus_ids_frozenset:
            args.output_fasta.write('>{}\n{}\n'.format(header, format_seq(sequence, 80)))


