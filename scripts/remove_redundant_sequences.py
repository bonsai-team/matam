#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
remove_redundant_sequences.py

Description: Remove sequences entirely included in bigger sequences

  remove_redundant_sequences.py -i input.fa -o output.fa

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016
Last Modified: 2016-06-14
Licence: GNU GPL 3.0

Copyright 2016 Pierre Pericard

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

    parser = argparse.ArgumentParser(description='Remove redundant sequences.')
    parser.add_argument('-i', '--input_fasta',
                        metavar='input',
                        type=argparse.FileType('r'),
                        default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta',
                        metavar='output',
                        type=argparse.FileType('w'),
                        default='-',
                        help='ouput fasta file')
    args = parser.parse_args()

    sequences_list = [(h, s) for h, s in read_fasta_file_handle(args.input_fasta)]

    sequences_list.sort(key=lambda x: len(x[1]))

    sequences_to_keep_list = list()

    for i in range(len(sequences_list)-1):
        short_seq = sequences_list[i][1]
        #~ print(short_seq)
        to_keep = True
        for j in range(i+1, len(sequences_list)):
            long_seq = sequences_list[j][1]
            if short_seq in long_seq:
                to_keep = False
                break
        if to_keep:
            sequences_to_keep_list.append(sequences_list[i])
    sequences_to_keep_list.append(sequences_list[-1])

    for header, seq in sequences_to_keep_list:
        args.output_fasta.write('>{0}\n{1}\n'.format(header, format_seq(seq)))
