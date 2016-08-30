#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

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

    parser = argparse.ArgumentParser(description='Filter a fasta file based on sequence name.')
    parser.add_argument('-i', '--input_fasta', metavar='input',
                        type=argparse.FileType('r'), required=True,
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta', metavar='output',
                        type=argparse.FileType('w'), required=True,
                        help='ouput fasta file')
    args = parser.parse_args()

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        sequence_id = header.split()[0]
        args.output_fasta.write(">{0}\n{1}\n".format(sequence_id, format_seq(sequence)))

