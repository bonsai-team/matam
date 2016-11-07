#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import sys
import os
import argparse
import string
import re
import random

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

    parser = argparse.ArgumentParser(description='Replace all the Ns by random nucleotides.')
    parser.add_argument('-i', '--input_fasta', metavar='input',
                        type=argparse.FileType('r'), default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta', metavar='output',
                        type=argparse.FileType('w'), default='-',
                        help='ouput fasta file')
    parser.add_argument('-n', '--max_consec_n',
                        metavar = 'MAXN',
                        type = int,
                        default = 5,
                        help = 'Maximum number of consecutive Ns allowed in the sequence. '
                               'Default is %(default)s')
    args = parser.parse_args()

    rejected_seq_num = 0

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        n_matches = re.findall(r'[^ACGT]+', sequence.upper())
        max_length_n = 0
        if n_matches:
            max_length_n = max((len(m) for m in n_matches))
        if max_length_n <= args.max_consec_n:
            sequence = re.sub(r'[^ACGT]', 'A', sequence.upper())
            args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
        else:
            rejected_seq_num += 1
            #~ print('\nSeq is rejected ({} max consecutive Ns)'.format(max_length_n))
            #~ print(">{0}\n{1}\n".format(header, format_seq(sequence)))
    
    sys.stderr.write('{} sequences were rejected\n'.format(rejected_seq_num))
