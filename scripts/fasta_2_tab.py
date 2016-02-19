#! /usr/bin/python -u
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


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Convert a fasta file in a tabulated file.')
    parser.add_argument('-i', '--input_fasta', metavar='input', 
                        type=argparse.FileType('r', 0), default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_tab', metavar='output',
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput tab file')
    parser.add_argument('--length', action='store_true',
                        help='Add a column with the sequence length')
    args = parser.parse_args()
    
    for header, sequence in read_fasta_file_handle(args.input_fasta):
        header = re.sub(r'\s+', ' ', header)
        args.output_tab.write('{0}\t{1}'.format(header, sequence))
        if args.length:
            args.output_tab.write('\t{2}'.format(len(sequence)))
        args.output_tab.write('\n')


