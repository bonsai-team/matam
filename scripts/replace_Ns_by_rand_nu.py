#! /usr/bin/python -u
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
    for i in xrange(0, len(seq), linereturn):
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
    args = parser.parse_args()
    
    ACGT = ['A', 'C', 'G', 'T']
    
    for header, sequence in read_fasta_file_handle(args.input_fasta):
        sequence = re.sub(r'[^ACGT]', lambda x: random.choice(ACGT), sequence.upper())
        args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))


