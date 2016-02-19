#! /usr/bin/python -u
# -*- coding: utf-8 -*-

"""

"""

import sys
import os
import argparse
import string
import re


def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in xrange(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Convert a tabulated sequence file in a fasta file.')
    parser.add_argument('-i', '--input_tab', metavar='input', 
                        type=argparse.FileType('r', 0), default='-',
                        help='input tab file')
    parser.add_argument('-o', '--output_fasta', metavar='output',
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput fasta file')
    parser.add_argument('--header', metavar='INT',
                        type=int, default=1,
                        help='Column of the header (first one by default)')
    parser.add_argument('--sequence', metavar='INT',
                        type=int, default=2,
                        help='Column of the sequence (second one by default)')
    args = parser.parse_args()
    
    # Convert column num to indices
    args.header -= 1
    args.sequence -= 1
    
    # Read tab and output fasta
    for line in args.input_tab:
        l = line.strip()
        if l:
            tab = l.split('\t')
            header = tab[args.header]
            sequence = tab[args.sequence]
            args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))




