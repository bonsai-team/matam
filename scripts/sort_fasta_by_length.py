#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sort_fasta_by_length

Description: Sort a FastA file by sequence length (increasing by default)

  sort_fasta_by_length.py -i input.fa -o output.fa
  sort_fasta_by_length.py < input.fa > output.fa

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2015
Last Modified: 2016-04-13
Licence: GNU GPL 3.0

Copyright 2015-2016 Pierre Pericard

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

import argparse


def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in (l.strip() for l in fasta_file_handle if l.strip()):
        if line[0] == '>':
            # Yield the last read header and sequence
            if sequence_nb:
                yield (header, ''.join(seqlines))
                seqlines = list()
            # Get header
            header = line[1:]
            sequence_nb += 1
        else:
            # Concatenate sequence
            seqlines.append(line)
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

    # Initiate argument parser
    parser = argparse.ArgumentParser(description='Sort a fasta file by increasing sequence length.')

    # -i / --input_fasta
    parser.add_argument('-i', '--input_fasta',
                        action='store',
                        metavar='INFILE',
                        type=argparse.FileType('r'),
                        default='-',
                        help="Input fasta file. "
                             "Default is <stdin>")

    # -o / --output_fasta
    parser.add_argument('-o', '--output_fasta',
                        action='store',
                        metavar='OUTFILE',
                        type=argparse.FileType('w'),
                        default='-',
                        help="Ouput fasta file. "
                             "Default is <stdout>")

    # -r / --reverse
    parser.add_argument('-r', '--reverse',
                        action='store_true',
                        help='Sort by decreasing length')

    # Parse arguments from command line
    args = parser.parse_args()

    # Variables initialization
    seq_list = list()

    # Load all sequences
    for header, sequence in read_fasta_file_handle(args.input_fasta):
        seq_list.append((header, sequence))

    # Sort sequences by length
    if args.reverse:
        seq_list.sort(key=lambda x: len(x[1]), reverse=True)
    else:
        seq_list.sort(key=lambda x: len(x[1]))

    # Write sorted sequences to output file
    for header, sequence in seq_list:
        args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))

