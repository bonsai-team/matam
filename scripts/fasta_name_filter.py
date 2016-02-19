#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FastaNameFilter

Description: Filter a fasta file based on a string to find in the
               sequences headers, or given a file with a list of id

  fastaNameFilter.py -i input.fa -o output.fa -s "stringtofind"
  fastaNameFilter.py -i input.fa -o output.fa -f sequencesnames.ids

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2014
Last Modified: 2016-01-13
Licence: GNU GPL 3.0

Copyright 2014-2016 Pierre Pericard

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
    for i in xrange(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Filter a fasta file based on sequence name.')
    parser.add_argument('-i', '--input_fasta', metavar='input', 
                        type=argparse.FileType('r', 0), default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta', metavar='output', 
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput fasta file')
    parser.add_argument('-s', '--stringtofind', metavar='string',
                        default='\&~"#\{\(\[-\|_\^@\)\]=\}', 
                        help='String to filter on')
    parser.add_argument('-f', '--fileids', metavar='file',
                        type=argparse.FileType('r'),
                        help='File with ids')
    args = parser.parse_args()
    
    ids_set = set()
    ids_list = list()
    
    if args.fileids:
        for line in args.fileids:
            ids_list.append(line.strip())
        ids_set = set(ids_list)
    
    tofind = re.compile(args.stringtofind, flags=re.IGNORECASE)
    
    for header, sequence in read_fasta_file_handle(args.input_fasta):
        seq_id = header.split()[0]
        if (tofind.search(header) or (seq_id in ids_set)):
            args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))

