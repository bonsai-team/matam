#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FastQ to FastA

Description : Convert a FastQ file to a FastA file
  
  # Running examples:
  
  fastq_to_fasta.py -i input.fastq -o output.fasta
  fastq_to_fasta.py < input.fastq > output.fasta

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-04-13
Last Modified: 2016-04-13
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

import argparse


def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    
    # Variables initialization
    line_count = 0
    header = ''
    seq = ''
    qual = ''
    
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        line_count += 1
        if line_count % 4 == 1:
            if header:
                # Yield previous sequence
                yield header, seq, qual
            # Get read complete header (sequence id and description)
            header = line[1:]
        elif line_count % 4 == 2:
            # Get read sequence
            seq = line
        elif line_count % 4 == 0:
            # Get read quality
            qual = line
    
    # Yield the last sequence
    yield header, seq, qual
    
    # Close input file
    fastq_file_handle.close()


def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


if __name__ == '__main__':
    
    # Argument parser initialization
    parser = argparse.ArgumentParser(description='Convert a fastq file to a fasta file')
    
    # -i / --input_fastq
    parser.add_argument('-i', '--input_fastq',
                        action='store',
                        metavar='INFASTQ', 
                        type=argparse.FileType('r'),
                        default='-',
                        help="Input fastq file. "
                             "Default to <stdin>")
    
    # -o / --output_fasta
    parser.add_argument('-o', '--output_fasta',
                        action='store',
                        metavar='OUTFASTA', 
                        type=argparse.FileType('w'),
                        default='-',
                        help="Output fasta file. "
                             "Default to <stdout>")
    
    # Parse arguments from command line
    args = parser.parse_args()
    
    
    # Read fastq file and write fasta sequences
    for header, seq, qual in read_fastq_file_handle(args.input_fastq):
        args.output_fasta.write(">{0}\n".format(header))
        args.output_fasta.write("{0}\n".format(format_seq(seq, 80)))
    
