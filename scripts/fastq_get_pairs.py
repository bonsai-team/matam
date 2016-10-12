#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fastq_get_pairs

Description: Retrieve paired and singleton reads from a single fastq file

  fastq_get_pairs.py -i input.fq

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-08-25
Last Modified: 2016-08-25
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


def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    # Variables initialization
    count = 0
    header = ''
    seq = ''
    qual = ''
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        count += 1
        if count % 4 == 1:
            if header:
                yield header, seq, qual
            header = line[1:]
        elif count % 4 == 2:
            seq = line
        elif count % 4 == 0:
            qual = line
    yield header, seq, qual
    # Close input file
    fastq_file_handle.close()


def buffer_paired_reads(fastq_fh):
    """
    """
    previous_read_id = ''
    read_buffer = list()
    # Reading each read in fastq file
    for header, seq, qual in read_fastq_file_handle(fastq_fh):
        #~ read_id = header.split()[0]
        read_id = header.split()[0][:-2]
        # Yield read buffer
        if read_id != previous_read_id:
            if previous_read_id:
                yield read_buffer
                read_buffer = list()
        # Append read into read buffer
        read = (header, seq, qual)
        read_buffer.append(read)
        # Store previous read id
        previous_read_id = read_id
    # Yield last read buffer
    yield read_buffer
    # Close input file
    fastq_fh.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Retrieve paired and singleton reads from a fastq file')
    parser.add_argument('-i', '--input_fastq',
                        metavar='INFQ',
                        type=argparse.FileType('r'),
                        default='-',
                        help='input fastq file')
    parser.add_argument('-b', '--basename',
                        metavar='BASENAME',
                        type=str,
                        help='Output files basename. '
                             'Default will be based on input file name')
    args = parser.parse_args()

    if not args.basename:
        input_basename = 'stdin'
        if args.input_fastq.name != '<stdin>':
            input_filename = os.path.basename(args.input_fastq.name)
            input_basename, input_extension = os.path.splitext(input_filename)
        args.basename = input_basename

    output_left_filename = args.basename + '.left' + input_extension
    out_left_fh = open(output_left_filename, 'w')

    output_right_filename = args.basename + '.right' + input_extension
    out_right_fh = open(output_right_filename, 'w')

    output_singleton_filename = args.basename + '.singleton' + input_extension
    out_single_fh = open(output_singleton_filename, 'w')

    previous_read_id = ''
    reads_buffer = list()

    for read_buffer in buffer_paired_reads(args.input_fastq):
        if len(read_buffer) == 2:
            read_left = read_buffer[0]
            out_left_fh.write('@{0}\n{1}\n+\n{2}\n'.format(*read_left))
            read_right = read_buffer[1]
            out_right_fh.write('@{0}\n{1}\n+\n{2}\n'.format(*read_right))
        elif len(read_buffer) == 1:
            read = read_buffer[0]
            out_single_fh.write('@{0}\n{1}\n+\n{2}\n'.format(*read))
