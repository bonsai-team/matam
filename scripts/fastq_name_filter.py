#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fastq_name_filter

Description: Filter a fastq file based on a string to find in the
               sequences headers

  fastq_name_filter.py -i input.fq -o output.fq -s "stringtofind"

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
            header = line[1:].split()[0]
        elif count % 4 == 2:
            seq = line
        elif count % 4 == 0:
            qual = line
    # yield last fastq sequence
    yield header, seq, qual
    # Close input file
    fastq_file_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter a fastq file based on sequence name.')
    parser.add_argument('-i', '--input_fastq', metavar='input',
                        type=argparse.FileType('r'), default='-',
                        help='input fastq file')
    parser.add_argument('-o', '--output_fastq', metavar='output',
                        type=argparse.FileType('w'), default='-',
                        help='ouput fastq file')
    parser.add_argument('-s', '--stringtofind', metavar='string',
                        type=str, help='String to filter on')
    parser.add_argument('-f', '--fileids', metavar='file',
                        type=argparse.FileType('r'),
                        help='File with ids')
    args = parser.parse_args()

    if not args.stringtofind and not args.fileids:
        parser.print_help()
        raise Exception('Either a string or an id file has to be supplied')

    if args.fileids:
        ids_list = list()
        for line in args.fileids:
            ids_list.append(line.strip())
        ids_set = frozenset(ids_list)
        for header, seq, qual in read_fastq_file_handle(args.input_fastq):
            read_id = header.split()[0]
            if read_id in ids_set:
                args.output_fastq.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
    else:
        tofind = re.compile(args.stringtofind, flags=re.IGNORECASE)
        for header, seq, qual in read_fastq_file_handle(args.input_fastq):
            if tofind.search(header):
                args.output_fastq.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))

