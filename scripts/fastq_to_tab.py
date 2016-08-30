#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fastq_to_tab

Description: Convert a fastq file to a tab file

  fastq_to_tab.py -i input.fq -o output.tab

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-04-12
Last Modified: 2016-04-12
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert a fastq file to a tab file')
    parser.add_argument('-i', '--input_fastq', metavar='INFASTQ',
                        type=argparse.FileType('r', 0), default='-',
                        help='input fastq file')
    parser.add_argument('-o', '--output_tab', metavar='OUTTAB',
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput tab file')
    args = parser.parse_args()

    for header, seq, qual in read_fastq_file_handle(args.input_fastq):
        args.output_tab.write('{0}\t{1}\t{2}\n'.format(header, seq, qual))
