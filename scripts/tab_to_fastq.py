#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
tab_to_fastq

Description: Convert a tab file to a fastq file

  tab_to_fastq.py -i input.tab -o output.fq

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert a tab file to a fastq file')
    parser.add_argument('-i', '--input_tab', metavar='INTAB',
                        type=argparse.FileType('r', 0), default='-',
                        help='input tab file')
    parser.add_argument('-o', '--output_fastq', metavar='OUTFASTQ',
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput fastq file')
    args = parser.parse_args()

    for tab in (l.strip().split('\t') for l in args.input_tab if l.strip()):
        header= tab[0]
        seq = tab[1]
        qual = tab[2]
        args.output_fastq.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
