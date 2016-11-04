#! /usr/bin/python -u
# -*- coding: utf-8 -*-

"""

"""

import sys
import os
import argparse
import string
import re

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Restore the .1 and .2 names to a fastq file missing them (e.g. when dowloaded with SRA toolkit)')
    parser.add_argument('-i', '--input_fastq', metavar='input',
                        type=argparse.FileType('r'),
                        default='-',
                        help='input interleaved fastq file')
    parser.add_argument('-o', '--output_fastq', metavar='output',
                        type=argparse.FileType('w'),
                        default='-',
                        help='ouput fastq file')
    args = parser.parse_args()

    count = 0
    read_num = 0

    for line in args.input_fastq:
        line = line.strip()
        if line:
            count += 1
            if count % 4 == 1:
                read_num += 1
                header_tab = line[1:].split()

                seqid = header_tab[0]
                parity = ((read_num + 1) % 2) + 1

                args.output_fastq.write('@{0}.{1}'.format(seqid, parity))

                if len(header_tab) > 1:
                    description_line = ' '.join(header_tab[1:])
                    args.output_fastq.write(' {0}'.format(description_line))

                args.output_fastq.write('\n')

            elif count % 4 == 3:
                args.output_fastq.write('+\n')
            else:
                args.output_fastq.write('{0}\n'.format(line))



