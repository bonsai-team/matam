#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re

if __name__ == '__main__':

    # Arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_sam', metavar='input.sam',
                        type=argparse.FileType('r', 0), default='-',
                        help='input sam file')
    parser.add_argument('-o', '--output_sam', metavar='output.sam',
                        type=argparse.FileType('w', 0), default='-',
                        help='output sam file')
    parser.add_argument('--geometric', action='store_true',
                        help='Uses an adaptative geometric threshold (straigth treshold by default)')
    parser.add_argument('-t', '--threshold', metavar='PCT',
                        type=float, default=0.9,
                        help='Threshold will be either threshold*max score (straigth mode) or threshold*previous score (geometric mode)')
    args = parser.parse_args()

    #
    blank = re.compile(r'^\s*$')
    score_pattern = re.compile(r"\s+AS:i:(\d+)\s+")
    previous_read_id = ""
    previous_score = 0
    threshold = 0.0
    to_write = True

    for line in args.input_sam:
        if not blank.match(line):
            tab = line.split()
            read_id = tab[0]
            score = score_pattern.search(line).group(1)
            if read_id != previous_read_id:
                to_write = True
                if not args.geometric:
                    threshold = float(score) * args.threshold
            else:
                if args.geometric:
                    threshold = float(previous_score) * (args.threshold)
                if (float(score) < threshold):
                    to_write = False

            previous_read_id = read_id
            previous_score = score

            #~ args.output_sam.write('{0}\t{1}\t{2}\t{3}\n'.format(read_id, score, threshold, to_write))

            if to_write:
                args.output_sam.write("{0}\n".format(line.strip()))

