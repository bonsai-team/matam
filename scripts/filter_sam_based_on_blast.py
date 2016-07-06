#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_sam
    parser.add_argument('-i', '--input_sam',
                        metavar='INSAM',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input sam file')
    # -b / --input_blast
    parser.add_argument('-b', '--input_blast',
                        metavar='INBLAST',
                        type=argparse.FileType('r'),
                        required=True,
                        help='Input blast tab file. ')
    # -o / --output_sam
    parser.add_argument('-o', '--output_sam',
                        metavar='OUTSAM',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Ouput filtered sam file')

    args = parser.parse_args()

    query_subject_tuple_list = list()

    # Reading blast tab file
    for tab in (l.split() for l in args.input_blast if l.split()):

        # Parse blast tab
        query_id = tab[0]
        subject_id = tab[1]

        # Store query-subject tuple
        query_subject_tuple_list.append((query_id, subject_id))

    # Convert query-subject tuple list in frozenset
    query_subject_tuple_frozenset = frozenset(query_subject_tuple_list)

    # Read and filter sam input
    for line in (l for l in args.input_sam if l.split()):
        if line[0] != '@':
            tab = line.split()
            query_id = tab[0]
            subject_id = tab[2]

            if (query_id, subject_id) in query_subject_tuple_frozenset:
                args.output_sam.write(line)
