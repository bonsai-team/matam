#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get the length of every sequences from a fasta file.')
    parser.add_argument('-i', '--input_fasta',
                        metavar='INFA',
                        type=argparse.FileType('r'),
                        default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_tab',
                        metavar='OUTTAB',
                        type=argparse.FileType('w'),
                        default='-',
                        help='ouput tab file (seqid TAB length)')
    args = parser.parse_args()

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        seq_id = header.split()[0]
        args.output_tab.write("{0}\t{1}\n".format(seq_id, len(sequence)))
