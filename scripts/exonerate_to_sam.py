#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
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


def reverse_complement(sequence):
    """
    Reverse complement a DNA sequence
    """
    # Define the reverse complement dict
    rev_comp_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                     'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    # Perform reverse complementation
    rev_comp_seq = ''.join([rev_comp_dict[nt] for nt in sequence][::-1])
    #
    return rev_comp_seq


def tab_to_sam(query_id, subject_id, reverse_complemented,
               query_start, query_end, subject_start, subject_end,
               exo_cigar, query_seq):
    """
    """
    sam_tab = list()

    # QNAME
    sam_tab.append(query_id)

    # FLAG
    flag = 0
    if reverse_complemented:
        flag = 16
    sam_tab.append(str(flag))

    # RNAME
    sam_tab.append(subject_id)

    # POS # 1-based leftmost ref position of the first mapped read nucl
    sam_tab.append(str(subject_start + 1))

    # MAPQ
    sam_tab.append('255')

    # CIGAR
    cigar = str()

    if query_start > 0:
        cigar += '{0}S'.format(query_start)

    exo_cigar_tab = exo_cigar.split()

    for i in range(0, len(exo_cigar_tab), 2):
        label = exo_cigar_tab[i]
        label_count = int(exo_cigar_tab[i+1])
        cigar += '{0}{1}'.format(label_count, label)

    if query_end < len(query_seq):
        cigar += '{0}S'.format(len(query_seq) - query_end)

    sam_tab.append(cigar)

    # RNEXT
    sam_tab.append('*')

    # PNEXT
    sam_tab.append('0')

    # TLEN
    sam_tab.append('0')

    # SEQ
    sam_tab.append(query_seq)

    # QUAL
    sam_tab.append('*')

    #
    return '\t'.join(sam_tab)


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_exo
    parser.add_argument('-i', '--input_exo',
                        metavar='INEXO',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input exonerate tab file')
    # -o / --output_sam
    parser.add_argument('-o', '--output_sam',
                        metavar='OUTSAM',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output sam file')
    # -q / --queries
    parser.add_argument('-q', '--queries',
                        metavar='QUERIES',
                        type=argparse.FileType('r'),
                        required=True,
                        help='Queries fasta file')

    args = parser.parse_args()

    # Get queries sequences
    queries_seq_dict = dict()
    for header, seq in read_fasta_file_handle(args.queries):
        seqid = header.split()[0]
        queries_seq_dict[seqid] = seq

    # Reading exonerate tab file
    for tab in (l.split('\t') for l in args.input_exo if l.strip()):

        # Parse exonerate tab
        query_id = tab[0]
        subject_id = tab[1]

        # Exonerate uses an in-between coordinate system.
        # see: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual
        query_start = int(tab[2])
        query_end = int(tab[3])
        subject_start = int(tab[4])
        subject_end = int(tab[5])

        exo_cigar = tab[6]

        # Get query seq
        query_seq = queries_seq_dict[query_id]
        query_len = len(query_seq)

        # Deal with reverse complementation
        reverse_complemented = False
        if query_end < query_start:
            reverse_complemented = True
            # Reverse query start and end
            query_start = query_len - query_start
            query_end = query_len - query_end
            # Reverse complement seq
            query_seq = reverse_complement(query_seq)

        # Write sam output
        sam_line = tab_to_sam(query_id, subject_id, reverse_complemented,
                              query_start, query_end,
                              subject_start, subject_end,
                              exo_cigar, query_seq)

        args.output_sam.write('{0}\n'.format(sam_line))
