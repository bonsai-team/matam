#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse


def read_tab_file_handle_sorted(tab_file_handle, factor_index=0):
    """
    Parse a tab file (sorted by a column) and return a generator
    """
    previous_factor_id = ''
    factor_tab_list = list()
    # Reading tab file
    for line in tab_file_handle:
        l = line.strip()
        if l:
            tab = l.split()
            current_factor = tab[factor_index]
            # Yield the previous factor tab list
            if current_factor != previous_factor_id:
                if previous_factor_id:
                    yield factor_tab_list
                    factor_tab_list = list()
            factor_tab_list.append(tab)
            previous_factor_id = current_factor
    # Yield the last tab list
    yield factor_tab_list
    # Close tab file
    tab_file_handle.close()


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


def parse_cigar(cigar):
    """
    Parse a CIGAR string and return a list of (operation, count) tuples
    """
    cigar_tab = list()
    count_str = ''
    for c in cigar:
        if c.isdigit():
            count_str += c
        else:
            operation = c
            count = 1
            if count_str:
                count = int(count_str)
            if cigar_tab and operation == cigar_tab[-1][0]:
                cigar_tab[-1][1] += count
            else:
                cigar_tab.append((operation, count))
            count_str = ''
    #
    return cigar_tab


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_sam
    parser.add_argument('-i', '--input_sam',
                        metavar='INSAM',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input sam file, sorted by reference and position')
    # -o / --output_sam
    parser.add_argument('-o', '--output_sam',
                        metavar='OUTSAM',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output filtered sam file')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar='REF',
                        type=argparse.FileType('r'),
                        required=True,
                        help='References fasta file')
    # -c / --cov_threshold
    parser.add_argument('-c', '--cov_threshold',
                        metavar='COV',
                        type=float,
                        default=50,
                        help='Identity threshold. '
                             'Default is %(default)s')

    args = parser.parse_args()

    # Get ref seqs and initialize positions count
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(args.references):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()

    # Reading sam file reference by reference
    for alignment_tabs_list in read_tab_file_handle_sorted(args.input_sam, 2):
        ref_id = alignment_tabs_list[0][2]
        ref_len = len(ref_seq_dict[ref_id])
        total_aligned_nt = 0
        for alignment_tab in alignment_tabs_list:
            cigar = alignment_tab[5]
            for operation, count in parse_cigar(cigar):
                if operation == 'M':
                    total_aligned_nt += count
        ref_coverage = total_aligned_nt/ref_len
        #~ sys.stderr.write('{}\t{}\n'.format(ref_id, ref_coverage))
        if ref_coverage <= args.cov_threshold:
            for alignment_tab in alignment_tabs_list:
                args.output_sam.write('{}\n'.format('\t'.join(alignment_tab)))
