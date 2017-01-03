#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import logging


# Create logger
logger = logging.getLogger(__name__)

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
                        metavar = 'INSAM',
                        type = argparse.FileType('r'),
                        default = '-',
                        help = 'Input sam file, sorted by reference and position')
    # -o / --output_cov
    parser.add_argument('-o', '--output_cov',
                        metavar = 'OUTCOV',
                        type = argparse.FileType('w'),
                        default = '-',
                        help = 'Output ref coverage tab file')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar = 'REF',
                        type = argparse.FileType('r'),
                        required = True,
                        help = 'References fasta file')
    # --histo
    parser.add_argument('--histo',
                        metavar = 'FILE',
                        type = str,
                        default = 'ref_coverage_histogram.svg',
                        help = 'References coverage histogram')
    # -v / --verbose
    parser.add_argument('-v', '--verbose',
                        action = 'store_true',
                        help = 'Increase verbosity')
    # --debug
    parser.add_argument('--debug',
                        action = 'store_true',
                        help = 'Debug mode')

    args = parser.parse_args()

    # Set logging
    # create console handler
    ch = logging.StreamHandler()
    #
    if args.debug:
        logger.setLevel(logging.DEBUG)
        # create formatter for debug level
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    else:
        if args.verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.WARNING)
        # create default formatter
        formatter = logging.Formatter('%(levelname)s - %(message)s')
    # add the formatter to the console handler
    ch.setFormatter(formatter)
    # add the handler to logger
    logger.addHandler(ch)

    # Get ref seqs and initialize positions count
    logger.info('Reading reference sequences')
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(args.references):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()

    #
    logger.info('Reading through input sam file')
    #
    coverage_list = list()
    frequency_count_dict = defaultdict(int)

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
        # Compute reference coverage
        ref_coverage = total_aligned_nt/ref_len
        #
        frequency_count_dict[ref_coverage] += 1
        coverage_list.append(ref_coverage)

    # Writing ref coverage tab
    logger.info('Writing output ref coverage tab file')
    sorted_ref_coverage_frequency_tuple_tuple = tuple(sorted(frequency_count_dict.items(), key=lambda x: x[0]))
    for frequency_tuple in sorted_ref_coverage_frequency_tuple_tuple:
        args.output_cov.write('{}\t{}\n'.format(*frequency_tuple))

    # Plotting histogram
    logger.info('Plotting histogram')
    plt.hist(coverage_list, range=(0,500), bins=250)
    plt.title("Ref Coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.ylim([0,200])
    plt.savefig(arg.histo)
    plt.show()