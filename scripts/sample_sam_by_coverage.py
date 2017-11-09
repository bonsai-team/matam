#!/usr/bin/env python3

import sys
import re
import random
import argparse
import numpy as np


def tab_list_group_by(tab_list, factor_index=0):
    """
    Read a sorted tab list and group by the given colum index
    """
    previous_factor_id = ''
    factor_tab_list = list()
    # Reading tab file
    for tab in tab_list:
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


def read_tab_file_handle_sorted(tab_file_handle, factor_index=0):
    """
    Parse a tab file (sorted by a column) and return a generator
    """
    for factor_tab_list in tab_list_group_by((l.split() for l in tab_file_handle if l.strip()), factor_index):
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


def get_alignment_length_on_ref(cigar):
    """
    """
    return sum((count for operation, count in parse_cigar(cigar) if operation in ('M','D','N','=','X')))


def compute_ref_coverage(ref_sam_tab_list, ref_length):
    """
    """
    # Init ref coverage list
    ref_coverage_list = np.zeros(ref_length, dtype=np.int)
    # Read alignments on the ref
    for alignment_tab in ref_sam_tab_list:
        starting_pos = int(alignment_tab[3]) - 1 #SAM positions are 1-based
        cigar = alignment_tab[5]
        alignment_length_on_ref = get_alignment_length_on_ref(cigar)
        end_pos = starting_pos + alignment_length_on_ref - 1 #We want the last mapped nucleotide, not the first free one
        ref_coverage_list[starting_pos:end_pos+1] += 1
    #
    return ref_coverage_list


def sample_by_depth(sam_handler, fasta_ref_handler, threshold, sampled_out_sam_handler):
    """
    sample the SAM file based on the coverage at each pos
    """
    # Read reference sequences and store them in a dict
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(fasta_ref_handler):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()
    # Read SAM by ref
    # OPTI: instead of storing all SAM records for each ref in mem, we could read twice the file block corresponding to this ref
    for ref_sam_tab_list in read_tab_file_handle_sorted((l for l in sam_handler if l[0]!='@'), 2):
        ref_id = ref_sam_tab_list[0][2]
        ref_length = len(ref_seq_dict[ref_id])
        # Compute ref coverage list
        ref_coverage_list = compute_ref_coverage(ref_sam_tab_list, ref_length)
        #
        for pos_sam_tab_list in tab_list_group_by(ref_sam_tab_list, 3):
            starting_pos = int(pos_sam_tab_list[0][3]) - 1 #SAM positions are 1-based
            ref_coverage = ref_coverage_list[starting_pos]
            if ref_coverage <= threshold:
                for alignment_tab in pos_sam_tab_list:
                    print('{}'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler)
            else:
                num_alignments_to_remove = ref_coverage - threshold
                random_alignments_order_indices = [i for i in range(len(pos_sam_tab_list))]
                random.shuffle(random_alignments_order_indices)
                #
                for i in random_alignments_order_indices:
                    alignment_tab = pos_sam_tab_list[i]
                    if num_alignments_to_remove > 0:
                        starting_pos = int(alignment_tab[3]) - 1 #SAM positions are 1-based
                        cigar = alignment_tab[5]
                        alignment_length_on_ref = get_alignment_length_on_ref(cigar)
                        end_pos = starting_pos + alignment_length_on_ref - 1 #We want the last mapped nucleotide, not the first free one
                        if min(ref_coverage_list[starting_pos:end_pos+1]) <= threshold:
                            print('{}'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler)
                        else:
                            ref_coverage_list[starting_pos:end_pos+1] -= 1
                            num_alignments_to_remove -= 1
                    else:
                        print('{}'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler)


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
                        type=int,
                        default=50,
                        help='Identity threshold. '
                             'Default is %(default)s')
    #
    args = parser.parse_args()
    #
    random.seed()
    sample_by_depth(args.input_sam, args.references, args.cov_threshold, args.output_sam)