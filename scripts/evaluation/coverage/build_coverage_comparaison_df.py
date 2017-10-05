#!/usr/bin/env python3

import os
import sys
import re
import random
import argparse
import math
import logging

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', '..')
sys.path.append(SCRIPTS_DIR)

from sample_sam_by_coverage import *

logger = logging.getLogger(__name__)


def write_rows(sam_handler, ref_seq_dict, tag, out_handler):
    for alignment_tabs_list in read_tab_file_handle_sorted((l for l in sam_handler if l[0]!='@'), 2):
        ref_id = alignment_tabs_list[0][2]
        ref_len = len(ref_seq_dict[ref_id])
        ref_coverage_list = compute_ref_coverage(alignment_tabs_list, ref_len)
        for pos, cov in enumerate(ref_coverage_list):
            row = (tag, ref_id, str(pos), str(cov))
            print('\t'.join(row), file=out_handler)


def make_dataframe(sam1_handler, sam2_handler, fasta_ref_handler, output_df_handler):
    dataframe = []
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(fasta_ref_handler):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()

    header = ('sampling', 'reference', 'position', 'coverage')
    print('\t'.join(header), file=output_df_handler)
    write_rows(sam1_handler, ref_seq_dict, "before", output_df_handler)
    write_rows(sam2_handler, ref_seq_dict, "after", output_df_handler)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # Arguments parsing
    parser = argparse.ArgumentParser(description='Build a table to compare the coverage of two SAM files')
    parser.add_argument('-s1', '--input_sam1',
                        type=argparse.FileType('r'),
                        required=True,
                        help='Non filtered sam, sorted by reference and position')
    parser.add_argument('-s2', '--input_sam2',
                        type=argparse.FileType('r'),
                        required=True,
                        help='filtered_sam, sorted by reference and position')
    parser.add_argument('-r', '--references',
                        type=argparse.FileType('r'),
                        required=True,
                        help='References fasta file')
    parser.add_argument('-o', '--output_df',
                        type=argparse.FileType('w'),
                        help="Tabulated file: header=('sampling', 'reference', 'position', 'coverage')",
                        default='-')


    args = parser.parse_args()
    logger.info('Generate comparaison table')
    make_dataframe(args.input_sam1, args.input_sam2, args.references, args.output_df)
    logger.info('Done')


