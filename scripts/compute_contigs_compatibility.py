#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import time
import logging

# Create logger
logger = logging.getLogger(__name__)

# Get program filename, filepath and directory
program_filepath = os.path.realpath(sys.argv[0])
program_dirpath, program_filename = os.path.split(program_filepath)
program_name, program_extension = os.path.splitext(program_filename)


class SamAlignment():
    """
    Class for a SAM alignment
    Positions are all 0-based
    """

    def __init__(self, sam_tab):
        """
        Constructor from a list from a SAM line
        """
        # Parse SAM line
        self.query_id = sam_tab[0]
        self.flag = int(sam_tab[1])
        self.subject_id = sam_tab[2]
        self.first_pos = int(sam_tab[3]) # 1-based leftmost mapping position of the first matching base
        self.mapping_qual = int(sam_tab[4])
        self.cigar = sam_tab[5]
        self.rnext = sam_tab[6]
        self.pnext = sam_tab[7]
        self.tlen = sam_tab[8]
        self.query_seq = sam_tab[9].upper()
        self.qual = sam_tab[10]

        # Compute additional variables
        # Is read reverse complemented ?
        self.reverse_complemented = bool(self.flag & 0x10)  # Bitwise AND. True if flag == 16
        self.subject_start = self.first_pos - 1 # 0-based
        self.query_length = len(self.query_seq)

        # Parse CIGAR
        self.cigar_tab = self.parse_cigar(self.cigar)

        # Compute subject_end, query_start, query_end
        tmp_cigar_tab = self.cigar_tab[:]
        # Deal with soft clipping and get the query start (0-based)
        self.query_start = 0
        self.left_softclip_num = 0
        self.overhang_num = 0
        if tmp_cigar_tab[0][0] == 'S':
            count = tmp_cigar_tab[0][1]
            self.query_start += count
            self.left_softclip_num += count
            self.overhang_num += count
            del tmp_cigar_tab[0]
        # Parse CIGAR
        self.query_end = self.query_start - 1
        self.subject_end = self.subject_start - 1
        self.indel_num = 0
        self.matches_mismatches_num = 0
        self.right_softclip_num = 0
        for operation, count in tmp_cigar_tab:
            if operation == 'M':
                self.matches_mismatches_num += count
                self.subject_end += count
                self.query_end += count
            elif operation == 'I':
                self.query_end += count
                self.indel_num += count
            elif operation == 'D':
                self.subject_end += count
                self.indel_num += count
            elif operation == 'S':
                self.right_softclip_num += count
                self.overhang_num += count

    @staticmethod
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


def read_tab_file_handle_sorted(tab_file_handle, factor_index=0):
    """
    Parse a tab file (sorted by a column) and return a generator
    """
    previous_factor_id = ''
    factor_tab_list = list()
    # Reading tab file
    for tab in (l.split() for l in tab_file_handle if l.strip()):
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


def return_compatibility_status(alignment_i, alignment_j):
    """
    Return compatibility status between 2 SamAlignments:
      - 0 : no overlap
      - 1 : compatible overlap (with no error)
      - 2 : incompatible overlap (at least one error)
    """
    ov_start_on_ref = alignment_j.subject_start
    ov_end_on_ref = min(alignment_i.subject_end, alignment_j.subject_end)
    overlap_size = ov_end_on_ref - ov_start_on_ref

    if overlap_size <= 0:
        return 0

    # By default, alignments are incompatible
    compatibility_status = 2

    ov_start_on_sbj_i = alignment_i.query_start + (alignment_j.subject_start - alignment_i.subject_start)
    ov_end_on_sbj_i = ov_start_on_sbj_i + overlap_size

    ov_start_on_sbj_j = alignment_j.query_start
    ov_end_on_sbj_j = ov_start_on_sbj_j + overlap_size

    # TO DO: extend the overlap with the soft-clipped sequences

    ov_seq_i = alignment_i.query_seq[ov_start_on_sbj_i:ov_end_on_sbj_i+1]
    ov_seq_j = alignment_j.query_seq[ov_start_on_sbj_j:ov_end_on_sbj_j+1]

    # If the two contigs have an identical sequence on the overlap, then
    # the alignments are compatible
    if ov_seq_i == ov_seq_j:
        compatibility_status = 1

    return compatibility_status


def compute_contigs_compatibility_matrix(sam_alignments_list):
    """
    Given a list of SamAlignments, sorted by position on the ref,
    returns a compatibility matrix encoded as follow
      - 0 : no overlap
      - 1 : compatible overlap (with no error)
      - 2 : incompatible overlap (at least one error)
    """
    contigs_num = len(sam_alignments_list)
    contigs_compatibility_matrix = [[0 for x in range(contigs_num)] for y in range(contigs_num)]
    # Scan all contigs
    for i in range(contigs_num - 1):
        alignment_i = sam_alignments_list[i]
        # Scan all following contigs
        for j in range(i+1, contigs_num):
            alignment_j = sam_alignments_list[j]
            # Compute the compatibility status between the 2 SamAlignments
            compatibility_status = return_compatibility_status(alignment_i, alignment_j)
            # Stop when we know that contigs j wont overlap with
            # contig i anymore
            if compatibility_status == 0:
                break
            # Save the compatibility status in the matrix
            contigs_compatibility_matrix[i][j] = compatibility_status
    #
    return contigs_compatibility_matrix


def compute_bin_list(compatibility_matrix):
    """
    """
    bin_list = list()
    # Initialise the first bin
    bin_list.append([0])
    # Scan all contigs
    for j in range(1, len(compatibility_matrix)):
        #
        in_a_bin = False
        # Scan through all bins
        for b in bin_list:
            #
            overlap_with_b = False
            is_incompatible_with_b = False
            #
            for i in b:
                # Get compatibility status. i is always < j
                compatibility_status = compatibility_matrix[i][j]
                # If i overlaps with j
                if compatibility_status:
                    overlap_with_b = True
                    # Test for direct incompatibility
                    if compatibility_status == 2:
                        # There is at least one mismatch in the overlap between i and j
                        is_incompatible_with_b = True
                        # No need to test for other contigs in that bin
                        # we known it cannot belong to this one
                        break
                    elif compatibility_status == 1:
                        # i and j overlaps with no error
                        # Test for secondary compatibility
                        for x in range(j+1, len(compatibility_matrix)):
                            status_i = compatibility_matrix[i][x]
                            status_j = compatibility_matrix[j][x]
                            if status_i and status_j:
                                if status_i != status_j:
                                    is_incompatible_with_b = True
                                    break
                        if is_incompatible_with_b:
                            break
            if overlap_with_b and not is_incompatible_with_b:
                b.append(j)
                in_a_bin = True
                break
        #
        if not in_a_bin:
            # If SamAlignment i is not compatible with any existing bin,
            # then create a new bin
            bin_list.append([j])
    #
    return bin_list


class DefaultHelpParser(argparse.ArgumentParser):
    """
    This is a slightly modified argparse parser to display the full help
    on parser error instead of only usage
    """
    def error(self, message):
        sys.stderr.write('\nError: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def parse_arguments():
    """
    Parse the command line, and check if arguments are correct
    """
    # Initiate argument parser
    parser = DefaultHelpParser(description='Program description',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))

    # Main parameters
    group_main = parser.add_argument_group('Main parameters')
    # -i / --input_sam
    group_main.add_argument('-i', '--input_sam',
                            action = 'store',
                            metavar = 'INSAM',
                            type = argparse.FileType('r'),
                            default = '-',
                            help = 'Input sam file, sorted by subject and position')
    # -o / --output_sam
    group_main.add_argument('-o', '--output_sam',
                            action = 'store',
                            metavar = 'OUTSAM',
                            type = argparse.FileType('w'),
                            default = '-',
                            help = 'Output sam file')
    # -v / --verbose
    group_main.add_argument('-v', '--verbose',
                            action = 'store_true',
                            help = 'Increase verbosity')

    # Debug
    group_debug = parser.add_argument_group('Debug parameters')
    # --debug
    group_debug.add_argument('--debug',
                             action = 'store_true',
                             help = 'Output debug infos')

    args = parser.parse_args()

    #
    return args


if __name__ == '__main__':

    # Set global t0
    global_t0_wall = time.time()

    # Init global error code
    global_error_code = 0

    # Arguments parsing
    args = parse_arguments()

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

    #~ # Print intro infos
    #~ if args.verbose:
        #~ print_intro(args)

    #
    #~ sys.stderr.write('\n')

    ########
    # Start

    # Write header start
    args.output_sam.write('@HD\tVN:1.0\n')

    for tab_list in read_tab_file_handle_sorted(args.input_sam, 2):
        # Get reference id
        ref_id = tab_list[0][2]
        query_num = len(tab_list)

        # Parse sam file for one reference and generate a SamAlignment object for each line
        sam_alignments_list = list()
        for sam_tab in tab_list:
            sam_alignment = SamAlignment(sam_tab)
            #~ print(sam_alignment.__dict__)
            sam_alignments_list.append(sam_alignment)
        #~ print(sam_alignments_list)

        # Compute contigs compatibility matrix
        contigs_compatibility_matrix = compute_contigs_compatibility_matrix(sam_alignments_list)

        #~ for line in contigs_compatibility_matrix:
            #~ print(line)
        #~ print()

        # If there is no incompatible contigs, no need to split sam file
        if not any(2 in t for t in contigs_compatibility_matrix):
            args.output_sam.write('@SQ\tSN:{0}\tLN:9999999999\n'.format(ref_id))
            for tab in tab_list:
                args.output_sam.write('{0}\n'.format('\t'.join(tab)))
        # Else we compute sets of compatible contigs
        else:
            # Compute bins of compatible contigs
            bin_list = compute_bin_list(contigs_compatibility_matrix)
            #~ print(bin_list)
            # Output SAM lines
            for b in range(len(bin_list)):
                args.output_sam.write('@SQ\tSN:{0}_{1}\tLN:9999999999\n'.format(ref_id, b))
                for i in bin_list[b]:
                    tab = tab_list[i]
                    # Modify ref column
                    tab[2] += '_{0}'.format(b)
                    # Output SAM
                    args.output_sam.write('{0}\n'.format('\t'.join(tab)))

        #~ print()

    #######
    # Exit

    logger.info('-- Program complete --')

    if global_error_code > 0:
        logger.warning('Problems might have happened during program execution. Please check log above')
    else:
        logger.debug('Execution completed in {0:.2f} seconds'.format(time.time() - global_t0_wall))

    #~ sys.stderr.write('\n')
