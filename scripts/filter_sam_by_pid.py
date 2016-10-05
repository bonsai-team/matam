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
                        help='Input sam file')
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
    # -t / --id_threshold
    parser.add_argument('-t', '--id_threshold',
                        metavar='ID',
                        type=float,
                        default=0.99,
                        help='Identity threshold')
    # -s / --soft_clip_as_mismatch
    parser.add_argument('-s', '--soft_clip_as_mismatch',
                        action = 'store_true',
                        help = 'Soft-clipped nucleotides are considered '
                               'as mismatches')

    args = parser.parse_args()

    # Get ref seqs and initialize positions count
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(args.references):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()

    # Reading sam file
    for tab in (l.split() for l in args.input_sam if l.strip()):
        # Parse sam tab
        query_id = tab[0]
        flag = int(tab[1])
        subject_id = tab[2]
        first_pos = int(tab[3]) # 1-based leftmost mapping position of the first matching base
        mapping_qual = int(tab[4])
        cigar = tab[5]
        rnext = tab[6]
        pnext = tab[7]
        tlen = tab[8]
        query_seq = tab[9].upper()
        qual = tab[10]

        # Compute additional variables
        # Is read reverse complemented ?
        reverse_complemented = flag & 0x10  # Bitwise AND. True if flag == 16
        subject_seq = ref_seq_dict[subject_id]
        subject_start = first_pos - 1 # 0-based
        query_length = len(query_seq)

        cigar_tab = parse_cigar(cigar)

        # Init variables
        query_start = 0
        overhang_num = 0
        indel_num = 0
        matches_num = 0
        mismatches_num = 0

        # Deal with soft clipping and get the query start (0-based)
        if cigar_tab[0][0] == 'S':
            count = cigar_tab[0][1]
            query_start += count
            overhang_num += count
            if args.soft_clip_as_mismatch:
                mismatches_num += count
            del cigar_tab[0]

        # Parse CIGAR
        query_end = query_start - 1
        subject_end = subject_start - 1
        for operation, count in cigar_tab:
            if operation == 'M':
                # Compute the number of matches on this block
                local_matches_num = sum((query_seq[query_end + 1 + i] == subject_seq[subject_end + 1 + i] for i in range(0, count)))
                matches_num += local_matches_num
                mismatches_num += count - local_matches_num
                subject_end += count
                query_end += count
            elif operation == 'I':
                query_end += count
                indel_num += count
            elif operation == 'D':
                subject_end += count
                indel_num += count
            elif operation == 'S':
                overhang_num += count
                mismatches_num += count

        identity = float(matches_num) / (matches_num + mismatches_num + indel_num)

        if identity >= args.id_threshold:
            args.output_sam.write('{0}\n'.format('\t'.join(tab)))
