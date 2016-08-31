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
                        required=True,
                        help='Input sam file, sorted by subject. '
                             'Only best alignments are expected.')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar='REF',
                        type=argparse.FileType('r'),
                        required=True,
                        help='References fasta file (optional). '
                             'If used, the blast file should only '
                             'contain the best alignments')

    args = parser.parse_args()

    # Get ref seqs and initialize positions count
    ref_seq_dict = dict()
    ref_positions_count_dict = dict()
    ref_metrics_dict = dict()
    total_ref_length = 0
    for header, seq in read_fasta_file_handle(args.references):
        seqid = header.split()[0]
        total_ref_length += len(seq)
        ref_seq_dict[seqid] = seq.upper()
        ref_positions_count_dict[seqid] = [0 for x in range(len(seq))]
        # aligned_contigs_num, aligned_contigs_length, matches_num
        # mismatches_num, indel_num, overhang_num
        ref_metrics_dict[seqid] = [0, 0, 0, 0, 0, 0]

    # Variables initialization
    previous_query_id = ''
    total_aligned_contigs_num = 0
    total_aligned_contigs_length = 0
    total_matches_num = 0
    total_mismatches_num = 0
    total_indel_num = 0
    total_overhang_num = 0

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

        # Deal with soft clipping and get the query start (0-based)
        query_start = 0
        overhang_num = 0
        if cigar_tab[0][0] == 'S':
            count = cigar_tab[0][1]
            query_start += count
            overhang_num += count
            del cigar_tab[0]

        # Parse CIGAR
        query_end = query_start - 1
        subject_end = subject_start - 1
        indel_num = 0
        matches_num = 0
        mismatches_num = 0
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

        # Store metrics
        if query_id != previous_query_id:
            total_aligned_contigs_num += 1
            total_aligned_contigs_length += query_length
            total_matches_num += matches_num
            total_mismatches_num += mismatches_num
            total_indel_num += indel_num
            total_overhang_num += overhang_num

        # Store metrics for each ref
        ref_metrics_dict[subject_id][0] += 1
        ref_metrics_dict[subject_id][1] += query_length
        ref_metrics_dict[subject_id][2] += matches_num
        ref_metrics_dict[subject_id][3] += mismatches_num
        ref_metrics_dict[subject_id][4] += indel_num
        ref_metrics_dict[subject_id][5] += overhang_num

        #
        ref_positions_count = ref_positions_count_dict[subject_id]
        for i in range(subject_start, subject_end + 1):
            ref_positions_count[i] += 1

        # Store previous subject id and score
        previous_query_id = query_id

    # Final stats
    total_leven_distance = total_mismatches_num + total_indel_num + total_overhang_num
    errors_num_per_kbp = total_leven_distance * 1000.0 / total_aligned_contigs_length

    total_covered_positions_count = 0
    coverage_count_list = [0 for i in range(11)]
    ref_stats_dict = dict()
    for ref_id in ref_seq_dict:
        ref_positions_count = ref_positions_count_dict[ref_id]
        covered_positions_count = 0
        for pos_coverage in ref_positions_count:
            if pos_coverage > 0:
                covered_positions_count += 1
            if pos_coverage >= 10:
                coverage_count_list[10] += 1
            else:
                coverage_count_list[pos_coverage] += 1
        covered_positions_percent = covered_positions_count * 100.0 / len(ref_positions_count)
        mean_coverage = sum(ref_positions_count)/float(len(ref_positions_count))
        median_coverage = sorted(ref_positions_count, key=int)[len(ref_positions_count)//2]
        ref_stats_dict[ref_id] = (covered_positions_percent, mean_coverage, median_coverage)
        total_covered_positions_count += covered_positions_count

    max_coverage = 0
    percent_coverage_list = [0.0 for i in range(11)]
    for i in range(11):
        coverage_count = coverage_count_list[i]
        if coverage_count > 0:
            max_coverage = i
        coverage_percent = coverage_count * 100.0 / total_ref_length
        percent_coverage_list[i] = coverage_percent

    total_ref_coverage = total_covered_positions_count * 100.0 / total_ref_length

    # Output
    #~ sys.stdout.write('Total ref length     = {0}\n\n'.format(total_ref_length))

    sys.stdout.write('1 MATCH PER READ:\n')
    sys.stdout.write('\tTotal aligned contigs num  = {0}\n'.format(total_aligned_contigs_num))
    sys.stdout.write('\tTotal aligned contigs len  = {0}\n\n'.format(total_aligned_contigs_length))

    sys.stdout.write('\tTotal matches num    = {0}\n'.format(total_matches_num))
    sys.stdout.write('\tTotal mismatches num = {0}\n'.format(total_mismatches_num))
    sys.stdout.write('\tTotal indel num      = {0}\n'.format(total_indel_num))
    sys.stdout.write('\tTotal overhang num   = {0}\n\n'.format(total_overhang_num))

    sys.stdout.write('\tTotal leven distance = {0}\n'.format(total_leven_distance))
    sys.stdout.write('\tAssembly error rate  = {0:.2f} errors / kbp\n\n'.format(errors_num_per_kbp))

    sys.stdout.write('ALL BEST MATCHES:\n')
    sys.stdout.write('\tTotal ref coverage   = {0:.2f}%\n'.format(total_ref_coverage))
    sys.stdout.write('\t\tCov')
    for i in range(max_coverage + 1):
        sys.stdout.write('\t{0}'.format(i))
    sys.stdout.write('+\n\t\t%align')
    for i in range(max_coverage + 1):
        sys.stdout.write('\t{0:.2f}%'.format(percent_coverage_list[i]))
    sys.stdout.write('\n\n')

    sys.stdout.write('\tPer-reference Stats:\n')
    sys.stdout.write('\t\tRefID\tCov. Pos. %\tMean Cov.\tMedian Cov.'
                     '\t#Align. Contigs\tAlign. Contigs Lgth\t#Matches'
                     '\t#Mismatches\t#Indel\t#Overhang\n')
    for ref_id, stats_list in sorted(ref_stats_dict.items(), key=lambda x: x[1][1], reverse=True):
        # Store metrics for each ref
        rl = ref_metrics_dict[ref_id]
        sys.stdout.write('\t\t{0}\t{1:.2f}%\t{2:.2f}\t{3}'.format(ref_id, stats_list[0], stats_list[1], stats_list[2]))
        sys.stdout.write('\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(rl[0], rl[1], rl[2], rl[3], rl[4], rl[5]))
