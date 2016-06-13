#!/usr/bin/env python
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


def analyse_exo_cigar(exo_cigar):
    """
    
    """
    exo_cigar_tab = exo_cigar.split()
    
    matches_or_mismatches = 0
    indel = 0
    
    for i in xrange(0, len(exo_cigar_tab), 2):
        label = exo_cigar_tab[i]
        label_count = int(exo_cigar_tab[i+1])
        
        if label == 'M':
            matches_or_mismatches += label_count
        elif label == 'D' or label == 'I':
            indel += label_count
    
    return matches_or_mismatches, indel


def tab_to_sam(query_id, subject_id, reverse_complement, query_start, 
               query_end, subject_start, subject_end, exo_cigar, 
               contig_seq_dict):
    """
    """
    sam_tab = list()
    
    query_seq = contig_seq_dict[query_id]
    
    # QNAME
    sam_tab.append(query_id)
    
    # FLAG
    flag = 0
    if reverse_complement:
        flag = 16
    sam_tab.append(str(flag))
    
    # RNAME
    sam_tab.append(subject_id)
    
    # POS
    sam_tab.append(str(subject_start + 1))
    
    # MAPQ
    sam_tab.append('255')
    
    # CIGAR
    cigar = str()
    
    if query_start > 0:
        cigar += '{0}S'.format(query_start)
    
    exo_cigar_tab = exo_cigar.split()
    
    for i in xrange(0, len(exo_cigar_tab), 2):
        label = exo_cigar_tab[i]
        label_count = int(exo_cigar_tab[i+1])
        cigar += '{0}{1}'.format(label_count, label)
    
    if query_end > len(query_seq):
        cigar += '{0}S'.format(len(query_seq) - query_end)
    
    sam_tab.append(cigar)
    
    # RNEXT
    sam_tab.append('*')
    
    # PNEXT
    sam_tab.append('0')
    
    # TLEN
    sam_tab.append('0')
    
    # SEQ
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    if reverse_complement:
        sam_tab.append(revcompl(query_seq))
    else:
        sam_tab.append(query_seq)
    
    # QUAL
    sam_tab.append('*')
    
    #
    return '\t'.join(sam_tab)


if __name__ == '__main__':
    
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_tab
    parser.add_argument('-i', '--input_tab', 
                        metavar='INEXOTAB', 
                        type=argparse.FileType('r'), 
                        required=True,
                        help='Input exonerate tab file, sorted by subject. '
                             'Only best alignments are expected.')
    # -c / --contigs
    parser.add_argument('-c', '--contigs', 
                        metavar='CONTIGS', 
                        type=argparse.FileType('r'), 
                        required=True,
                        help='Assembly contigs fasta file')
    # -r / --references
    parser.add_argument('-r', '--references', 
                        metavar='REF', 
                        type=argparse.FileType('r'), 
                        required=True,
                        help='References fasta file (optional). '
                             'If used, the blast file should only '
                             'contain the best alignments')
    # -s / --output_sam
    parser.add_argument('-s', '--output_sam', 
                        metavar='OUTSAM', 
                        type=argparse.FileType('w'), 
                        required=False,
                        help='Optional output sam file')
    
    args = parser.parse_args()
    
    
    # Get contigs length
    contig_length_dict = dict()
    contig_seq_dict = dict()
    for header, seq in read_fasta_file_handle(args.contigs):
        seqid = header.split()[0]
        contig_seq_dict[seqid] = seq
        contig_length_dict[seqid] = len(seq)
    
    
    # Get ref length and initialize positions count
    ref_length_dict = dict()
    ref_positions_count_dict = dict()
    total_ref_length = 0
    for header, seq in read_fasta_file_handle(args.references):
        seqid = header.split()[0]
        total_ref_length += len(seq)
        ref_length_dict[seqid] = len(seq)
        ref_positions_count_dict[seqid] = [0 for x in xrange(len(seq))]
    
    
    # Variables initialization
    previous_query_id = ''
    total_aligned_contigs_num = 0
    total_aligned_contigs_length = 0
    total_matches_num = 0
    total_mismatches_num = 0
    total_indel_num = 0
    total_overhang_num = 0
    
    
    # Reading exonerate tab file
    for tab in (l.split('\t') for l in args.input_tab if l.strip()):
        
        # Parse exonerate tab
        query_id = tab[0]
        subject_id = tab[1]
        
        query_start = int(tab[2])
        query_end = int(tab[3])
        subject_start = int(tab[4])
        subject_end = int(tab[5])
        
        exo_cigar = tab[6]
        mismatches_num = int(tab[7])
        
        # Compute additional metrics
        reverse_complement = False
        if query_end < query_start:
            reverse_complement = True
            # Reverse query start and end if reversed complemented
            tmp = query_start
            query_start = query_end
            query_end = tmp
        
        matches_or_mismatches, indel_num = analyse_exo_cigar(exo_cigar)
        matches_num = matches_or_mismatches - mismatches_num
        
        alignment_on_query_length = query_end - query_start
        
        query_length = contig_length_dict[query_id]
        
        overhang_num = 0
        
        if query_start > 0:
            overhang_num += query_start
        if query_end < query_length:
            overhang_num += query_length - query_end
        
        # Store metrics
        if query_id != previous_query_id:
            total_aligned_contigs_num += 1
            total_aligned_contigs_length += query_length
            total_matches_num += matches_num
            total_mismatches_num += mismatches_num
            total_indel_num += indel_num
            total_overhang_num += overhang_num
        
        #
        ref_positions_count = ref_positions_count_dict[subject_id]
        for i in xrange(subject_start-1, subject_end):
            ref_positions_count[i] += 1
        
        # Write optional sam output
        if args.output_sam:
            sam_line = tab_to_sam(query_id, subject_id, reverse_complement,
                                  query_start, query_end,
                                  subject_start, subject_end,
                                  exo_cigar, contig_seq_dict)
            args.output_sam.write('{0}\n'.format(sam_line))
        
        # Store previous subject id and score
        previous_query_id = query_id
    
    # Final stats
    total_leven_distance = total_mismatches_num + total_indel_num + total_overhang_num
    errors_num_per_kbp = total_leven_distance * 1000.0 / total_aligned_contigs_length
    
    total_covered_positions_count = 0
    coverage_count_list = [0 for i in xrange(11)]
    for ref_id, ref_length in ref_length_dict.items():
        covered_positions_count = 0
        ref_positions_count = ref_positions_count_dict[ref_id]
        for pos_coverage in ref_positions_count:
            if pos_coverage > 0:
                covered_positions_count += 1
            if pos_coverage >= 10:
                coverage_count_list[10] += 1
            else:
                coverage_count_list[pos_coverage] += 1
        total_covered_positions_count += covered_positions_count
    
    max_coverage = 0
    percent_coverage_list = [0.0 for i in xrange(11)]
    for i in xrange(11):
        coverage_count = coverage_count_list[i]
        if coverage_count > 0:
            max_coverage = i
        coverage_percent = coverage_count * 100.0 / total_ref_length
        percent_coverage_list[i] = coverage_percent
    
    total_ref_coverage = total_covered_positions_count * 100.0 / total_ref_length
    
    # Output
    sys.stdout.write('Total aligned contigs num  = {0}\n'.format(total_aligned_contigs_num))
    sys.stdout.write('Total aligned contigs len  = {0}\n\n'.format(total_aligned_contigs_length))
    
    sys.stdout.write('Total ref length     = {0}\n\n'.format(total_ref_length))
    
    sys.stdout.write('Total matches num    = {0}\n'.format(total_matches_num))
    sys.stdout.write('Total mismatches num = {0}\n'.format(total_mismatches_num))
    sys.stdout.write('Total indel num      = {0}\n'.format(total_indel_num))
    sys.stdout.write('Total overhang num   = {0}\n\n'.format(total_overhang_num))
    
    sys.stdout.write('Total leven distance = {0}\n'.format(total_leven_distance))
    sys.stdout.write('Assembly error rate  = {0:.2f} errors / kbp\n\n'.format(errors_num_per_kbp))
    
    sys.stdout.write('Total ref coverage   = {0:.2f}%\n'.format(total_ref_coverage))
    sys.stdout.write('\tCov')
    for i in xrange(max_coverage + 1):
        sys.stdout.write('\t{0}'.format(i))
    sys.stdout.write('+\n\t%align')
    for i in xrange(max_coverage + 1):
        sys.stdout.write('\t{0:.2f}%'.format(percent_coverage_list[i]))
    sys.stdout.write('\n')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

