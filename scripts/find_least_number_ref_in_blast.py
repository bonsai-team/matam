#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict


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


if __name__ == '__main__':
    
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_tab', metavar='BLAST', 
                        type=argparse.FileType('r'), default='-',
                        help='Input blast tab file')
    parser.add_argument('-o', '--output_tab', metavar='BLAST', 
                        type=argparse.FileType('w'), default='-',
                        help='Ouput filtered blast tab file')
    args = parser.parse_args()
    
    #
    subject_alignments_dict = defaultdict(list)
    subject_query_num_dict = defaultdict(int)
    subject_pident_length_dict = defaultdict(list)
    subject_length_dict = dict()
    
    #
    for tab in [l.split() for l in args.input_tab]:
        
        # Blast tab parsing
        #~ qseqid = tab[0]
        sseqid = tab[1]
        pident = float(tab[2])
        length = int(tab[3])
        #~ mismatch = int(tab[4])
        #~ gapopen = int(tab[5])
        #~ qstart = int(tab[6])
        #~ qend = int(tab[7])
        sstart = int(tab[8])
        send = int(tab[9])
        evalue = float(tab[10])
        #~ bitscore = float(tab[11])
        qlen = int(tab[12])
        slen = int(tab[13])
        revcomp = (sstart > send)
        
        #
        if sseqid not in subject_length_dict:
            subject_length_dict[sseqid] = slen
        
        #
        #~ if length >= (0.8 * qlen) and pident >= 80.0:
        if revcomp:
            subject_alignments_dict[sseqid].append((send, sstart))
        else:
            subject_alignments_dict[sseqid].append((sstart, send))
        
        subject_query_num_dict[sseqid] += 1
        
        subject_pident_length_dict[sseqid].append((pident, length))
    
    
    # Find query_num min and max
    min_query_num = min(subject_query_num_dict.itervalues())
    max_query_num = max(subject_query_num_dict.itervalues())
    
    print "min max query_num: ", min_query_num, max_query_num
    
    # Sort the alignments by starting position
    for sseqid in subject_alignments_dict:
        subject_alignments_dict[sseqid].sort(key=lambda x: x[0])
    
    # Compute the subjects coverage
    subject_coverage_length_dict = dict()
    
    for sseqid, alignment_list in subject_alignments_dict.items():
        slen = subject_length_dict[sseqid]
        
        interval_start_pos = 0
        interval_end_pos = 0
        coverage_length = 0
        
        for start_pos, end_pos in alignment_list:
            #~ print start_pos, end_pos, interval_start_pos, interval_end_pos, coverage_length
            if start_pos > interval_end_pos:
                if interval_end_pos:
                    #Add the interval length to the coverage length
                    coverage_length += (interval_end_pos - interval_start_pos) + 1
                interval_start_pos = start_pos
            if end_pos > interval_end_pos:
                interval_end_pos = end_pos
        coverage_length += (interval_end_pos - interval_start_pos) + 1
        coverage_length_percent = coverage_length * 100.0 / slen
        
        #~ print sseqid, slen, coverage_length, coverage_length_percent, alignment_list
        
        subject_coverage_length_dict[sseqid] = coverage_length_percent
    
    min_coverage_length = min(subject_coverage_length_dict.itervalues())
    max_coverage_length = max(subject_coverage_length_dict.itervalues())
    
    print "min max coverage_length :", min_coverage_length, max_coverage_length
    
    # Compute the average pident for each subject
    subject_average_pident_dict = dict()
    
    for sseqid, pident_length_list in subject_pident_length_dict.items():
        sum_length = 0
        sum_weighted_pident = 0.0
        for pident, length in pident_length_list:
            sum_length += length
            sum_weighted_pident += length * pident
        average_pident = float(sum_weighted_pident) / sum_length
        subject_average_pident_dict[sseqid] = average_pident
    
        #~ print sseqid, sum_length, sum_weighted_pident, average_pident, pident_length_list
    
    min_average_pident = min(subject_average_pident_dict.itervalues())
    max_average_pident = max(subject_average_pident_dict.itervalues())
    
    print "min max average_pident :", min_average_pident, max_average_pident
    
    # Compute the sorting factor
    
    subject_sorting_factor_dict = dict()
    
    for sseqid, query_num in subject_query_num_dict.items():
        coverage_length_percent = subject_coverage_length_dict[sseqid]
        average_pident = subject_average_pident_dict[sseqid]
        
        #~ print sseqid, query_num, coverage_length_percent, average_pident
        
        # No sorting factor: 101, 99.23%, 99.79% 
        
        #~ sorting_factor = query_num # 16, 93.36%, 88.52%
        #~ sorting_factor = coverage_length_percent # 85, 99.05%, 99.37%
        #~ sorting_factor = average_pident # 48, 86.37%, 93.13%
        #~ sorting_factor = query_num * average_pident # 14, 93.08%, 89.70%
        #~ sorting_factor = query_num * coverage_length_percent # 16, 93.31%, 88.52%
        #~ sorting_factor = coverage_length_percent * average_pident # 15, 90.32%, 92.17%
        #~ sorting_factor = query_num * coverage_length_percent * average_pident # 14, 93.12%, 89.70%
        
        norm_query_num = float(query_num - min_query_num)/max_query_num
        norm_coverage_length_percent = float(coverage_length_percent - min_coverage_length)/max_coverage_length
        norm_average_pident = float(average_pident - min_average_pident)/max_average_pident
        
        norm_sorting_factor = norm_query_num * norm_coverage_length_percent * norm_average_pident # 9, 90.46%, 93.03%
        #~ norm_sorting_factor = norm_query_num * norm_coverage_length_percent**1.2 * norm_average_pident**0.8 # 11, 91.03%, 93.09%
        #~ norm_sorting_factor = norm_query_num * norm_average_pident # 9, 90.46%, 92.98%
        #~ norm_sorting_factor = norm_coverage_length_percent * norm_average_pident # 14, 90.58%, 93.24%
        #~ norm_sorting_factor = norm_coverage_length_percent * norm_query_num # 16, 93.31%, 88.52%
        
        #~ subject_sorting_factor_dict[sseqid] = sorting_factor
        subject_sorting_factor_dict[sseqid] = norm_sorting_factor
    
    # Reset input file reading
    args.input_tab.seek(0)
    
    # Variable initialisation
    subject_list = list()
    pident_length_list = list()
    aligned_length = 0
    total_query_length = 0
    
    # Sort and find the best alignment
    for tab_list in read_tab_file_handle_sorted(args.input_tab):
        #~ tab_list.sort(key = lambda x: subject_coverage_length_dict[x[1]], reverse = True)
        #~ tab_list.sort(key = lambda x: subject_sorting_factor_dict[x[1]], reverse = True)
        tab_list.sort(key = lambda x: (-subject_sorting_factor_dict[x[1]], float(x[10])))
        
        best_tab = tab_list[0]
        
        # Blast tab parsing
        #~ qseqid = best_tab[0]
        sseqid = best_tab[1]
        pident = float(best_tab[2])
        length = int(best_tab[3])
        #~ mismatch = int(best_tab[4])
        #~ gapopen = int(best_tab[5])
        qstart = int(best_tab[6])
        qend = int(best_tab[7])
        sstart = int(best_tab[8])
        send = int(best_tab[9])
        evalue = float(best_tab[10])
        #~ bitscore = float(best_tab[11])
        qlen = int(best_tab[12])
        slen = int(best_tab[13])
        revcomp = (sstart > send)
        
        subject_list.append(sseqid)
        pident_length_list.append((pident, length))
        
        total_query_length += qlen
        aligned_length += (qend - qstart)
        
        #~ print best_tab
        #~ print sseqid, (pident, length), qlen, total_query_length, aligned_length
        
        args.output_tab.write('{0}\n'.format('\t'.join(tab_list[0])))
    
    # Compute final stats
    uniq_subject_set = set(subject_list)
    aligned_length_percent = float(aligned_length) * 100.0 / total_query_length
    
    # Compute global average pident
    sum_length = 0
    sum_weighted_pident = 0.0
    for pident, length in pident_length_list:
        sum_length += length
        sum_weighted_pident += length * pident
    global_average_pident = float(sum_weighted_pident) / sum_length
    
    # Print final stats
    print "Final subject number : ", len(uniq_subject_set)
    print "Contig length aligned : {0:.2f}%".format(aligned_length_percent)
    print "Average pident : {0:.2f}%".format(global_average_pident)
    #~ print '{0}, {1:.2f}%, {2:.2f}%'.format(len(uniq_subject_set),
                                           #~ aligned_length_percent,
                                           #~ global_average_pident)
    
    
    


