#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


if __name__ == '__main__':
    
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_tab
    parser.add_argument('-i', '--input_tab', 
                        metavar='INBLAST', 
                        type=argparse.FileType('r'), 
                        default='-',
                        help='Input blast tab file. ' 
                             'Assuming sorted by subject, ' 
                             'eval, score, and qcov')
    # -c / --contigs
    parser.add_argument('-c', '--contigs', 
                        metavar='CONTIGS', 
                        type=argparse.FileType('r'), 
                        required=True,
                        help='Assembly contigs fasta file')
    # -o / --output_tab
    parser.add_argument('-o', '--output_tab', 
                        metavar='OUTBLAST', 
                        type=argparse.FileType('w'), 
                        default='-',
                        help='Ouput filtered blast tab file')
    # -p / --threshold_percent
    parser.add_argument('-p', '--threshold_percent',
                        metavar='FLOAT',
                        type=float,
                        default=1.0,
                        help='Score threshold percent')
    
    args = parser.parse_args()
    
    
    # Get contigs length
    contig_length_dict = dict()
    for header, seq in read_fasta_file_handle(args.contigs):
        seqid = header.split()[0]
        contig_length_dict[seqid] = len(seq)
    
    
    # Variables initialization
    previous_query_id = ''
    #~ previous_score = 0
    best_score = -1
    to_write = False
    
    # Reading blast tab file
    for tab in (l.split() for l in args.input_tab if l.split()):
        
        # Parse blast tab
        query_id = tab[0]
        percent_id = float(tab[2])
        query_start = int(tab[6])
        query_end = int(tab[7])
        blast_score = int(tab[11])
        
        query_length = contig_length_dict[query_id]
        
        query_coverage = (query_end - query_start + 1) * 100.0 / query_length
        
        # Compute filtering score
        score = percent_id * query_coverage
        
        #
        if previous_query_id != query_id:
            to_write = True
            best_score = blast_score
        else:
            if blast_score < best_score * args.threshold_percent:
                to_write = False
        
        #~ if query_id == '83':
            #~ print tab, to_write
        
        # Write the tab line if needed
        if to_write:
            args.output_tab.write('{0}\n'.format('\t'.join(tab)))
        
        # Store previous query id and score
        previous_query_id = query_id
        #~ previous_score = blast_score
    

