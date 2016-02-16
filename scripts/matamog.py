#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re

def parse_arguments():
    #
    parser = argparse.ArgumentParser(description='')
    #
    group_main = parser.add_argument_group('Main Parameters')
    group_main.add_argument('-i', '--input_fastx', metavar='FASTX', 
                            type=argparse.FileType('r', 0), default='-',
                            help='Input reads file (fasta or fastq format)')
    group_main.add_argument('-d', '--ref_db', metavar='FASTA',
                            type=argparse.FileType('r', 0),
                            help='Reference database (fasta format, with Silva taxo headers)')
    group_main.add_argument('-o', '--output_contigs', metavar='OUTPUT', 
                            type=str, default='DEFAULTNAME',
                            help='Output contigs file (fasta format)')
    group_main.add_argument('-s', '--steps', metavar='INT',
                            type=int, nargs='+', default=[1,2,3,4,5],
                            help='Steps to execute')
    #
    group_perf = parser.add_argument_group('Performance')
    group_perf.add_argument('--cpu', metavar='INT',
                            type=int, default=3,
                            help='Max number of CPU to use')
    #
    group_db = parser.add_argument_group('Ref DB pre-processing (Step 0)')
    group_db.add_argument('--clustering_id_threshold', metavar='REAL',
                          type=float, default=0.95,
                          help='Identity threshold for clustering')
    #
    group_mapping = parser.add_argument_group('Mapping (Step 1)')
    group_mapping.add_argument('--best', metavar='INT',
                               type=int, default=10,
                               help='Get up to --best good alignments per read')
    #
    group_filt = parser.add_argument_group('Alignment Filtering (Step 2)')
    group_filt.add_argument('--score_threshold', metavar='REAL',
                            type=float, default=0.9,
                            help='Score threshold (real between 0 and 1)')
    group_filt.add_argument('--geometric_mode', action='store_true',
                            help='Use geometric mode filtering')
    #
    group_ovg = parser.add_argument_group('Overlap Graph Building (Step 3)')
    group_ovg.add_argument('--min_identity', metavar='REAL',
                           type=float, default=1.0,
                           help='Minimum identity of an overlap between 2 reads')
    group_ovg.add_argument('--min_overlap_length', metavar='INT',
                           type=int, default=50,
                           help='Minimum length of an overlap')
    group_ovg.add_argument('--multi', action='store_true',
                           help='Use multi-ref mode')
    #
    group_gcomp = parser.add_argument_group('Graph Compaction & Contig Identification (Step 4)')
    group_gcomp.add_argument('--min_read_node', metavar='INT',
                             type=int, default=1,
                             help='Minimum number of read to keep a node')
    group_gcomp.add_argument('--min_overlap_edge', metavar='INT',
                             type=int, default=1,
                             help='Minimum number of overlap to keep an edge')
    #
    group_lca = parser.add_argument_group('LCA Labelling (Step 5)')
    group_lca.add_argument('--quorum', metavar='FLOAT',
                           type=float, default=0.51,
                           help='Quorum for LCA computing')
    #
    group_ass = parser.add_argument_group('Contig Assembly (Step 6)')
    
    args = parser.parse_args()
    #
    return args
    

if __name__ == '__main__':
    
    # Arguments parsing
    args = parse_arguments()
    
    steps_set = frozenset(args.steps)
    print steps_set
    
    # STEP 1: Mappping
    if 1 in steps_set:
        sys.stdout.write('## Mapping step (1):\n')
        
        
        
    
    
    exit(0)

















