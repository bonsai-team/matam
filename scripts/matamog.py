#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re

matamog_dir = os.path.realpath(sys.argv[0])[:-18]
sumaclust_bin = matamog_dir + '/bin/sumaclust'
extract_taxo_bin = matamog_dir + '/bin/extract_taxo_from_fasta.py'
clean_name_bin = matamog_dir + '/bin/fasta_clean_name.py'


class DefaultHelpParser(argparse.ArgumentParser):
    """
    This is a slightly modified argparse parser to display the full help
    on parser error instead of only usage
    """
    def error(self, message):
        sys.stderr.write('\nerror: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def parse_arguments():
    """
    Parse the command line, and check if arguments are correct
    """
    #
    parser = DefaultHelpParser(description='Matamog',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))
    #
    group_main = parser.add_argument_group('Main Parameters')
    group_main.add_argument('-i', '--input_fastx', metavar='FASTX', 
                            type=argparse.FileType('r', 0), default='-',
                            help='Input reads file (fasta or fastq format)')
    group_main.add_argument('-d', '--ref_db', metavar='FASTA',
                            type=argparse.FileType('r', 0), required=True,
                            help='Reference database (fasta format, with Silva taxo headers)')
    group_main.add_argument('-o', '--output_contigs', metavar='OUTPUT', 
                            type=str, default='DEFAULTNAME',
                            help='Output contigs file (fasta format)')
    group_main.add_argument('-s', '--steps', metavar='INT',
                            type=int, nargs='+', default=[1,2,3,4,5,6,7],
                            help='Steps to execute')
    #
    group_perf = parser.add_argument_group('Performance')
    group_perf.add_argument('--cpu', metavar='INT',
                            type=int, default=3,
                            help='Max number of CPU to use')
    group_perf.add_argument('--max_memory', metavar='INT',
                            type=int, default=-1,
                            help='Maximum memory to use (in MBi). Default is auto')
    #
    group_db = parser.add_argument_group('Ref DB pre-processing (Step 1)')
    group_db.add_argument('--clustering_id_threshold', metavar='REAL',
                          type=float, default=0.95,
                          help='Identity threshold for clustering')
    #
    group_mapping = parser.add_argument_group('Mapping (Step 2)')
    group_mapping.add_argument('--best', metavar='INT',
                               type=int, default=10,
                               help='Get up to --best good alignments per read')
    #
    group_filt = parser.add_argument_group('Alignment Filtering (Step 3)')
    group_filt.add_argument('--score_threshold', metavar='REAL',
                            type=float, default=0.9,
                            help='Score threshold (real between 0 and 1)')
    group_filt.add_argument('--geometric_mode', action='store_true',
                            help='Use geometric mode filtering')
    #
    group_ovg = parser.add_argument_group('Overlap Graph Building (Step 4)')
    group_ovg.add_argument('--min_identity', metavar='REAL',
                           type=float, default=1.0,
                           help='Minimum identity of an overlap between 2 reads')
    group_ovg.add_argument('--min_overlap_length', metavar='INT',
                           type=int, default=50,
                           help='Minimum length of an overlap')
    group_ovg.add_argument('--multi', action='store_true',
                           help='Use multi-ref mode')
    #
    group_gcomp = parser.add_argument_group('Graph Compaction & Contig Identification (Step 5)')
    group_gcomp.add_argument('--min_read_node', metavar='INT',
                             type=int, default=1,
                             help='Minimum number of read to keep a node')
    group_gcomp.add_argument('--min_overlap_edge', metavar='INT',
                             type=int, default=1,
                             help='Minimum number of overlap to keep an edge')
    #
    group_lca = parser.add_argument_group('LCA Labelling (Step 6)')
    group_lca.add_argument('--quorum', metavar='FLOAT',
                           type=float, default=0.51,
                           help='Quorum for LCA computing')
    #
    group_ass = parser.add_argument_group('Contig Assembly (Step 7)')
    #
    args = parser.parse_args()
    
    # Arguments checking
    
    # check if steps are in proper range
    steps_set = set(args.steps)
    if len(steps_set - set(range(1, 8))) > 0:
        raise Exception("steps not in range [1,7]")
    
    
    
    #
    return args


def print_intro(args):
    """
    Print the introduction
    """
    
    sys.stdout.write("""
#################################
           MATAMOG
#################################
    \n""")
    
    #~ sys.stdout.write('PARAM: Executable  = {0}\n'.format(os.path.abspath(os.path.dirname(os.path.realpath(sys.argv[0])))))
    sys.stdout.write('PARAM: Matamog dir = {0}\n'.format(matamog_dir))
    sys.stdout.write('PARAM: Input Fastx = {0}\n'.format(args.input_fastx.name))
    sys.stdout.write('PARAM: Ref DB      = {0}\n'.format(args.ref_db.name))
    sys.stdout.write('PARAM: Steps       = {0}\n'.format(args.steps))
    sys.stdout.write('\n')
    
    return 0


if __name__ == '__main__':
    
    # Arguments parsing
    args = parse_arguments()
    
    #
    print_intro(args)
    
    #
    steps_set = frozenset(args.steps)
    
    # STEP 1: Ref DB pre-processing
    if 1 in steps_set:
        sys.stdout.write('## Ref DB pre-processing step (1):\n')
        
        #~ sumaclust_command_line = sumaclust_bin + ' '
        
        
    
    
    exit(0)

















