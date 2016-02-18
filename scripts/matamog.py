#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
import subprocess


matamog_dir = os.path.realpath(sys.argv[0])[:-19]
sumaclust_bin = matamog_dir + '/bin/sumaclust'
extract_taxo_bin = matamog_dir + '/bin/extract_taxo_from_fasta.py'
clean_name_bin = matamog_dir + '/bin/fasta_clean_name.py'
name_filter_bin = matamog_dir + '/bin/fasta_name_filter.py'
indexdb_bin = matamog_dir + '/bin/indexdb_rna'
sortmerna_bin = matamog_dir + '/bin/sortmerna'
filter_score_bin = matamog_dir + '/bin/filter_score_multialign.py'
ovgraphbuild_bin = matamog_dir + '/bin/ovgraphbuild'

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
                            type=str, required=True,
                            help='Input reads file (fasta or fastq format)')
    group_main.add_argument('-d', '--ref_db', metavar='FASTA',
                            type=str, required=True,
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
    group_mapping.add_argument('--min_lis', metavar='INT',
                               type=int, default=10,
                               help=argparse.SUPPRESS)
    group_mapping.add_argument('--evalue', metavar='REAL',
                               type=float, default=1e-10,
                               help='Max e-value to keep an alignment for (default: 10e-10)')
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
                           help='Quorum for LCA computing. Has to be between 0.51 and 1')
    #
    group_ass = parser.add_argument_group('Contig Assembly (Step 7)')
    #
    args = parser.parse_args()
    
    # Arguments checking
    
    # check if steps are in proper range
    steps_set = set(args.steps)
    if len(steps_set - set(range(1, 8))) > 0:
        raise Exception("steps not in range [1,7]")
    
    if args.clustering_id_threshold < 0 or args.clustering_id_threshold > 1:
        raise Exception("clustering id threshold not in range [0,1]")
    
    if args.score_threshold < 0 or args.score_threshold > 1:
        raise Exception("score threshold not in range [0,1]")
    
    if args.min_identity < 0 or args.min_identity > 1:
        raise Exception("min identity not in range [0,1]")
    
    if args.quorum < 0 or args.quorum > 1:
        raise Exception("quorum not in range [0.51,1]")
    
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
    sys.stdout.write('PARAM: Input Fastx = {0}\n'.format(args.input_fastx))
    sys.stdout.write('PARAM: Ref DB      = {0}\n'.format(args.ref_db))
    sys.stdout.write('PARAM: Steps       = {0}\n'.format(args.steps))
    sys.stdout.write('\n')
    
    return 0


if __name__ == '__main__':
    
    # Arguments parsing
    args = parse_arguments()
    
    #
    print_intro(args)
    
    #
    ref_db_filename = args.ref_db.split('/')[-1]
    ref_db_basename = '.'.join(ref_db_filename.split('.')[:-1])
    
    input_fastx_filename = args.input_fastx.split('/')[-1]
    input_fastx_basename = '.'.join(input_fastx_filename.split('.')[:-1])
    
    #
    steps_set = frozenset(args.steps)
    
    ###############################
    # STEP 1: Ref DB pre-processing
    cluster_id_int = int(args.clustering_id_threshold * 100)
    clustered_ref_db_basename = 'NR{0}'.format(cluster_id_int)
    sortmerna_index_directory = 'sortmerna_index'
    sortmerna_index_basename = sortmerna_index_directory + '/' + clustered_ref_db_basename
    
    if 1 in steps_set:
        sys.stdout.write('## Ref DB pre-processing step (1):\n\n')
        
        # Extract taxo from ref DB and sort by ref id
        
        ref_db_taxo_filename = ref_db_basename + '.taxo.tab'
        
        # Clean fasta headers
        
        # Convert Us in Ts
        
        # Trim Ns from both sides
        
        # Filter out small sequences
        
        # Option: Either filter out seq with Ns or replace Ns with random nucl
        
        # Sort sequences by decreasing length
        
        cleaned_ref_db_filename = ref_db_basename + '.cleaned.fasta'
        
        ## FOR EACH kingdom IN ['archaea', 'bacteria', 'eukaryota']
        ## Clutering with sumaclust
        
        sumaclust_filename = ref_db_basename + '.sumaclust_' 
        sumaclust_filename += '{0}'.format(cluster_id_int)
        sumaclust_filename += '.fasta'
        
        sumaclust_cmd_line = sumaclust_bin + ' -t ' + '{0:.2f}'.format(args.clustering_id_threshold)
        sumaclust_cmd_line += ' -p ' + str(args.cpu) + ' -F ' + sumaclust_filename
        sumaclust_cmd_line += ' ' + cleaned_ref_db_filename
        
        sys.stdout.write('CMD: {0}\n\n'.format(sumaclust_cmd_line))
        
        ## Extracting centroids
        
        filter_cmd_line = name_filter_bin + ' -s "cluster_center=True" -i '
        filter_cmd_line += sumaclust_filename + ' -o tmp.fasta'
        
        sys.stdout.write('CMD: {0}\n'.format(filter_cmd_line))
        
        ## Concatenate kingdom centroids
        
        sumaclust_kingdom_filename = ref_db_basename + '.sumaclust_'
        sumaclust_kingdom_filename += '{0}'.format(cluster_id_int)
        sumaclust_kingdom_filename += '_by_kingdom.fasta'
        
        # Clean fasta headers
        clean_name_cmd_line = clean_name_bin + ' -i ' + sumaclust_kingdom_filename
        clean_name_cmd_line += ' -o ' + clustered_ref_db_basename + '.fasta'
        
        sys.stdout.write('CMD: {0}\n'.format(clean_name_cmd_line))
        subprocess.call(clean_name_cmd_line, shell=True)
        
        # SortMeRNA Ref DB indexing
        try:
            if not os.path.exists(sortmerna_index_directory):
                os.makedirs(sortmerna_index_directory)
        except OSError:
            sys.stderr.write("\nERROR: {0} cannot be created\n\n".format(sortmerna_index_directory))
            raise
        
        indexdb_cmd_line = indexdb_bin + ' -v --ref ' + clustered_ref_db_basename
        indexdb_cmd_line += '.fasta,' + sortmerna_index_basename
        
        sys.stdout.write('CMD: {0}\n'.format(indexdb_cmd_line))
        subprocess.call(indexdb_cmd_line, shell=True)
        
    ######################################
    # STEP 2: Reads mapping against Ref DB
    sortme_output_basename = input_fastx_basename + '.S2.1_vs_' + clustered_ref_db_basename
    sortme_output_basename += '_b' + str(args.best) + '_m' + str(args.min_lis)
    
    if 2 in steps_set:
        sys.stdout.write('## Mapping step (2):\n\n')
        
        # Set SortMeRNA command line
        sortmerna_cmd_line = sortmerna_bin + ' --ref ' + clustered_ref_db_basename
        sortmerna_cmd_line += '.fasta,' + sortmerna_index_basename + ' --reads '
        sortmerna_cmd_line += args.input_fastx + ' --aligned ' + sortme_output_basename
        sortmerna_cmd_line += ' --fastx --sam --blast "1 cigar qcov" --log --best '
        sortmerna_cmd_line += str(args.best) + ' --min_lis ' + str(args.min_lis) 
        sortmerna_cmd_line += ' -e {0:.2e}'.format(args.evalue)
        sortmerna_cmd_line += ' -a ' + str(args.cpu) + ' -v'
        
        # Run SortMeRNA
        sys.stdout.write('CMD: {0}\n'.format(sortmerna_cmd_line))
        subprocess.call(sortmerna_cmd_line, shell=True)
    
    #############################
    # STEP 3: Alignment Filtering
    score_threshold_int = int(args.score_threshold * 100)
    sam_filt_basename = sortme_output_basename + '.scr_filt_'
    if args.geometric_mode:
        sam_filt_basename += 'geo_'
    sam_filt_basename += str(score_threshold_int) + 'pct'
    read_ref_taxo_basename = sam_filt_basename + '.read_ref_taxo'
    
    if 3 in steps_set:
        sys.stdout.write('## Alignment filtering step (3):\n\n')
        
        # Filtering scores command line
        filter_score_cmd_line = 'cat ' + sortme_output_basename + '.sam'
        filter_score_cmd_line += ' | grep -v "^@" | sort -k 1,1V -k 12,12Vr'
        filter_score_cmd_line += ' | ' + filter_score_bin + ' -t ' + str(args.score_threshold)
        if args.geometric_mode:
            filter_score_cmd_line += ' --geometric'
        filter_score_cmd_line += ' > ' + sam_filt_basename + '.sam'
        
        # Run scores filtering
        sys.stdout.write('CMD: {0}\n'.format(filter_score_cmd_line))
        subprocess.call(filter_score_cmd_line, shell=True)
        
        # Generate read ref taxo file
        read_taxo_cmd_line = 'cat ' + sam_filt_basename + '.sam | cut -f1,3 | sort -k2,2'
        read_taxo_cmd_line += ' | join -12 -21 - ' + ref_db_basename + '.taxo.tab'
        read_taxo_cmd_line += ' | sort -k2,2 | awk "{print \$2,\$1,\$3}" | sed "s/ /\\t/g" > '
        read_taxo_cmd_line += read_ref_taxo_basename + '.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(read_taxo_cmd_line))
        subprocess.call(read_taxo_cmd_line, shell=True)
    
    ################################
    # STEP 4: Overlap Graph Building
    min_identity_int = int(args.min_identity * 100)
    ovgraphbuild_basename = sam_filt_basename + '.ovgb_i' + str(min_identity_int)
    ovgraphbuild_basename += '_o' + str(args.min_overlap_length) + '_'
    if args.multi:
        ovgraphbuild_basename += 'multi'
    else:
        ovgraphbuild_basename += 'single'
    
    if 4 in steps_set:
        sys.stdout.write('## Overlap Graph building step (4):\n\n')
        
        # Ovgraphbuild command line
        ovgraphbuild_cmd_line = ovgraphbuild_bin + ' -v --debug'
        ovgraphbuild_cmd_line += ' -i ' + str(args.min_identity)
        ovgraphbuild_cmd_line += ' -m ' + str(args.min_overlap_length)
        if args.multi:
            ovgraphbuild_cmd_line += ' --multi_ref'
        ovgraphbuild_cmd_line += ' --asqg --csv --output_basename '
        ovgraphbuild_cmd_line += ovgraphbuild_basename
        ovgraphbuild_cmd_line += ' -s ' + sam_filt_basename + '.sam'
        
        # Run ovgraphbuild
        sys.stdout.write('CMD: {0}\n'.format(ovgraphbuild_cmd_line))
        subprocess.call(ovgraphbuild_cmd_line, shell=True)
    
    
    exit(0)

















