#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
import subprocess

matam_bin = os.path.realpath(sys.argv[0])
matam_dir = matam_bin[:-17]
sumaclust_bin = matam_dir + '/bin/sumaclust'
extract_taxo_bin = matam_dir + '/bin/extract_taxo_from_fasta.py'
clean_name_bin = matam_dir + '/bin/fasta_clean_name.py'
fasta_name_filter_bin = matam_dir + '/bin/fasta_name_filter.py'
fasta_length_filter_bin = matam_dir + '/bin/fasta_length_filter.py'
fastq_name_filter_bin = matam_dir + '/bin/fastq_name_filter.py'
indexdb_bin = matam_dir + '/bin/indexdb_rna'
sortmerna_bin = matam_dir + '/bin/sortmerna'
sga_bin = matam_dir + '/bin/sga'
filter_score_bin = matam_dir + '/bin/filter_score_multialign.py'
ovgraphbuild_bin = matam_dir + '/bin/ovgraphbuild'
componentsearch_jar = matam_dir + '/bin/ComponentSearch.jar'
compute_lca_bin = matam_dir + '/bin/compute_lca_from_tab.py'
compute_stats_lca_bin = matam_dir + '/bin/compute_stats_from_lca.py'
compute_compressed_graph_stats_bin = matam_dir + '/bin/compute_compressed_graph_stats.py'
replace_Ns_bin = matam_dir + '/bin/replace_Ns_by_rand_nu.py'
sort_fasta_bin = matam_dir + '/bin/sort_fasta_by_length.py'
sga_assemble_bin = matam_dir + '/bin/sga_assemble.py'
find_least_bin = matam_dir + '/bin/find_least_number_ref_in_blast.py'


def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in (l.strip() for l in fasta_file_handle if l.strip()):
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


def read_tab_file_handle_sorted(tab_file_handle, factor_index):
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


def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in xrange(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


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
    parser = DefaultHelpParser(description='matam',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))
    
    # Main parameters
    group_main = parser.add_argument_group('Main Parameters')
    # -i / --input_fastx
    group_main.add_argument('-i', '--input_fastx', 
                            action='store',
                            metavar='FASTX', 
                            type=str, 
                            required=True,
                            help='Input reads file (fasta or fastq format)')
    # -d / --ref_db
    group_main.add_argument('-d', '--ref_db', 
                            action='store',
                            metavar='FASTA',
                            type=str, 
                            required=True,
                            help="Reference database "
                                 "(fasta format, with Silva taxo headers)")
    # -o / --output_contigs
    group_main.add_argument('-o', '--output_contigs', 
                            action='store',
                            metavar='OUTPUT', 
                            type=str, 
                            default='DEFAULTNAME',
                            help='Output contigs file (fasta format)')
    # -s / --steps
    group_main.add_argument('-s', '--steps', 
                            action='store',
                            metavar='INT',
                            type=int, 
                            nargs='+', 
                            default=[2,3,4,5,6,7,8],
                            help="Steps to execute. "
                                 "Default is %(default)s")
    # --db_wkdir
    group_main.add_argument('--db_wkdir',
                            action='store',
                            metavar='DBDIR',
                            type=str,
                            default='.',
                            help='Working directory for the database.'
                                 ' Default is cwd')
    # --wkdir
    group_main.add_argument('--wkdir',
                            action='store',
                            metavar='WKDIR',
                            type=str,
                            default='.',
                            help='Working directory starting from the alignment step.'
                                 ' Default is cwd')
    # -v / --verbose
    group_main.add_argument('-v', '--verbose',
                            action='count',
                            default=0,
                            help='Set verbose level (from 0 to 3) using -vvv or -v -v.')
    
    # Performance parameters
    group_perf = parser.add_argument_group('Performance')
    # --cpu
    group_perf.add_argument('--cpu', 
                            action='store',
                            metavar='CPU',
                            type=int, 
                            default=3,
                            help='Max number of CPU to use. '
                                 'Default is %(default)s cpu')
    # --max_memory
    group_perf.add_argument('--max_memory', 
                            action='store',
                            metavar='MAXMEM',
                            type=int, 
                            default=4000,
                            help="Maximum memory to use (in MBi). "
                                 "Default is %(default)sMBi")
    
    # Step 0: Ref DB pre-processing
    group_db = parser.add_argument_group('Ref DB pre-processing (Step 0)')
    # --remove_Ns
    group_db.add_argument('--remove_Ns', 
                          action='store_true',
                          help="Remove sequences with Ns. "
                               "Default is replacing Ns with random nucleotides")
    
    # Step 1: Clustering Ref DB
    group_clust = parser.add_argument_group('Clustering Ref DB (Step 1)')
    # --clustering_id_threshold
    group_clust.add_argument('--clustering_id_threshold',
                             action='store',
                             metavar='REAL',
                             type=float, 
                             default=0.95,
                             help='Identity threshold for clustering. '
                                  'Default is %(default)s')
    # --kingdoms
    group_clust.add_argument('--kingdoms', 
                             action='store',
                             metavar='STR',
                             type=str, 
                             nargs='+', 
                             default=['archaea','bacteria','eukaryota'],
                             help='Kingdoms to clusterize the DB for. '
                                  'Default is %(default)s')
    # --index_only
    group_clust.add_argument('--index_only', 
                             action='store_true',
                             help=argparse.SUPPRESS)
    
    # Step 2: Mapping
    group_mapping = parser.add_argument_group('Mapping (Step 2)')
    # --best
    group_mapping.add_argument('--best', 
                               action='store',
                               metavar='INT',
                               type=int, 
                               default=10,
                               help='Get up to --best good alignments per read. '
                                    'Default is %(default)s')
    # --min_lis
    group_mapping.add_argument('--min_lis', 
                               action='store',
                               metavar='INT',
                               type=int, 
                               default=10,
                               help=argparse.SUPPRESS)
    # --evalue
    group_mapping.add_argument('--evalue', 
                               action='store',
                               metavar='REAL',
                               type=float, 
                               default=1e-5,
                               help='Max e-value to keep an alignment for. '
                                    'Default is %(default)s')
    
    # Step 3: Alignment Filtering
    group_filt = parser.add_argument_group('Alignment Filtering (Step 3)')
    # --score_threshold
    group_filt.add_argument('--score_threshold', 
                            action='store',
                            metavar='REAL',
                            type=float, 
                            default=0.9,
                            help='Score threshold (real between 0 and 1). '
                                 'Default is %(default)s')
    # --straight_mode
    group_filt.add_argument('--straight_mode', 
                            action='store_true',
                            help='Use straight mode filtering. '
                                 'Default is geometric mode')
    
    # Step 4: Overlap Graph Building
    group_ovg = parser.add_argument_group('Overlap Graph Building (Step 4)')
    # --min_identity
    group_ovg.add_argument('--min_identity', 
                           action='store',
                           metavar='REAL',
                           type=float, 
                           default=1.0,
                           help='Minimum identity of an overlap between 2 reads. '
                                'Default is %(default)s')
    # --min_overlap_length
    group_ovg.add_argument('--min_overlap_length', 
                           action='store',
                           metavar='INT',
                           type=int, 
                           default=50,
                           help='Minimum length of an overlap. '
                                'Default is %(default)s')
    # --multi
    group_ovg.add_argument('--multi', 
                           action='store_true',
                           help='Use multi-ref mode')
    
    # Step 5: Graph Compaction & Components Identification
    group_gcomp = parser.add_argument_group('Graph Compaction & Components Identification (Step 5)')
    # -N / --min_read_node
    group_gcomp.add_argument('-N', '--min_read_node', 
                             action='store',
                             metavar='INT',
                             type=int, 
                             default=2,
                             help='Minimum number of read to keep a node. '
                                  'Default is %(default)s')
    # -E / --min_overlap_edge
    group_gcomp.add_argument('-E', '--min_overlap_edge', 
                             action='store',
                             metavar='INT',
                             type=int, 
                             default=20,
                             help='Minimum number of overlap to keep an edge. '
                                  'Default is %(default)s')
    
    # Step 6: LCA Labelling
    group_lca = parser.add_argument_group('LCA Labelling (Step 6)')
    # --quorum
    group_lca.add_argument('--quorum', 
                           action='store',
                           metavar='FLOAT',
                           type=float, 
                           default=0.51,
                           help='Quorum for LCA computing. Has to be between 0.51 and 1. '
                                'Default is %(default)s')
    
    # Step 7: Contig Assembly
    group_ass = parser.add_argument_group('Contig Assembly (Step 7)')
    
    # Step 8: Post Assembly Stats
    group_stats = parser.add_argument_group('Post Assembly Stats (Step 8)')
    # -bt / --blast_task
    group_stats.add_argument('-bt', '--blast_task',
                             action='store',
                             metavar='TASK',
                             type=str,
                             choices=['blastn', 'megablast'],
                             default='blastn',
                             help='Blast task (blastn or megablast). '
                                  'Default is %(default)s')
    # -be / --blast_evalue
    group_stats.add_argument('-be', '--blast_evalue',
                             action='store',
                             metavar='FLOAT',
                             type=float,
                             default=1e-5,
                             help='Blast evalue. '
                                  'Default is %(default)s')
    
    # Debug parameters
    group_debug = parser.add_argument_group('Debug')
    # --simulate_only
    group_debug.add_argument('--simulate_only', 
                             action='store_true',
                             help='Only output pipeline commands without executing them')
    # --test_dataset
    group_debug.add_argument('--test_dataset', 
                             action='store_true',
                             help='Compute additional stats when using a test dataset')
    
    #
    args = parser.parse_args()
    
    # Arguments checking
    
    # check if steps are in proper range
    steps_set = set(args.steps)
    if len(steps_set - set(range(0, 9))) > 0:
        parser.print_help()
        raise Exception("steps not in range [0,8]")
    
    if args.clustering_id_threshold < 0 or args.clustering_id_threshold > 1:
        parser.print_help()
        raise Exception("clustering id threshold not in range [0,1]")
    
    if args.score_threshold < 0 or args.score_threshold > 1:
        parser.print_help()
        raise Exception("score threshold not in range [0,1]")
    
    if args.min_identity < 0 or args.min_identity > 1:
        parser.print_help()
        raise Exception("min identity not in range [0,1]")
    
    if args.quorum < 0 or args.quorum > 1:
        parser.print_help()
        raise Exception("quorum not in range [0.51,1]")
    
    # Get absolute path for all arguments
    args.input_fastx = os.path.abspath(args.input_fastx)
    args.ref_db = os.path.abspath(args.ref_db)
    if args.output_contigs != 'DEFAULTNAME':
        args.output_contigs = os.path.abspath(args.output_contigs)
    args.db_wkdir = os.path.abspath(args.db_wkdir)
    args.wkdir = os.path.abspath(args.wkdir)
    
    #
    return args


def print_intro(args):
    """
    Print the introduction
    """
    
    sys.stdout.write("""
#################################
             MATAM
#################################\n\n""")
    
    # MATAM bin
    sys.stdout.write("""CMD: \
{matam} \
""".format(matam=matam_bin))
    
    # Verbose
    if args.verbose > 0:
        sys.stdout.write('-')
        for i in xrange(args.verbose):
            sys.stdout.write('v')
        sys.stdout.write(' ')
    
    # Debug
    if args.simulate_only:
        sys.stdout.write('--simulate_only ')
    if args.test_dataset:
        sys.stdout.write('--test_dataset ')
    
    # Performance
    sys.stdout.write("""\
--cpu {cpu} --max_memory {memory} \
""".format(cpu=args.cpu,
           memory=args.max_memory))
    
    # Step 0
    if args.remove_Ns:
        sys.stdout.write('--remove_Ns ')
    
    # Step 1
    sys.stdout.write("""\
--clustering_id_threshold {clustid} --kingdoms \
""".format(clustid=args.clustering_id_threshold))
    
    for kingdom in args.kingdoms:
        sys.stdout.write('{king} '.format(king=kingdom))
    
    if args.index_only:
        sys.stdout.write('--index_only ')
    
    # Step 2
    sys.stdout.write("""\
--best {best} --min_lis {min_lis} --evalue {evalue} \
""".format(best=args.best,
           min_lis=args.min_lis,
           evalue=args.evalue))
    
    # Step 3
    sys.stdout.write("""\
--score_threshold {scr_th} \
""".format(scr_th=args.score_threshold))

    if args.straight_mode:
        sys.stdout.write('--straight_mode')
    
    # Step 4
    sys.stdout.write("""\
--min_identity {min_id} --min_overlap_length {min_ov_lgth} \
""".format(min_id=args.min_identity,
           min_ov_lgth=args.min_overlap_length))
    
    if args.multi:
        sys.stdout.write('--multi ')
    
    # Step 5
    sys.stdout.write("""\
--min_read_node {N} --min_overlap_edge {E} \
""".format(N=args.min_read_node,
           E=args.min_overlap_edge))
    
    # Step 6
    sys.stdout.write("""\
--quorum {quorum} \
""".format(quorum=args.quorum))
    
    # Step 8
    sys.stdout.write("""\
--blast_task {btask} \
--blast_evalue {beval} \
""".format(btask=args.blast_task,
           beval=args.blast_evalue))
    
    # Main parameters
    sys.stdout.write("""\
--db_wkdir {dw} --wkdir {w} \
--steps \
""".format(dw=args.db_wkdir,
           w=args.wkdir))
    
    for step in args.steps:
        sys.stdout.write('{0} '.format(step))
    
    sys.stdout.write("""\
--input_fastx {i} --ref_db {d} --output_contigs {o} \n\n\
""".format(i=args.input_fastx,
           d=args.ref_db,
           o=args.output_contigs))
    
    # Important parameters, to write very evidently
    sys.stdout.write('PARAM: Input Fastx = {0}\n'.format(args.input_fastx))
    sys.stdout.write('PARAM: Ref DB      = {0}\n'.format(args.ref_db))
    sys.stdout.write('PARAM: Steps       = {0}\n'.format(args.steps))
    sys.stdout.write('\n')
    
    return 0


if __name__ == '__main__':

    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    
    # Arguments parsing
    args = parse_arguments()
    
    #
    print_intro(args)
    
    #
    ref_db_filename = args.ref_db.split('/')[-1]
    ref_db_basename = '.'.join(ref_db_filename.split('.')[:-1])
    
    input_fastx_filename = args.input_fastx.split('/')[-1]
    input_fastx_basename = '.'.join(input_fastx_filename.split('.')[:-1])
    input_fastx_extension = input_fastx_filename.split('.')[-1]
    
    #
    steps_set = frozenset(args.steps)
    
    ###############################
    # STEP 0: Ref DB pre-processing
    
    ref_db_taxo_filename = ref_db_basename + '.taxo.tab'
    ref_db_taxo_filepath = args.db_wkdir + '/' + ref_db_taxo_filename
    
    cleaned_ref_db_basename = ref_db_basename
    if args.remove_Ns:
        cleaned_ref_db_basename += '_noNs'
    else:
        cleaned_ref_db_basename += '_rdNs'
    cleaned_ref_db_filename = cleaned_ref_db_basename + '.cleaned.fasta'
    cleaned_ref_db_filepath = args.db_wkdir + '/' + cleaned_ref_db_filename
    
    if 0 in steps_set:
        sys.stdout.write('## Ref DB pre-processing step (0):\n\n')
        
        # Extract taxo from ref DB and sort by ref id
        cmd_line = extract_taxo_bin + ' -i ' + args.ref_db + ' | sort -k1,1 > '
        cmd_line += ref_db_taxo_filepath
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # Trim Ns from both sides
        # Filter out small sequences
        
        # Convert Us in Ts
        # Option: Either filter out seq with Ns or replace Ns with random nucl
        cmd_line = 'cat ' + args.ref_db
        cmd_line += ' | sed "/^>/!s/U/T/g" | sed "/^>/!s/u/t/g"'
        if not args.remove_Ns:
            cmd_line += ' | ' + replace_Ns_bin
        cmd_line += ' | ' + sort_fasta_bin + ' --reverse > '
        cmd_line += cleaned_ref_db_filepath
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
    
    
    ###########################
    # STEP 1: Ref DB clustering
    
    clustering_id_threshold_int = int(args.clustering_id_threshold * 100)
    
    sumaclust_kingdom_basename = ref_db_basename + '.sumaclust_'
    sumaclust_kingdom_basename += '{0}'.format(clustering_id_threshold_int) + '_by_kingdom'
    sumaclust_kingdom_filename = sumaclust_kingdom_basename + '.fasta'
    sumaclust_kingdom_filepath = args.db_wkdir + '/' + sumaclust_kingdom_filename
    
    clustered_ref_db_basename = cleaned_ref_db_basename + '_NR{0}'.format(clustering_id_threshold_int)
    clustered_ref_db_filename = clustered_ref_db_basename + '.fasta'
    clustered_ref_db_filepath = args.db_wkdir + '/' + clustered_ref_db_filename
    
    sortmerna_index_directory = args.db_wkdir + '/' + 'sortmerna_index'
    sortmerna_index_basepath = sortmerna_index_directory + '/' + clustered_ref_db_basename
    
    blast_db_directory = args.db_wkdir + '/' + 'blastdb'
    blast_db_basepath = blast_db_directory + '/' + clustered_ref_db_basename
    
    try:
        os.remove(sumaclust_kingdom_filepath)
    except OSError:
        pass
    
    if 1 in steps_set:
        sys.stdout.write('## Ref DB clustering step (1):\n\n')
        
        if not args.index_only:
            for kingdom in args.kingdoms:
                
                cleaned_ref_db_kingdom_basename = cleaned_ref_db_basename + '.' + kingdom
                cleaned_ref_db_kingdom_filename = cleaned_ref_db_kingdom_basename + '.fasta'
                cleaned_ref_db_kingdom_filepath = args.db_wkdir + '/' + cleaned_ref_db_kingdom_filename
                
                sumaclust_basename = cleaned_ref_db_kingdom_basename + '.sumaclust_' 
                sumaclust_basename += '{0}'.format(clustering_id_threshold_int)
                sumaclust_filename = sumaclust_basename + '.fasta'
                sumaclust_filepath = args.db_wkdir + '/' + sumaclust_filename
                
                sumaclust_centroids_filename = sumaclust_basename + '.centroids.fasta'
                sumaclust_centroids_filepath = args.db_wkdir + '/' + sumaclust_centroids_filename
                
                # Extracting kingdoms fasta files
                cmd_line = fasta_name_filter_bin + ' -i ' + cleaned_ref_db_filepath
                cmd_line += ' -s \' ' + kingdom + '\' > ' # !! need to be a space before the kingdom
                cmd_line += cleaned_ref_db_kingdom_filepath
                
                sys.stdout.write('CMD: {0}\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
                # Clutering with sumaclust
                # sumaclust -l is used to simulate semi-global alignment, since SumaClust is a global aligner
                sumaclust_cmd_line = sumaclust_bin + ' -l -t ' + '{0:.2f}'.format(args.clustering_id_threshold)
                sumaclust_cmd_line += ' -p ' + str(args.cpu)
                sumaclust_cmd_line += ' ' + cleaned_ref_db_kingdom_filepath
                sumaclust_cmd_line += ' > ' + sumaclust_filepath
                
                sys.stdout.write('CMD: {0}\n'.format(sumaclust_cmd_line))
                if not args.simulate_only:
                    subprocess.call(sumaclust_cmd_line, shell=True)
                
                ## Extracting centroids
                filter_cmd_line = fasta_name_filter_bin + ' -s "cluster_center=True" -i '
                filter_cmd_line += sumaclust_filepath + ' -o '
                filter_cmd_line += sumaclust_centroids_filepath
                
                sys.stdout.write('CMD: {0}\n'.format(filter_cmd_line))
                if not args.simulate_only:
                    subprocess.call(filter_cmd_line, shell=True)
                
                ## Concatenate kingdom centroids
                cmd_line = 'cat ' + sumaclust_centroids_filepath + ' >> '
                cmd_line += sumaclust_kingdom_filepath
                
                sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
            # Clean fasta headers
            clean_name_cmd_line = clean_name_bin + ' -i ' + sumaclust_kingdom_filepath
            clean_name_cmd_line += ' -o ' + clustered_ref_db_filepath
            
            sys.stdout.write('CMD: {0}\n\n'.format(clean_name_cmd_line))
            if not args.simulate_only:
                subprocess.call(clean_name_cmd_line, shell=True)
        
        # SortMeRNA Ref DB indexing
        try:
            if not os.path.exists(sortmerna_index_directory):
                os.makedirs(sortmerna_index_directory)
        except OSError:
            sys.stderr.write("\nERROR: {0} cannot be created\n\n".format(sortmerna_index_directory))
            raise
        
        indexdb_cmd_line = indexdb_bin + ' -v --ref ' + clustered_ref_db_filepath
        indexdb_cmd_line += ',' + sortmerna_index_basepath
        
        sys.stdout.write('CMD: {0}\n\n'.format(indexdb_cmd_line))
        if not args.simulate_only:
            subprocess.call(indexdb_cmd_line, shell=True)
        sys.stdout.write('\n')
        
        # Blast Ref DB indexing
        try:
            if not os.path.exists(blast_db_directory):
                os.makedirs(blast_db_directory)
        except OSError:
            sys.stderr.write("\nERROR: {0} cannot be created\n\n".format(blast_db_directory))
            raise
        
        cmd_line = 'makeblastdb -in ' + clustered_ref_db_filepath
        cmd_line += ' -dbtype nucl -out ' + blast_db_basepath
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ######################################
    # STEP 2: Reads mapping against Ref DB
    
    sortme_output_basename = input_fastx_basename 
    sortme_output_basename += '.sortmerna_vs_' + clustered_ref_db_basename
    sortme_output_basename += '_b' + str(args.best) + '_m' + str(args.min_lis)
    
    sortme_output_basepath = args.wkdir + '/' + sortme_output_basename
    
    sortme_output_fastx_filepath = sortme_output_basepath + '.' + input_fastx_extension
    
    if 2 in steps_set:
        sys.stdout.write('## Mapping step (2):\n\n')
        
        # Set SortMeRNA command line
        sortmerna_cmd_line = sortmerna_bin + ' --ref ' + clustered_ref_db_filepath
        sortmerna_cmd_line += ',' + sortmerna_index_basepath + ' --reads '
        sortmerna_cmd_line += args.input_fastx + ' --aligned ' + sortme_output_basepath
        sortmerna_cmd_line += ' --fastx --sam --blast "1 cigar qcov" --log --best '
        sortmerna_cmd_line += str(args.best) + ' --min_lis ' + str(args.min_lis) 
        sortmerna_cmd_line += ' -e {0:.2e}'.format(args.evalue)
        sortmerna_cmd_line += ' -a ' + str(args.cpu) + ' -v'
        
        # Run SortMeRNA
        sys.stdout.write('CMD: {0}\n\n'.format(sortmerna_cmd_line))
        if not args.simulate_only:
            subprocess.call(sortmerna_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    #############################
    # STEP 3: Alignment Filtering
    
    score_threshold_int = int(args.score_threshold * 100)
    
    sam_filt_basename = sortme_output_basename + '.scr_filt_'
    if args.straight_mode:
        sam_filt_basename += 'str_'
    else:
        sam_filt_basename += 'geo_'
    sam_filt_basename += str(score_threshold_int) + 'pct'
    sam_filt_filename = sam_filt_basename + '.sam'
    sam_filt_filepath = args.wkdir + '/' + sam_filt_filename
    
    sortme_output_filepath = sortme_output_basepath + '.sam'
    
    if 3 in steps_set:
        sys.stdout.write('## Alignment filtering step (3):\n\n')
        
        # Filtering scores command line
        filter_score_cmd_line = 'cat ' + sortme_output_filepath
        filter_score_cmd_line += ' | grep -v "^@" | sort -k 1,1V -k 12,12Vr'
        filter_score_cmd_line += ' | ' + filter_score_bin + ' -t ' + str(args.score_threshold)
        if not args.straight_mode:
            filter_score_cmd_line += ' --geometric'
        filter_score_cmd_line += ' > ' + sam_filt_filepath
        
        sys.stdout.write('CMD: {0}\n'.format(filter_score_cmd_line))
        if not args.simulate_only:
            subprocess.call(filter_score_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ################################
    # STEP 4: Overlap Graph Building
    
    min_identity_int = int(args.min_identity * 100)
    
    ovgraphbuild_basename = sam_filt_basename + '.ovgb_i' + str(min_identity_int)
    ovgraphbuild_basename += '_o' + str(args.min_overlap_length) + '_'
    if args.multi:
        ovgraphbuild_basename += 'multi'
    else:
        ovgraphbuild_basename += 'single'
    ovgraphbuild_basepath = args.wkdir + '/' + ovgraphbuild_basename
    
    ovgraphbuild_nodes_csv_filepath = ovgraphbuild_basepath + '.nodes.csv'
    ovgraphbuild_edges_csv_filepath = ovgraphbuild_basepath + '.edges.csv'
    
    if 4 in steps_set:
        sys.stdout.write('## Overlap Graph building step (4):\n\n')
        
        # Ovgraphbuild command line
        ovgraphbuild_cmd_line = ovgraphbuild_bin + ' -v --debug'
        ovgraphbuild_cmd_line += ' -i ' + str(args.min_identity)
        ovgraphbuild_cmd_line += ' -m ' + str(args.min_overlap_length)
        if args.multi:
            ovgraphbuild_cmd_line += ' --multi_ref'
        ovgraphbuild_cmd_line += ' --asqg --csv --output_basename '
        ovgraphbuild_cmd_line += ovgraphbuild_basepath
        ovgraphbuild_cmd_line += ' -r ' + clustered_ref_db_filepath
        ovgraphbuild_cmd_line += ' -s ' + sam_filt_filepath
        
        # Run ovgraphbuild
        sys.stdout.write('CMD: {0}\n\n'.format(ovgraphbuild_cmd_line))
        if not args.simulate_only:
            subprocess.call(ovgraphbuild_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    #######################################################
    # STEP 5: Graph Compaction & Components Identification
    
    componentsearch_basename = ovgraphbuild_basename + '.ctgs'
    componentsearch_basename += '_N' + str(args.min_read_node)
    componentsearch_basename += '_E' + str(args.min_overlap_edge)
    componentsearch_basepath = args.wkdir + '/' + componentsearch_basename
    
    if 5 in steps_set:
        sys.stdout.write('## Graph Compaction & Components Identification step (5):\n\n')
        
        # componentsearch command line
        componentsearch_cmd_line = 'java -Xmx' + str(args.max_memory) + 'M -cp "'
        componentsearch_cmd_line += componentsearch_jar + '" main.Main'
        componentsearch_cmd_line += ' -N ' + str(args.min_read_node)
        componentsearch_cmd_line += ' -E ' + str(args.min_overlap_edge)
        componentsearch_cmd_line += ' -b ' + componentsearch_basepath
        componentsearch_cmd_line += ' -n ' + ovgraphbuild_nodes_csv_filepath
        componentsearch_cmd_line += ' -e ' + ovgraphbuild_edges_csv_filepath
        
        # Run componentsearch
        sys.stdout.write('CMD: {0}\n\n'.format(componentsearch_cmd_line))
        if not args.simulate_only:
            subprocess.call(componentsearch_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    #######################
    # STEP 6: LCA Labelling
    
    quorum_int = int(args.quorum * 100)
    
    labelled_nodes_basename = componentsearch_basename + '.nodes_contracted'
    labelled_nodes_basename += '.component_lca' + str(quorum_int) + 'pct'
    
    components_lca_filename = componentsearch_basename + '.component_lca' + str(quorum_int) + 'pct.tab'
    components_lca_filepath = args.wkdir + '/' + components_lca_filename
    
    stats_filename = labelled_nodes_basename + '.stats'
    stats_filepath = args.wkdir + '/' + stats_filename
    
    if 6 in steps_set:
        sys.stdout.write('## LCA Labelling (6):\n\n')
        
        #
        cmd_line = 'tail -n +2 ' + componentsearch_basepath + '.components.csv'
        cmd_line += ' | sed "s/;/\\t/g" | sort -k2,2 > '
        cmd_line += componentsearch_basepath + '.components.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'tail -n +2 ' + ovgraphbuild_nodes_csv_filepath
        cmd_line += ' | sed "s/;/\\t/g" | sort -k1,1 > '
        cmd_line += ovgraphbuild_basepath + '.nodes.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        read_id_metanode_component_filepath = componentsearch_basepath + '.read_id_metanode_component.tab'
        complete_taxo_filepath = componentsearch_basepath + '.read_metanode_component_taxo.tab'
        
        cmd_line = 'join -a1 -e"NULL" -o "1.2,0,2.3,2.1" -11 -22 ' 
        cmd_line += ovgraphbuild_basepath + '.nodes.tab '
        cmd_line += componentsearch_basepath + '.components.tab '
        cmd_line += '| sed "s/ /\\t/g" | sort -k1,1  > '
        cmd_line += read_id_metanode_component_filepath
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + sam_filt_filepath + ' | cut -f1,3 | sort -k2,2'
        cmd_line += ' | join -12 -21 - ' + ref_db_taxo_filepath
        cmd_line += ' | sort -k2,2 | awk "{print \$2,\$3}" | sed "s/ /\\t/g" '
        cmd_line += ' | join -11 -21 ' + read_id_metanode_component_filepath
        cmd_line += ' - | sed "s/ /\\t/g" | cut -f2-5 > '
        cmd_line += complete_taxo_filepath
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + complete_taxo_filepath + ' | sort -k3,3 -k1,1 | '
        cmd_line += compute_lca_bin + ' -t 4 -f 3 -g 1 -m ' + str(args.quorum)
        cmd_line += ' -o ' + components_lca_filepath
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # Compressed graph Stats
        cmd_line = compute_compressed_graph_stats_bin + ' --nodes_contracted '
        cmd_line += componentsearch_basepath + '.nodes_contracted.csv --edges_contracted '
        cmd_line += componentsearch_basepath + '.edges_contracted.csv --components_lca '
        cmd_line += components_lca_filepath
        if args.test_dataset:
            cmd_line += ' --test_dataset'
            cmd_line += ' --species_taxo ' + args.db_wkdir + '/16sp.taxo.tab'
            cmd_line += ' --read_node_component ' + read_id_metanode_component_filepath
        cmd_line += ' -o ' + stats_filepath 
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ##########################
    # STEP 7: Components Assembly
    
    assembly_program = 'sga'
    
    components_assembly_basename = componentsearch_basename + '.'
    components_assembly_basename += assembly_program + '_by_component'
    components_assembly_filename = components_assembly_basename + '.fasta'
    components_assembly_filepath = args.wkdir + '/' + components_assembly_filename
    
    components_assembly_log_filename = components_assembly_basename + '.log'
    components_assembly_log_filepath = args.wkdir + '/' + components_assembly_log_filename
    
    contig_min_length = 500
    
    output_contigs_basename = components_assembly_basename + '.min_' 
    output_contigs_basename += str(contig_min_length) + 'bp'
    output_contigs_filename = output_contigs_basename + '.fasta'
    output_contigs_filepath = args.wkdir + '/' + output_contigs_filename
    
    if args.output_contigs == 'DEFAULTNAME':
        args.output_contigs = output_contigs_filepath
    
    assembly_wkdir = componentsearch_basepath + '.' + assembly_program
        
    if 7 in steps_set:
        sys.stdout.write('## Contigs Assembly (7):\n\n')
        
        # Generate the read_id-->component_id file
        component_read_filepath = componentsearch_basepath + '.component_read.tab'
        
        cmd_line = 'cat ' + componentsearch_basepath + '.components.tab'
        cmd_line += ' | join -12 -21 - ' + ovgraphbuild_basepath + '.nodes.tab'
        cmd_line += ' | awk \' {print $2"\\t"$4}\' | sort -k1,1n > ' 
        cmd_line += component_read_filepath
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # Convert input fastq to tab
        
        # Reading components LCA and storing them in a dict
        sys.stdout.write('INFO: Reading components LCA assignment from {0}\n\n'.format(components_lca_filepath))
        
        component_lca_dict = dict()
        
        with open(components_lca_filepath, 'r') as component_lca_fh:
            component_lca_dict = {t[0]: t[1] for t in (l.split() for l in component_lca_fh) if len(t)==2}
        
        # Prepare assembly wkdir and log file
        if os.path.exists(components_assembly_log_filepath):
            os.remove(components_assembly_log_filepath)
        
        try:
            if not os.path.exists(assembly_wkdir):
                os.makedirs(assembly_wkdir)
        except OSError:
            sys.stderr.write("\nERROR: {0} cannot be created\n\n".format(assembly_wkdir))
            raise
        
        # Open output contigs file
        components_assembly_fh = open(components_assembly_filepath, 'w')
        
        # Begin reading read-->component file (sorted by component)
        contig_count = 0
        component_count = 0
        
        with open(component_read_filepath, 'r') as component_read_fh:
            # We will be working in the assembly wkdir because assembly tools 
            # generate lots of files we dont want in our matam wkdir
            previous_wkdir = os.getcwd()
            os.chdir(assembly_wkdir)
            
            # Assembling all reads of every component, one at a time
            for component_tab_list in read_tab_file_handle_sorted(component_read_fh, 0):
                component_count += 1
                component_id = component_tab_list[0][0]
                
                # To prevent assembly of singleton reads
                if component_id == 'NULL':
                    continue
                
                # Starting component assembly
                sys.stdout.write('\rINFO: Assembling component #{0}'.format(component_count))
                
                # Cleaning previous assembly
                if os.path.exists('contigs.fa'):
                    os.remove('contigs.fa')
                if os.path.exists('tmp'):
                    subprocess.call('rm -rf tmp', shell=True)
                
                # Write the component reads ids
                if not args.simulate_only:
                    with open('reads_single_component.ids', 'w') as wfh:
                        for tab in component_tab_list:
                            read_id = tab[1]
                            wfh.write('{0}\n'.format(read_id))
                
                # Generate a fastq file with this component reads
                cmd_line = fastq_name_filter_bin + ' -f reads_single_component.ids'
                cmd_line += ' -i ' + sortme_output_fastx_filepath + ' -o reads_single_component.fq'
                
                #~ sys.stdout.write('CMD: {0}\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
                # Assemble those reads with SGA
                cmd_line = 'echo "component #' + component_id + '" >> '
                cmd_line += components_assembly_log_filepath + ' && '
                cmd_line += sga_assemble_bin + ' -i reads_single_component.fq'
                cmd_line += ' -o contigs.fa --sga_bin ' + sga_bin
                cmd_line += ' --cpu ' + str(args.cpu)
                cmd_line += ' --tmp_dir tmp'
                cmd_line += ' >> ' + components_assembly_log_filepath + ' 2>&1'
                
                #~ sys.stdout.write('CMD: {0}\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
                # Concatenate the component contigs in the output contigs file
                if not args.simulate_only:
                    component_lca = ''
                    if component_id in component_lca_dict:
                        component_lca = component_lca_dict[component_id]
                    with open('contigs.fa', 'r') as contigs_fh:
                        for header, seq in read_fasta_file_handle(contigs_fh):
                            if len(seq):
                                contig_count += 1
                                components_assembly_fh.write('>{0} component={1} '.format(contig_count, component_id))
                                components_assembly_fh.write('lca={0}\n{1}\n'.format(component_lca, format_seq(seq)))
            
            # Return to previous directory
            os.chdir(previous_wkdir)
        
        # Close assembly contigs file
        components_assembly_fh.close()
        
        # Filter assembly by length
        cmd_line = fasta_length_filter_bin + ' -m ' + str(contig_min_length)
        cmd_line += ' -i ' + components_assembly_filepath
        cmd_line += ' -o ' + args.output_contigs
        
        sys.stdout.write('\n\nCMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        sys.stdout.write('INFO: Assembly contigs in {0}\n\n'.format(args.output_contigs))
    
    #############################
    # STEP 8: Post assembly Stats
    
    max_target_seqs = 10000
    
    blast_output_basename = components_assembly_basename + '.' + args.blast_task + '_vs_'
    blast_output_basename += clustered_ref_db_basename + '.max_target_seqs_'
    blast_output_basename += str(max_target_seqs)
    blast_output_basepath = args.wkdir + '/' + blast_output_basename
    
    blast_output_filename = blast_output_basename + '.tab'
    blast_output_filepath = args.wkdir + '/' + blast_output_filename
    
    if 8 in steps_set:
        sys.stdout.write('## Post assembly Stats (8):\n\n')
        
        #
        cmd_line = 'blastn -query ' + components_assembly_filepath
        cmd_line += ' -task ' + args.blast_task + ' -db ' + blast_db_basepath
        cmd_line += ' -out ' + blast_output_filepath
        cmd_line += ' -evalue ' + str(args.blast_evalue)
        cmd_line += ' -outfmt "6 std qlen slen" -dust "no"'
        cmd_line += ' -max_target_seqs ' + str(max_target_seqs)
        cmd_line += ' -num_threads ' + str(args.cpu)
            
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'sort -k2,2 ' + blast_output_filepath
        cmd_line += ' > ' + blast_output_basepath + '.sorted_subject.tab'
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = find_least_bin + ' -i ' + blast_output_filepath
        cmd_line += ' -s ' + blast_output_basepath + '.sorted_subject.tab'
        cmd_line += ' -o ' + blast_output_basepath + '.least_num_ref.tab'
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # When using a test dataset, remapping with blastn against the original sequences
        if args.test_dataset:
            # Indexing original sequences
            cmd_line = 'makeblastdb -in ' + args.db_wkdir + '/16sp.fasta -dbtype nucl '
            cmd_line += '-out ' + blast_db_directory + '/16sp'
            
            sys.stdout.write('\nCMD: {0}'.format(cmd_line))
            if not args.simulate_only:
                subprocess.call(cmd_line, shell=True)
            sys.stdout.write('\n')
            
            # Blast assembly contigs against original sequences
            
            test_blast_output_filename = components_assembly_basename + '.' + args.blast_task + '_vs_16sp'
            #~ test_blast_output_filename += '.max_target_seqs_' + str(max_target_seqs) + '.tab'
            test_blast_output_filename += '.max_target_seqs_1.tab'
            test_blast_output_filepath = args.wkdir + '/' + test_blast_output_filename
            
            cmd_line = 'blastn -query ' + components_assembly_filepath
            cmd_line += ' -task ' + args.blast_task + ' -db ' + blast_db_directory + '/16sp'
            cmd_line += ' -out ' + test_blast_output_filepath
            cmd_line += ' -evalue ' + str(args.blast_evalue)
            cmd_line += ' -outfmt "6 std qlen slen" -dust "no"'
            #~ cmd_line += ' -max_target_seqs ' + str(max_target_seqs)
            cmd_line += ' -max_target_seqs 1'
            cmd_line += ' -num_threads ' + str(args.cpu)
            
            sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
            if not args.simulate_only:
                subprocess.call(cmd_line, shell=True)
    
    exit(0)

