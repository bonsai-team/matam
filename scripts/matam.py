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
fastq_name_filter_bin = matam_dir + '/bin/fastq_name_filter.py'
indexdb_bin = matam_dir + '/bin/indexdb_rna'
sortmerna_bin = matam_dir + '/bin/sortmerna'
sga_bin = matam_dir + '/bin/sga'
filter_score_bin = matam_dir + '/bin/filter_score_multialign.py'
ovgraphbuild_bin = matam_dir + '/bin/ovgraphbuild'
contigsearch_jar = matam_dir + '/bin/ContigSearch.jar'
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
        sys.stderr.write('\nerror: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def parse_arguments():
    """
    Parse the command line, and check if arguments are correct
    """
    #
    parser = DefaultHelpParser(description='matam',
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
                            type=int, nargs='+', default=[0,1,2,3,4,5,6,7,8],
                            help='Steps to execute')
    #
    group_perf = parser.add_argument_group('Performance')
    group_perf.add_argument('--cpu', metavar='INT',
                            type=int, default=3,
                            help='Max number of CPU to use')
    group_perf.add_argument('--max_memory', metavar='INT',
                            type=int, default=4000,
                            help='Maximum memory to use (in MBi). Default is 4000MBi')
    #
    group_db = parser.add_argument_group('Ref DB pre-processing (Step 0)')
    group_db.add_argument('--remove_Ns', action='store_true',
                          help='Remove sequences with Ns. (Default is replacing Ns with random nucleotides)')
    #
    group_clust = parser.add_argument_group('Clustering Ref DB (Step 1)')
    group_clust.add_argument('--clustering_id_threshold', metavar='REAL',
                             type=float, default=0.95,
                             help='Identity threshold for clustering')
    group_clust.add_argument('--kingdoms', metavar='STR',
                             type=str, nargs='+', default=['archaea','bacteria','eukaryota'],
                             help='Kingdoms to clusterize the DB for. Default is: archaea bacteria eukaryota')
    group_clust.add_argument('--index_only', action='store_true',
                             help=argparse.SUPPRESS)
    #
    group_mapping = parser.add_argument_group('Mapping (Step 2)')
    group_mapping.add_argument('--best', metavar='INT',
                               type=int, default=10,
                               help='Get up to --best good alignments per read')
    group_mapping.add_argument('--min_lis', metavar='INT',
                               type=int, default=10,
                               help=argparse.SUPPRESS)
    group_mapping.add_argument('--evalue', metavar='REAL',
                               type=float, default=1e-5,
                               help='Max e-value to keep an alignment for (default: 1e-5)')
    #
    group_filt = parser.add_argument_group('Alignment Filtering (Step 3)')
    group_filt.add_argument('--score_threshold', metavar='REAL',
                            type=float, default=0.9,
                            help='Score threshold (real between 0 and 1)')
    group_filt.add_argument('--straight_mode', action='store_true',
                            help='Use straight mode filtering. Default is geometric mode')
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
    group_debug = parser.add_argument_group('Debug')
    group_debug.add_argument('--simulate_only', action='store_true',
                             help='Only output pipeline commands without executing them')
    group_debug.add_argument('--test_dataset', action='store_true',
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
    
    sys.stdout.write("""CMD: \
{matam} --cpu {cpu} --max_memory {memory} \
""".format(matam=matam_bin,
           cpu=args.cpu,
           memory=args.max_memory))
    
    if args.remove_Ns:
        sys.stdout.write('--remove_Ns ')
    
    sys.stdout.write("""\
--clustering_id_threshold {clustid} --kingdoms \
""".format(clustid=args.clustering_id_threshold))
    
    for kingdom in args.kingdoms:
        sys.stdout.write('{king} '.format(king=kingdom))
    
    if args.index_only:
        sys.stdout.write('--index_only ')
    
    sys.stdout.write("""\
--best {best} --min_lis {min_lis} --evalue {evalue} \
--score_threshold {scr_th} \
""".format(best=args.best,
           min_lis=args.min_lis,
           evalue=args.evalue,
           scr_th=args.score_threshold))
    
    if args.straight_mode:
        sys.stdout.write('--straight_mode')
    
    sys.stdout.write("""\
--min_identity {min_id} --min_overlap_length {min_ov_lgth} \
""".format(min_id=args.min_identity,
           min_ov_lgth=args.min_overlap_length))
    
    if args.multi:
        sys.stdout.write('--multi ')
    
    sys.stdout.write("""\
--min_read_node {N} --min_overlap_edge {E} \
--quorum {quorum} --steps \
""".format(N=args.min_read_node,
           E=args.min_overlap_edge,
           quorum=args.quorum))
    
    for step in args.steps:
        sys.stdout.write('{0} '.format(step))
    
    sys.stdout.write("""\
--input_fastx {i} --ref_db {d} --output_contigs {o} \n\n\
""".format(i=args.input_fastx,
           d=args.ref_db,
           o=args.output_contigs))
    
    #~ sys.stdout.write('PARAM: Executable  = {0}\n'.format(os.path.abspath(os.path.dirname(os.path.realpath(sys.argv[0])))))
    sys.stdout.write('PARAM: Matam dir = {0}\n'.format(matam_dir))
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
    
    #
    steps_set = frozenset(args.steps)
    
    ###############################
    # STEP 0: Ref DB pre-processing
    cleaned_ref_db_basename = ref_db_basename
    if args.remove_Ns:
        cleaned_ref_db_basename += '_noNs'
    else:
        cleaned_ref_db_basename += '_rdNs'
    ref_db_taxo_filename = ref_db_basename + '.taxo.tab'
    cleaned_ref_db_filename = cleaned_ref_db_basename + '.cleaned.fasta'
    
    if 0 in steps_set:
        sys.stdout.write('## Ref DB pre-processing step (0):\n\n')
        
        # Extract taxo from ref DB and sort by ref id
        cmd_line = extract_taxo_bin + ' -i ' + args.ref_db + ' | sort -k1,1 > '
        cmd_line += ref_db_taxo_filename
        
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
        cmd_line += cleaned_ref_db_filename
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
    
    
    ###########################
    # STEP 1: Ref DB clustering
    cluster_id_int = int(args.clustering_id_threshold * 100)
    clustered_ref_db_basename = cleaned_ref_db_basename + '_NR{0}'.format(cluster_id_int)
    sortmerna_index_directory = 'sortmerna_index'
    sortmerna_index_basename = sortmerna_index_directory + '/' + clustered_ref_db_basename
    blast_db_directory = 'blastdb'
    blast_db_basename = blast_db_directory + '/' + clustered_ref_db_basename
    
    sumaclust_kingdom_filename = ref_db_basename + '.sumaclust_'
    sumaclust_kingdom_filename += '{0}'.format(cluster_id_int)
    sumaclust_kingdom_filename += '_by_kingdom.fasta'
    
    try:
        os.remove(sumaclust_kingdom_filename)
    except OSError:
        pass
    
    if 1 in steps_set:
        sys.stdout.write('## Ref DB clustering step (1):\n\n')
        
        if not args.index_only:
            for kingdom in args.kingdoms:
                cleaned_ref_db_kingdom_basename = cleaned_ref_db_basename + '.' + kingdom
                
                # Extracting kingdoms fasta files
                cmd_line = fasta_name_filter_bin + ' -i ' + cleaned_ref_db_filename
                cmd_line += ' -s \' ' + kingdom + '\' > ' # !! need to be a space before the kingdom
                cmd_line += cleaned_ref_db_kingdom_basename + '.fasta'
                
                sys.stdout.write('CMD: {0}\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
                # Clutering with sumaclust
                sumaclust_basename = cleaned_ref_db_kingdom_basename + '.sumaclust_' 
                sumaclust_basename += '{0}'.format(cluster_id_int)
                
                # sumaclust -l is used to simulate semi-global alignment, since SumaClust is a global aligner
                sumaclust_cmd_line = sumaclust_bin + ' -l -t ' + '{0:.2f}'.format(args.clustering_id_threshold)
                sumaclust_cmd_line += ' -p ' + str(args.cpu)
                sumaclust_cmd_line += ' ' + cleaned_ref_db_kingdom_basename + '.fasta'
                sumaclust_cmd_line += ' > ' + sumaclust_basename + '.fasta'
                
                sys.stdout.write('CMD: {0}\n'.format(sumaclust_cmd_line))
                if not args.simulate_only:
                    subprocess.call(sumaclust_cmd_line, shell=True)
                
                ## Extracting centroids
                filter_cmd_line = fasta_name_filter_bin + ' -s "cluster_center=True" -i '
                filter_cmd_line += sumaclust_basename + '.fasta -o '
                filter_cmd_line += sumaclust_basename + '.centroids.fasta'
                
                sys.stdout.write('CMD: {0}\n'.format(filter_cmd_line))
                if not args.simulate_only:
                    subprocess.call(filter_cmd_line, shell=True)
                
                ## Concatenate kingdom centroids
                cmd_line = 'cat ' + sumaclust_basename + '.centroids.fasta >> '
                cmd_line += sumaclust_kingdom_filename
                
                sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
            # Clean fasta headers
            clean_name_cmd_line = clean_name_bin + ' -i ' + sumaclust_kingdom_filename
            clean_name_cmd_line += ' -o ' + clustered_ref_db_basename + '.fasta'
            
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
        
        indexdb_cmd_line = indexdb_bin + ' -v --ref ' + clustered_ref_db_basename
        indexdb_cmd_line += '.fasta,' + sortmerna_index_basename
        
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
        
        cmd_line = 'makeblastdb -in ' + clustered_ref_db_basename + '.fasta' 
        cmd_line += ' -dbtype nucl -out ' + blast_db_basename
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ######################################
    # STEP 2: Reads mapping against Ref DB
    sortme_output_basename = input_fastx_basename + '.sortmerna_vs_' + clustered_ref_db_basename
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
    
    if 3 in steps_set:
        sys.stdout.write('## Alignment filtering step (3):\n\n')
        
        # Filtering scores command line
        filter_score_cmd_line = 'cat ' + sortme_output_basename + '.sam'
        filter_score_cmd_line += ' | grep -v "^@" | sort -k 1,1V -k 12,12Vr'
        filter_score_cmd_line += ' | ' + filter_score_bin + ' -t ' + str(args.score_threshold)
        if not args.straight_mode:
            filter_score_cmd_line += ' --geometric'
        filter_score_cmd_line += ' > ' + sam_filt_basename + '.sam'
        
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
        ovgraphbuild_cmd_line += ' -r ' + clustered_ref_db_basename + '.fasta'
        ovgraphbuild_cmd_line += ' -s ' + sam_filt_basename + '.sam'
        
        # Run ovgraphbuild
        sys.stdout.write('CMD: {0}\n\n'.format(ovgraphbuild_cmd_line))
        if not args.simulate_only:
            subprocess.call(ovgraphbuild_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ##################################################
    # STEP 5: Graph Compaction & Contig Identification
    contigsearch_basename = ovgraphbuild_basename + '.ctgs'
    contigsearch_basename += '_N' + str(args.min_read_node)
    contigsearch_basename += '_E' + str(args.min_overlap_edge)
    
    if 5 in steps_set:
        sys.stdout.write('## Graph Compaction & Contig Identification step (5):\n\n')
        
        # ContigSearch command line
        contigsearch_cmd_line = 'java -Xmx' + str(args.max_memory) + 'M -cp "'
        contigsearch_cmd_line += contigsearch_jar + '" main.Main'
        contigsearch_cmd_line += ' -N ' + str(args.min_read_node)
        contigsearch_cmd_line += ' -E ' + str(args.min_overlap_edge)
        contigsearch_cmd_line += ' -b ' + contigsearch_basename
        contigsearch_cmd_line += ' -n ' + ovgraphbuild_basename + '.nodes.csv'
        contigsearch_cmd_line += ' -e ' + ovgraphbuild_basename + '.edges.csv'
        
        # Run ContigSearch
        sys.stdout.write('CMD: {0}\n\n'.format(contigsearch_cmd_line))
        if not args.simulate_only:
            subprocess.call(contigsearch_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    #######################
    # STEP 6: LCA Labelling
    quorum_int = int(args.quorum * 100)
    labelled_nodes_basename = contigsearch_basename + '.nodes_contracted'
    labelled_nodes_basename += '.component_lca' + str(quorum_int) + 'pct'
    components_lca_filename = contigsearch_basename + '.component_lca' + str(quorum_int) + 'pct.tab'
    
    if 6 in steps_set:
        sys.stdout.write('## LCA Labelling (6):\n\n')
        
        #
        cmd_line = 'tail -n +2 ' + contigsearch_basename + '.contigs.csv'
        cmd_line += ' | sed "s/;/\\t/g" | sort -k2,2 > '
        cmd_line += contigsearch_basename + '.contigs.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'tail -n +2 ' + ovgraphbuild_basename + '.nodes.csv'
        cmd_line += ' | sed "s/;/\\t/g" | sort -k1,1 > '
        cmd_line += ovgraphbuild_basename + '.nodes.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        read_id_metanode_component_filename = contigsearch_basename + '.read_id_metanode_component.tab'
        complete_taxo_filename = contigsearch_basename + '.read_metanode_component_taxo.tab'
        
        cmd_line = 'join -a1 -e"NULL" -o "1.2,0,2.3,2.1" -11 -22 ' 
        cmd_line += ovgraphbuild_basename + '.nodes.tab '
        cmd_line += contigsearch_basename + '.contigs.tab '
        cmd_line += '| sed "s/ /\\t/g" | sort -k1,1  > '
        cmd_line += read_id_metanode_component_filename
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + sam_filt_basename + '.sam | cut -f1,3 | sort -k2,2'
        cmd_line += ' | join -12 -21 - ' + ref_db_basename + '.taxo.tab'
        cmd_line += ' | sort -k2,2 | awk "{print \$2,\$3}" | sed "s/ /\\t/g" '
        cmd_line += ' | join -11 -21 ' + read_id_metanode_component_filename
        cmd_line += ' - | sed "s/ /\\t/g" | cut -f2-5 > '
        cmd_line += complete_taxo_filename
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        #~ cmd_line = 'cat ' + complete_taxo_filename + ' | sort -k2,2 | '
        #~ cmd_line += compute_lca_bin + ' -t 4 -f 2 -g 1 -m ' + str(args.quorum) + ' > '
        #~ cmd_line += contigsearch_basename + '.metanode_lca' + str(quorum_int) + 'pct.tab'
        #~ 
        #~ sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        #~ if not args.simulate_only:
            #~ subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + complete_taxo_filename + ' | sort -k3,3 -k1,1 | '
        cmd_line += compute_lca_bin + ' -t 4 -f 3 -g 1 -m ' + str(args.quorum)
        cmd_line += ' -o ' + components_lca_filename
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # Compressed graph Stats
        stats_filename = labelled_nodes_basename + '.stats'
        
        cmd_line = compute_compressed_graph_stats_bin + ' --nodes_contracted '
        cmd_line += contigsearch_basename + '.nodes_contracted.csv --edges_contracted '
        cmd_line += contigsearch_basename + '.edges_contracted.csv --components_lca '
        cmd_line += components_lca_filename
        if args.test_dataset:
            cmd_line += ' --test_dataset'
            cmd_line += ' --species_taxo 16sp.taxo.tab'
            cmd_line += ' --read_node_component ' + read_id_metanode_component_filename
        cmd_line += ' -o ' + stats_filename 
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ##########################
    # STEP 7: Contigs Assembly
    if args.output_contigs == 'DEFAULTNAME':
        args.output_contigs = contigsearch_basename + '.sga_by_component.fasta'
    sga_log_filename = contigsearch_basename + '.sga_by_component.log'
    sga_wkdir = contigsearch_basename + '.sga'
    
    if 7 in steps_set:
        sys.stdout.write('## Contigs Assembly (7):\n\n')
        
        # Generate the read_id-->component_id file
        contig_read_filename = contigsearch_basename + '.component_read.tab'
        
        cmd_line = 'cat ' + contigsearch_basename + '.contigs.tab'
        cmd_line += ' | join -12 -21 - ' + ovgraphbuild_basename + '.nodes.tab'
        cmd_line += ' | awk \' {print $2"\\t"$4}\' | sort -k1,1n > ' 
        cmd_line += contig_read_filename
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # Convert input fastq to tab
        
        # Reading components LCA and storing them in a dict
        sys.stdout.write('INFO: Reading components LCA assignment from {0}\n\n'.format(components_lca_filename))
        
        contig_lca_dict = dict()
        
        with open(components_lca_filename, 'r') as contig_lca_fh:
            contig_lca_dict = {t[0]: t[1] for t in (l.split() for l in contig_lca_fh) if len(t)==2}
        
        # Prepare sga wkdir and log file
        if os.path.exists(sga_log_filename):
            os.remove(sga_log_filename)
        
        try:
            if not os.path.exists(sga_wkdir):
                os.makedirs(sga_wkdir)
        except OSError:
            sys.stderr.write("\nERROR: {0} cannot be created\n\n".format(sga_wkdir))
            raise
        
        # Open output contigs file
        output_fh = open(args.output_contigs, 'w')
        
        # Begin reading read-->component file (sorted by component)
        contig_count = 0
        component_count = 0
        
        with open(contig_read_filename, 'r') as contig_read_fh:
            # We will be working in the sga wkdir because sga generates
            # lots of files we dont want in our matam wkdir
            os.chdir(sga_wkdir)
            
            # Assembling all reads of every component, one at a time
            for contig_tab_list in read_tab_file_handle_sorted(contig_read_fh, 0):
                component_count += 1
                component_id = contig_tab_list[0][0]
                
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
                    with open('reads_one_contig.ids', 'w') as wfh:
                        for tab in contig_tab_list:
                            read_id = tab[1]
                            wfh.write('{0}\n'.format(read_id))
                
                # Generate a fastq file with this component reads
                cmd_line = fastq_name_filter_bin + ' -f reads_one_contig.ids'
                cmd_line += ' -i ../' + sortme_output_basename + '.fastq'
                cmd_line += ' -o reads_one_contig.fq'
                
                #~ sys.stdout.write('CMD: {0}\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
                # Assemble those reads with SGA
                cmd_line = 'echo "component #' + component_id + '" >> ../'
                cmd_line += sga_log_filename + ' && '
                cmd_line += sga_assemble_bin + ' -i reads_one_contig.fq'
                cmd_line += ' -o contigs.fa --sga_bin ' + sga_bin
                cmd_line += ' --cpu ' + str(args.cpu)
                cmd_line += ' --tmp_dir tmp'
                cmd_line += ' >> ../' + sga_log_filename + ' 2>&1'
                
                #~ sys.stdout.write('CMD: {0}\n'.format(cmd_line))
                if not args.simulate_only:
                    subprocess.call(cmd_line, shell=True)
                
                # Concatenate the component contigs in the output contigs file
                if not args.simulate_only:
                    component_lca = ''
                    if component_id in contig_lca_dict:
                        component_lca = contig_lca_dict[component_id]
                    with open('contigs.fa', 'r') as sga_contigs_fh:
                        for header, seq in read_fasta_file_handle(sga_contigs_fh):
                            if len(seq):
                                contig_count += 1
                                output_fh.write('>{0} component={1} '.format(contig_count, component_id))
                                output_fh.write('lca={0}\n{1}\n'.format(component_lca, format_seq(seq)))
            
        sys.stdout.write('\n\nINFO: Assembly contigs in {0}\n\n'.format(args.output_contigs))
        
        # Return to upper directory
        os.chdir('..')
    
    #############################
    # STEP 8: Post assembly Stats
    contig_assembly_filename = args.output_contigs.split('/')[-1]
    contig_assembly_basename = '.'.join(contig_assembly_filename.split('.')[:-1])
    max_target_seqs = 100000
    #~ blast_task = 'megablast'
    blast_task = 'blastn'
    
    blast_output_basename = contig_assembly_basename + '.' + blast_task + '_vs_'
    blast_output_basename += clustered_ref_db_basename + '.max_target_seqs_'
    blast_output_basename += str(max_target_seqs)
    
    if 8 in steps_set:
        sys.stdout.write('## Post assembly Stats (8):\n\n')
        
        blast_output_filename = blast_output_basename + '.tab'
        
        #
        cmd_line = 'blastn -query ' + args.output_contigs
        cmd_line += ' -task ' + blast_task + ' -db ' + blast_db_basename
        cmd_line += ' -out ' + blast_output_filename + ' -evalue 1e-5'
        cmd_line += ' -outfmt "6 std qlen slen" -dust "no"'
        cmd_line += ' -max_target_seqs ' + str(max_target_seqs)
        cmd_line += ' -num_threads ' + str(args.cpu)
            
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'sort -k2,2 ' + blast_output_filename
        cmd_line += ' > ' + blast_output_basename + '.sorted_subject.tab'
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = find_least_bin + ' -i ' + blast_output_filename
        cmd_line += ' -s ' + blast_output_basename + '.sorted_subject.tab'
        cmd_line += ' -o ' + blast_output_basename + '.least_num_ref.tab'
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        if not args.simulate_only:
            subprocess.call(cmd_line, shell=True)
        
        # When using a test dataset, remapping with blastn against the original sequences
        if args.test_dataset:
            # Indexing original sequences
            cmd_line = 'makeblastdb -in 16sp.fasta -dbtype nucl '
            cmd_line += '-out ' + blast_db_directory + '/16sp'
            
            sys.stdout.write('CMD: {0}'.format(cmd_line))
            if not args.simulate_only:
                subprocess.call(cmd_line, shell=True)
            sys.stdout.write('\n')
            
            # Blast assembly contigs against original sequences
            
            test_blast_output_filename = contig_assembly_basename + '.' + blast_task + '_vs_16sp'
            #~ test_blast_output_filename += '.max_target_seqs_' + str(max_target_seqs) + '.tab'
            test_blast_output_filename += '.max_target_seqs_1.tab'
            
            cmd_line = 'blastn -query ' + args.output_contigs
            cmd_line += ' -task ' + blast_task + ' -db ' + blast_db_directory + '/16sp'
            cmd_line += ' -out ' + test_blast_output_filename + ' -evalue 1e-5'
            cmd_line += ' -outfmt "6 std qlen slen" -dust "no"'
            #~ cmd_line += ' -max_target_seqs ' + str(max_target_seqs)
            cmd_line += ' -max_target_seqs 1'
            cmd_line += ' -num_threads ' + str(args.cpu)
            
            sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
            if not args.simulate_only:
                subprocess.call(cmd_line, shell=True)
    
    exit(0)

