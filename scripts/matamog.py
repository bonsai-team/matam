#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
import subprocess

matamog_bin = os.path.realpath(sys.argv[0])
matamog_dir = matamog_bin[:-19]
sumaclust_bin = matamog_dir + '/bin/sumaclust'
extract_taxo_bin = matamog_dir + '/bin/extract_taxo_from_fasta.py'
clean_name_bin = matamog_dir + '/bin/fasta_clean_name.py'
name_filter_bin = matamog_dir + '/bin/fasta_name_filter.py'
indexdb_bin = matamog_dir + '/bin/indexdb_rna'
sortmerna_bin = matamog_dir + '/bin/sortmerna'
filter_score_bin = matamog_dir + '/bin/filter_score_multialign.py'
ovgraphbuild_bin = matamog_dir + '/bin/ovgraphbuild'
contigsearch_jar = matamog_dir + '/bin/ContigSearch.jar'
compute_lca_bin = matamog_dir + '/bin/compute_lca_from_tab.py'
compute_stats_lca_bin = matamog_dir + '/bin/compute_stats_from_lca.py'
replace_Ns_bin = matamog_dir + '/bin/replace_Ns_by_rand_nu.py'
sort_fasta_bin = matamog_dir + '/bin/sort_fasta_by_length.py'


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
                            type=int, nargs='+', default=[0,1,2,3,4,5,6,7],
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
    args = parser.parse_args()
    
    # Arguments checking
    
    # check if steps are in proper range
    steps_set = set(args.steps)
    if len(steps_set - set(range(0, 8))) > 0:
        raise Exception("steps not in range [0,7]")
    
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
#################################\n\n""")
    
    sys.stdout.write("""CMD: \
{matamog} --cpu {cpu} --max_memory {memory} \
""".format(matamog=matamog_bin,
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
    sys.stdout.write('PARAM: Matamog dir = {0}\n'.format(matamog_dir))
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
        subprocess.call(cmd_line, shell=True)
    
    ###########################
    # STEP 1: Ref DB clustering
    cluster_id_int = int(args.clustering_id_threshold * 100)
    clustered_ref_db_basename = cleaned_ref_db_basename + '_NR{0}'.format(cluster_id_int)
    sortmerna_index_directory = 'sortmerna_index'
    sortmerna_index_basename = sortmerna_index_directory + '/' + clustered_ref_db_basename
    
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
                cmd_line = name_filter_bin + ' -i ' + cleaned_ref_db_filename
                cmd_line += ' -s \' ' + kingdom + '\' > ' # !! need to be a space before the kingdom
                cmd_line += cleaned_ref_db_kingdom_basename + '.fasta'
                
                sys.stdout.write('CMD: {0}\n'.format(cmd_line))
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
                subprocess.call(sumaclust_cmd_line, shell=True)
                
                ## Extracting centroids
                filter_cmd_line = name_filter_bin + ' -s "cluster_center=True" -i '
                filter_cmd_line += sumaclust_basename + '.fasta -o '
                filter_cmd_line += sumaclust_basename + '.centroids.fasta'
                
                sys.stdout.write('CMD: {0}\n'.format(filter_cmd_line))
                subprocess.call(filter_cmd_line, shell=True)
                
                ## Concatenate kingdom centroids
                cmd_line = 'cat ' + sumaclust_basename + '.centroids.fasta >> '
                cmd_line += sumaclust_kingdom_filename
                
                sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
                subprocess.call(cmd_line, shell=True)
                
            # Clean fasta headers
            clean_name_cmd_line = clean_name_bin + ' -i ' + sumaclust_kingdom_filename
            clean_name_cmd_line += ' -o ' + clustered_ref_db_basename + '.fasta'
            
            sys.stdout.write('CMD: {0}\n\n'.format(clean_name_cmd_line))
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
        subprocess.call(indexdb_cmd_line, shell=True)
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
    read_ref_taxo_basename = sam_filt_basename + '.read_ref_taxo'
    
    if 3 in steps_set:
        sys.stdout.write('## Alignment filtering step (3):\n\n')
        
        # Filtering scores command line
        filter_score_cmd_line = 'cat ' + sortme_output_basename + '.sam'
        filter_score_cmd_line += ' | grep -v "^@" | sort -k 1,1V -k 12,12Vr'
        filter_score_cmd_line += ' | ' + filter_score_bin + ' -t ' + str(args.score_threshold)
        if not args.straight_mode:
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
        
        sys.stdout.write('CMD: {0}\n\n'.format(read_taxo_cmd_line))
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
        ovgraphbuild_cmd_line += ' -r ' + clustered_ref_db_basename + '.fasta'
        ovgraphbuild_cmd_line += ' -s ' + sam_filt_basename + '.sam'
        
        # Run ovgraphbuild
        sys.stdout.write('CMD: {0}\n\n'.format(ovgraphbuild_cmd_line))
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
        subprocess.call(contigsearch_cmd_line, shell=True)
        sys.stdout.write('\n')
    
    #######################
    # STEP 6: LCA Labelling
    quorum_int = int(args.quorum * 100)
    labelled_nodes_basename = contigsearch_basename + '.nodes_contracted'
    labelled_nodes_basename += '.metanode_contig_lca' + str(quorum_int) + 'pct'
    
    if 6 in steps_set:
        sys.stdout.write('## LCA Labelling (6):\n\n')
        
        #
        cmd_line = 'tail -n +2 ' + contigsearch_basename + '.contigs.csv'
        cmd_line += ' | sed "s/;/\\t/g" | sort -k2,2 | awk \'{print $2"\\t"$3"\\t"$1}\' > '
        cmd_line += contigsearch_basename + '.contigs.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'tail -n +2 ' + ovgraphbuild_basename + '.nodes.csv'
        cmd_line += ' | sed "s/;/\\t/g" | sort -k1,1 > '
        cmd_line += ovgraphbuild_basename + '.nodes.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        #
        complete_taxo_filename = contigsearch_basename + '.read_id_metanode_contig_ref_taxo.tab'
        
        cmd_line = 'join -11 -21 ' + ovgraphbuild_basename + '.nodes.tab '
        cmd_line += contigsearch_basename + '.contigs.tab '
        cmd_line += '| sort -k2,2 | awk \'{print $1"\\t"$2"\\t"$4"\\t"$5}\' | join -12 -21 - '
        cmd_line += read_ref_taxo_basename + '.tab | sed "s/ /\\t/g" > '
        cmd_line += complete_taxo_filename
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + complete_taxo_filename + ' | sort -k3,3 | '
        cmd_line += compute_lca_bin + ' -t 6 -f 3 -g 1 -m ' + str(args.quorum) + ' > '
        cmd_line += contigsearch_basename + '.metanode_lca' + str(quorum_int) + 'pct.tab'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + complete_taxo_filename + ' | sort -k4,4 | '
        cmd_line += compute_lca_bin + ' -t 6 -f 4 -g 1 -m ' + str(args.quorum) + ' > '
        cmd_line += contigsearch_basename + '.contig_lca' + str(quorum_int) + 'pct.tab'
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        #
        cmd_line = 'cat ' + contigsearch_basename + '.nodes_contracted.csv'
        cmd_line += ' | tail -n +2 | sed "s/;/\\t/g" | sort -k1,1 | join -11 -21 - '
        cmd_line += contigsearch_basename + '.metanode_lca' + str(quorum_int) + 'pct.tab'
        cmd_line += ' | sort -k4,4 | join -14 -21 - '
        cmd_line += contigsearch_basename + '.contig_lca' + str(quorum_int) + 'pct.tab'
        cmd_line += ' | sort -k4,4 | join -a1 -e"NULL" -o "1.2,1.3,0,1.1,2.2,1.5,1.6" -14 -21 - '
        cmd_line += '16sp.taxo.tab' + ' | sort -k3,3 | ' ## Filename in hard !!!!!!!
        cmd_line += 'awk \'BEGIN{print "Id Size Specie ContigId TrueTaxo NodeLCA ContigLCA"}{print $0}\''
        cmd_line += ' | sed "s/;/,/g" | sed "s/ /;/g" > '
        cmd_line += labelled_nodes_basename + '.csv'
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        # LCA Stats (only with a known test dataset)
        stats_filename = labelled_nodes_basename + '.stats'
        
        cmd_line = 'echo "## STATS LCA (NODE LEVEL):\n" > ' + stats_filename
        subprocess.call(cmd_line, shell=True)
        
        cmd_line = compute_stats_lca_bin + ' --header -l 6 -t 5 -i '
        cmd_line += labelled_nodes_basename + '.csv >> '
        cmd_line += stats_filename
        
        sys.stdout.write('CMD: {0}\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        
        cmd_line = 'echo "\n## STATS LCA (CONTIG LEVEL):\n" >> ' + stats_filename
        subprocess.call(cmd_line, shell=True)
        
        cmd_line = compute_stats_lca_bin + ' --header -l 7 -t 5 -i '
        cmd_line += labelled_nodes_basename + '.csv >> '
        cmd_line += stats_filename
        
        sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))
        subprocess.call(cmd_line, shell=True)
        sys.stdout.write('\n')
    
    ##########################
    # STEP 7: Contigs Assembly
    
    
    
    
    
    
    exit(0)


















