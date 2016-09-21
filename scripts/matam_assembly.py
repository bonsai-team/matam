#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import time
import logging


# Create logger
logger = logging.getLogger(__name__)

# Get program filename
program_filename = os.path.basename(sys.argv[0])

# Get MATAM root dir absolute path
matam_assembly_bin = os.path.realpath(sys.argv[0])
matam_bin_dir = os.path.dirname(matam_assembly_bin)
matam_root_dir = os.path.dirname(matam_bin_dir)
matam_db_dir = os.path.join(matam_root_dir, 'db')

# Set default ref db
default_ref_db = os.path.join(matam_db_dir, 'SILVA_123_SSURef_rdNs_NR95')

# Get all dependencies bin
matam_script_dir = os.path.join(matam_root_dir, 'scripts')
clean_name_bin = os.path.join(matam_script_dir, 'fasta_clean_name.py')
filter_score_bin = os.path.join(matam_script_dir, 'filter_score_multialign.py')
compute_lca_bin = os.path.join(matam_script_dir, 'compute_lca_from_tab.py')
compute_compressed_graph_stats_bin = os.path.join(matam_script_dir, 'compute_compressed_graph_stats.py')
sga_assemble_bin = os.path.join(matam_script_dir, 'sga_assemble.py')
remove_redundant_bin = os.path.join(matam_script_dir, 'remove_redundant_sequences.py')
fastq_name_filter_bin = os.path.join(matam_script_dir, 'fastq_name_filter.py')
evaluate_assembly_bin = os.path.join(matam_script_dir, 'evaluate_assembly.py')
get_best_matches_bin = os.path.join(matam_script_dir, 'get_best_matches_from_blast.py')
gener_scaff_blast_bin = os.path.join(matam_script_dir, 'generate_scaffolding_blast.py')
filter_sam_blast_bin = os.path.join(matam_script_dir, 'filter_sam_based_on_blast.py')
compute_contigs_compatibility_bin = os.path.join(matam_script_dir, 'compute_contigs_compatibility.py')
scaffold_contigs_bin = os.path.join(matam_script_dir, 'scaffold_contigs.py')
fasta_length_filter_bin = os.path.join(matam_script_dir, 'fasta_length_filter.py')

sortmerna_bin_dir = os.path.join(matam_root_dir, 'sortmerna')
sortmerna_bin = os.path.join(sortmerna_bin_dir, 'sortmerna')
indexdb_bin = os.path.join(sortmerna_bin_dir, 'indexdb_rna')

ovgraphbuild_bin_dir = os.path.join(matam_root_dir, 'ovgraphbuild', 'bin')
ovgraphbuild_bin = os.path.join(ovgraphbuild_bin_dir, 'ovgraphbuild')

componentsearch_bin_dir = os.path.join(matam_root_dir, 'componentsearch')
componentsearch_jar = os.path.join(componentsearch_bin_dir, 'ComponentSearch.jar')

assembler_bin_dir = os.path.join(matam_root_dir, 'sga', 'src', 'SGA')
assembler_name = 'sga'
assembler_bin = os.path.join(assembler_bin_dir, assembler_name)

# Define a null file handle
FNULL = open(os.devnull, 'w')


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


def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
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
    parser = DefaultHelpParser(description='MATAM assembly',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))

    # Main parameters
    group_main = parser.add_argument_group('Main parameters')
    # -i / --input_fastx
    group_main.add_argument('-i', '--input_fastx',
                            action = 'store',
                            metavar = 'FASTX',
                            type = str,
                            required = True,
                            help = 'Input reads file (fasta or fastq format)')
    # -d / --ref_db
    group_main.add_argument('-d', '--ref_db',
                            action = 'store',
                            metavar = 'DBPATH',
                            type = str,
                            help = 'MATAM ref db. '
                                   'Default is $MATAM_DIR/db/SILVA_123_SSURef_rdNs_NR95')
    # -o / --out_dir
    group_main.add_argument('-o', '--out_dir',
                            action = 'store',
                            metavar = 'OUTDIR',
                            type = str,
                            help = 'Output directory.'
                                   'Default will be "matam_assembly"')
    # -v / --verbose
    group_main.add_argument('-v', '--verbose',
                            action = 'store_true',
                            help = 'Increase verbosity')

    # Performance parameters
    group_perf = parser.add_argument_group('Performance')
    # --cpu
    group_perf.add_argument('--cpu',
                            action = 'store',
                            metavar = 'CPU',
                            type = int,
                            default = 1,
                            help = 'Max number of CPU to use. '
                                   'Default is %(default)s cpu')
    # --max_memory
    group_perf.add_argument('--max_memory',
                            action = 'store',
                            metavar = 'MAXMEM',
                            type = int,
                            default = 10000,
                            help = 'Maximum memory to use (in MBi). '
                                   'Default is %(default)s MBi')

    # Reads mapping parameters
    group_mapping = parser.add_argument_group('Read mapping')
    # --best
    group_mapping.add_argument('--best',
                               action = 'store',
                               metavar = 'INT',
                               type = int,
                               default = 10,
                               help = 'Get up to --best good alignments per read. '
                                      'Default is %(default)s')
    # --min_lis
    group_mapping.add_argument('--min_lis',
                               action = 'store',
                               metavar = 'INT',
                               type = int,
                               default = 10,
                               help = argparse.SUPPRESS)
    # --evalue
    group_mapping.add_argument('--evalue',
                               action = 'store',
                               metavar = 'REAL',
                               type = float,
                               default = 1e-5,
                               help = 'Max e-value to keep an alignment for. '
                                      'Default is %(default)s')

    # Alignment filtering parameters
    group_filt = parser.add_argument_group('Alignment filtering')
    # --score_threshold
    group_filt.add_argument('--score_threshold',
                            action = 'store',
                            metavar = 'REAL',
                            type = float,
                            default = 0.9,
                            help = 'Score threshold (real between 0 and 1). '
                                   'Default is %(default)s')
    # --straight_mode
    group_filt.add_argument('--straight_mode',
                            action = 'store_true',
                            help = 'Use straight mode filtering. '
                                   'Default is geometric mode')

    # Overlap-graph building parameters
    group_ovg = parser.add_argument_group('Overlap-graph building')
    # --min_identity
    group_ovg.add_argument('--min_identity',
                           action = 'store',
                           metavar = 'REAL',
                           type = float,
                           default = 1.0,
                           help = 'Minimum identity of an overlap between 2 reads. '
                                  'Default is %(default)s')
    # --min_overlap_length
    group_ovg.add_argument('--min_overlap_length',
                           action = 'store',
                           metavar = 'INT',
                           type = int,
                           default = 50,
                           help = 'Minimum length of an overlap. '
                                  'Default is %(default)s')

    # Graph compaction & Components identification
    group_gcomp = parser.add_argument_group('Graph compaction & Components identification')
    # -N / --min_read_node
    group_gcomp.add_argument('-N', '--min_read_node',
                             action = 'store',
                             metavar = 'INT',
                             type = int,
                             default = 2,
                             help = 'Minimum number of read to keep a node. '
                                    'Default is %(default)s')
    # -E / --min_overlap_edge
    group_gcomp.add_argument('-E', '--min_overlap_edge',
                             action = 'store',
                             metavar = 'INT',
                             type = int,
                             default = 10,
                             help = 'Minimum number of overlap to keep an edge. '
                                    'Default is %(default)s')

    # LCA labelling
    group_lca = parser.add_argument_group('LCA labelling')
    # --quorum
    group_lca.add_argument('--quorum',
                           action = 'store',
                           metavar = 'FLOAT',
                           type = float,
                           default = 0.51,
                           help = 'Quorum for LCA computing. Has to be between 0.51 and 1. '
                                  'Default is %(default)s')

    # Computing compressed graph stats

    # Contigs assembly

    # Scaffolding

    # Visualization

    # Advanced parameters
    group_adv = parser.add_argument_group('Advanced parameters')
    # --keep_tmp
    group_adv.add_argument('--keep_tmp',
                            action = 'store_true',
                            help = 'Do not remove tmp files')
    # --true_references
    # Fasta sequences of the known true references
    group_adv.add_argument('--true_references',
                           action = 'store',
                           type = str,
                           help = argparse.SUPPRESS)
    # --true_ref_taxo
    # Taxonomies of the true ref (represented by a 3 character id)
    # This is only needed when using a simulated dataset
    group_adv.add_argument('--true_ref_taxo',
                           action = 'store',
                           type = str,
                           help = argparse.SUPPRESS)
    # --debug
    group_adv.add_argument('--debug',
                            action = 'store_true',
                            help = 'Output debug infos')

    #
    args = parser.parse_args()

    # Arguments checking
    if args.score_threshold < 0 or args.score_threshold > 1:
        parser.print_help()
        raise Exception("score threshold not in range [0,1]")

    if args.min_identity < 0 or args.min_identity > 1:
        parser.print_help()
        raise Exception("min identity not in range [0,1]")

    if args.quorum < 0 or args.quorum > 1:
        parser.print_help()
        raise Exception("quorum not in range [0.51,1]")

    # Set debug parameters
    if args.debug:
        args.verbose = True
        args.keep_tmp = True

    # Set default ref db
    if not args.ref_db:
        args.ref_db = default_ref_db

    # Set default output dir
    if not args.out_dir:
        # args.out_dir = 'matam.{0}'.format(os.getpid())
        args.out_dir = 'matam_assembly'

    # Get absolute path for all arguments
    args.input_fastx = os.path.abspath(args.input_fastx)
    args.ref_db = os.path.abspath(args.ref_db)
    args.out_dir = os.path.abspath(args.out_dir)
    if args.true_references:
        args.true_references = os.path.abspath(args.true_references)
    if args.true_ref_taxo:
        args.true_ref_taxo = os.path.abspath(args.true_ref_taxo)

    #
    return args


def print_intro(args):
    """
    Print the introduction
    """

    sys.stderr.write("""
#################################
         MATAM assembly
#################################\n\n""")

    # Retrieve complete cmd line
    cmd_line = '{binpath} '.format(binpath=matam_assembly_bin)

    # Verbose
    if args.verbose:
        cmd_line += '--verbose '

    # Advanced parameters
    if args.debug:
        cmd_line += '--debug '

    if args.keep_tmp:
        cmd_line += '--keep_tmp '

    # Performance
    cmd_line += """--cpu {cpu} --max_memory {memory} \
""".format(cpu=args.cpu,
           memory=args.max_memory)

    # Read mapping
    cmd_line += '--best {0} '.format(args.best)
    cmd_line += '--evalue {0:.2e} '.format(args.evalue)

    # Alignment filtering
    cmd_line += '--score_threshold {0:.2f} '.format(args.score_threshold)
    if args.straight_mode:
        cmd_line += '--straight_mode '

    # Overlap-graph building
    cmd_line += '--min_identity {0:.2f} '.format(args.min_identity)
    cmd_line += '--min_overlap_length {0} '.format(args.min_overlap_length)

    # Graph compaction & Components identification
    cmd_line += '--min_read_node {0} '.format(args.min_read_node)
    cmd_line += '--min_overlap_edge {0} '.format(args.min_overlap_edge)

    # LCA labelling
    cmd_line += '--quorum {0:.2f} '.format(args.quorum)

    # Computing compressed graph stats

    # Contigs assembly

    # Scaffolding

    # Visualization

    # Main parameters
    cmd_line += '--out_dir {0}'.format(args.out_dir)
    cmd_line += '--ref_db {0} '.format(args.ref_db)
    cmd_line += '--input_fastx {0}'.format(args.input_fastx)

    # Print cmd line
    sys.stderr.write('CMD: {0}\n\n'.format(cmd_line))

    return 0


def rm_files(filepath_list):
    """
    Try to delete all files in filepath_list
    """
    for filepath in filepath_list:
        try:
            logger.debug('rm {0}'.format(filepath))
            os.remove(filepath)
        except:
            pass


if __name__ == '__main__':

    # Set global t0
    global_t0_wall = time.time()

    # Arguments parsing
    args = parse_arguments()

    # Print intro infos
    print_intro(args)

    # Init error code
    error_code = 0

    # Set logging
    # create console handler
    ch = logging.StreamHandler()
    #
    if args.debug:
        logger.setLevel(logging.DEBUG)
        # create formatter for debug level
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    else:
        if args.verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.WARNING)
        # create default formatter
        formatter = logging.Formatter('%(levelname)s - %(message)s')
    # add the formatter to the console handler
    ch.setFormatter(formatter)
    # add the handler to logger
    logger.addHandler(ch)

    # Init list of tmp files to delete at the end
    to_rm_filepath_list = list()

    ##############################################
    # Set all files and directories names + paths

    input_fastx_filepath = args.input_fastx
    input_fastx_filename = os.path.basename(input_fastx_filepath)
    input_fastx_basename, input_fastx_extension = os.path.splitext(input_fastx_filename)

    ref_db_basepath = args.ref_db
    ref_db_dir, ref_db_basename = os.path.split(ref_db_basepath)

    try:
        if not os.path.exists(args.out_dir):
            logger.debug('mkdir {0}'.format(args.out_dir))
            os.makedirs(args.out_dir)
    except OSError:
        logger.exception('Could not create output directory {0}'.format(args.out_dir))
        raise

    complete_ref_db_basename = ref_db_basename + '.complete'
    complete_ref_db_basepath = os.path.join(ref_db_dir, complete_ref_db_basename)
    complete_ref_db_filename = complete_ref_db_basename + '.fasta'
    complete_ref_db_filepath = os.path.join(ref_db_dir, complete_ref_db_filename)

    complete_ref_db_taxo_filename = complete_ref_db_basename + '.taxo.tab'
    complete_ref_db_taxo_filepath = os.path.join(ref_db_dir, complete_ref_db_taxo_filename)

    clustered_ref_db_basename = ref_db_basename + '.clustered'
    clustered_ref_db_basepath = os.path.join(ref_db_dir, clustered_ref_db_basename)
    clustered_ref_db_filename = clustered_ref_db_basename + '.fasta'
    clustered_ref_db_filepath = os.path.join(ref_db_dir, clustered_ref_db_filename)

    # Read mapping
    sortme_output_basename = input_fastx_basename
    sortme_output_basename += '.sortmerna_vs_' + ref_db_basename
    sortme_output_basename += '_b' + str(args.best) + '_m' + str(args.min_lis)
    sortme_output_basepath = os.path.join(args.out_dir, sortme_output_basename)

    sortme_output_fastx_filepath = sortme_output_basepath + input_fastx_extension
    sortme_output_sam_filepath = sortme_output_basepath + '.sam'

    # Alignment filtering
    score_threshold_int = int(args.score_threshold * 100)

    sam_filt_basename = sortme_output_basename + '.scr_filt_'
    if args.straight_mode:
        sam_filt_basename += 'str_'
    else:
        sam_filt_basename += 'geo_'
    sam_filt_basename += str(score_threshold_int) + 'pct'
    sam_filt_filename = sam_filt_basename + '.sam'
    sam_filt_filepath = os.path.join(args.out_dir, sam_filt_filename)

    # Overlap-graph building
    min_identity_int = int(args.min_identity * 100)

    ovgraphbuild_basename = sam_filt_basename + '.ovgb_i' + str(min_identity_int)
    ovgraphbuild_basename += '_o' + str(args.min_overlap_length)
    ovgraphbuild_basepath = os.path.join(args.out_dir, ovgraphbuild_basename)

    ovgraphbuild_nodes_csv_filepath = ovgraphbuild_basepath + '.nodes.csv'
    ovgraphbuild_edges_csv_filepath = ovgraphbuild_basepath + '.edges.csv'

    # Graph compaction & Components identification
    componentsearch_basename = ovgraphbuild_basename + '.ctgs'
    componentsearch_basename += '_N' + str(args.min_read_node)
    componentsearch_basename += '_E' + str(args.min_overlap_edge)
    componentsearch_basepath = os.path.join(args.out_dir, componentsearch_basename)

    contracted_nodes_basepath = componentsearch_basepath + '.nodes_contracted'
    contracted_nodes_filepath = contracted_nodes_basepath + '.csv'
    contracted_edges_filepath = componentsearch_basepath + '.edges_contracted.csv'

    # LCA labelling
    read_id_metanode_component_filepath = componentsearch_basepath + '.read_id_metanode_component.tab'
    complete_taxo_filepath = componentsearch_basepath + '.read_metanode_component_taxo.tab'

    quorum_int = int(args.quorum * 100)

    labelled_nodes_basename = componentsearch_basename + '.nodes_contracted'
    labelled_nodes_basename += '.component_lca' + str(quorum_int) + 'pct'

    components_lca_filename = componentsearch_basename + '.component_lca' + str(quorum_int) + 'pct.tab'
    components_lca_filepath = os.path.join(args.out_dir, components_lca_filename)

    # Computing compressed graph stats
    stats_filename = componentsearch_basename + '.graph.stats'
    stats_filepath = os.path.join(args.out_dir, stats_filename)

    # Contigs assembly
    component_read_filepath = componentsearch_basepath + '.component_read.tab'

    contigs_assembly_wkdir = componentsearch_basepath + '.' + assembler_name
    try:
        if not os.path.exists(contigs_assembly_wkdir):
            logger.debug('mkdir {0}'.format(contigs_assembly_wkdir))
            os.makedirs(contigs_assembly_wkdir)
    except OSError:
        logger.exception('Assembly directory {0} cannot be created'.format(contigs_assembly_wkdir))
        raise

    contigs_basename = componentsearch_basename + '.'
    contigs_basename += assembler_name + '_by_component'
    contigs_filename = contigs_basename + '.fasta'
    contigs_filepath = os.path.join(args.out_dir, contigs_filename)

    contigs_assembly_log_filename = contigs_basename + '.log'
    contigs_assembly_log_filepath = os.path.join(args.out_dir, contigs_assembly_log_filename)
    # Remove log file from previous assemblies
    # because we will only append to this file
    if os.path.exists(contigs_assembly_log_filepath):
        os.remove(contigs_assembly_log_filepath)

    contigs_symlink_basename = 'contigs'
    contigs_symlink_filename = contigs_symlink_basename + '.fasta'
    contigs_symlink_filepath = os.path.join(args.out_dir, contigs_symlink_filename)

    # Scaffolding
    scaff_sortme_output_basename = contigs_symlink_basename
    scaff_sortme_output_basename += '.sortmerna_vs_complete_' + ref_db_basename
    scaff_sortme_output_basename += '_num_align_0'
    scaff_sortme_output_basepath = os.path.join(args.out_dir, scaff_sortme_output_basename)

    scaff_sortme_output_blast_filename = scaff_sortme_output_basename + '.blast'
    scaff_sortme_output_blast_filepath = os.path.join(args.out_dir, scaff_sortme_output_blast_filename)

    scaff_sortme_output_sam_filename = scaff_sortme_output_basename + '.sam'
    scaff_sortme_output_sam_filepath = os.path.join(args.out_dir, scaff_sortme_output_sam_filename)

    best_only_blast_basename = scaff_sortme_output_blast_filename + '.best_only'
    best_only_blast_filename = best_only_blast_basename + '.tab'
    best_only_blast_filepath = os.path.join(args.out_dir, best_only_blast_filename)

    selected_best_only_blast_basename = best_only_blast_basename + '.selected'
    selected_best_only_blast_filename = selected_best_only_blast_basename + '.tab'
    selected_best_only_blast_filepath = os.path.join(args.out_dir, selected_best_only_blast_filename)

    selected_sam_filename = selected_best_only_blast_basename + '.sam'
    selected_sam_filepath = os.path.join(args.out_dir, selected_sam_filename)

    binned_sam_basename = selected_best_only_blast_basename + '.binned'
    binned_sam_filename = binned_sam_basename + '.sam'
    binned_sam_filepath = os.path.join(args.out_dir, binned_sam_filename)

    bam_filename = binned_sam_basename + '.bam'
    bam_filepath = os.path.join(args.out_dir, bam_filename)

    sorted_bam_basename = binned_sam_basename + '.sorted'
    sorted_bam_basepath = os.path.join(args.out_dir, sorted_bam_basename)
    sorted_bam_filename = sorted_bam_basename + '.bam'
    sorted_bam_filepath = os.path.join(args.out_dir, sorted_bam_filename)

    mpileup_filename = sorted_bam_basename + '.mpileup'
    mpileup_filepath = os.path.join(args.out_dir, mpileup_filename)

    scaffolds_basename = sorted_bam_basename + '.scaffolds'
    scaffolds_filename = scaffolds_basename + '.fasta'
    scaffolds_filepath = os.path.join(args.out_dir, scaffolds_filename)

    scaffolds_symlink_basename = 'scaffolds'
    scaffolds_symlink_filename = scaffolds_symlink_basename + '.fasta'
    scaffolds_symlink_filepath = os.path.join(args.out_dir, scaffolds_symlink_filename)

    scaffolds_NR_basename = scaffolds_symlink_basename + '.NR'
    scaffolds_NR_filename = scaffolds_NR_basename + '.fasta'
    scaffolds_NR_filepath = os.path.join(args.out_dir, scaffolds_NR_filename)

    large_scaffolds_NR_basename = scaffolds_NR_basename + '.min_' + str(500) + 'bp'
    large_scaffolds_NR_filename = large_scaffolds_NR_basename + '.fasta'
    large_scaffolds_NR_filepath = os.path.join(args.out_dir, large_scaffolds_NR_filename)

    ###############################
    # Reads mapping against ref db

    logger.info('Reads mapping against ref db')

    cmd_line = sortmerna_bin + ' --ref ' + clustered_ref_db_filepath
    cmd_line += ',' + clustered_ref_db_basepath + ' --reads '
    cmd_line += input_fastx_filepath + ' --aligned ' + sortme_output_basepath
    cmd_line += ' --fastx --sam --blast "1" --log --best '
    cmd_line += str(args.best) + ' --min_lis ' + str(args.min_lis)
    cmd_line += ' -e {0:.2e}'.format(args.evalue)
    cmd_line += ' -a ' + str(args.cpu)
    if args.verbose:
        cmd_line += ' -v '

    # Set t0
    t0_wall = time.time()

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)
    if args.verbose:
        sys.stderr.write('\n')

    # Output running time
    logger.debug('Reads mapping terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    #############################
    # Alignment filtering

    logger.info('Alignment filtering')

    cmd_line = 'cat ' + sortme_output_sam_filepath
    cmd_line += ' | grep -v "^@" | sort -k 1,1V -k 12,12nr'
    cmd_line += ' | ' + filter_score_bin + ' -t ' + str(args.score_threshold)
    if not args.straight_mode:
        cmd_line += ' --geometric'
    cmd_line += ' > ' + sam_filt_filepath

    # Set t0
    t0_wall = time.time()

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Output running time
    logger.debug('Alignment filtering terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Tag tmp files for removal
    to_rm_filepath_list.append(sortme_output_sam_filepath)
    to_rm_filepath_list.append(sortme_output_basepath + '.log')
    to_rm_filepath_list.append(sortme_output_basepath + '.blast')

    #########################
    # Overlap-graph building

    logger.info('Overlap-graph building')

    cmd_line = ovgraphbuild_bin
    cmd_line += ' -i ' + str(args.min_identity)
    cmd_line += ' -m ' + str(args.min_overlap_length)
    cmd_line += ' --csv --output_basename '
    cmd_line += ovgraphbuild_basepath
    cmd_line += ' -r ' + clustered_ref_db_filepath
    cmd_line += ' -s ' + sam_filt_filepath
    if args.verbose:
        cmd_line += ' -v'
    if args.debug:
        cmd_line += ' --debug'

    # Set t0
    t0_wall = time.time()

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Output running time
    logger.debug('Overlap-graph building terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Tag tmp files for removal
    to_rm_filepath_list.append(sam_filt_filepath)

    ###############################################
    # Graph compaction & Components identification

    logger.info('Graph compaction & Components identification')

    cmd_line = 'java -Xmx' + str(args.max_memory) + 'M -cp "'
    cmd_line += componentsearch_jar + '" main.Main'
    cmd_line += ' -N ' + str(args.min_read_node)
    cmd_line += ' -E ' + str(args.min_overlap_edge)
    cmd_line += ' -b ' + componentsearch_basepath
    cmd_line += ' -n ' + ovgraphbuild_nodes_csv_filepath
    cmd_line += ' -e ' + ovgraphbuild_edges_csv_filepath

    # Set t0
    t0_wall = time.time()

    logger.debug('CMD: {0}'.format(cmd_line))
    if args.verbose:
        error_code += subprocess.call(cmd_line, shell=True)
    else:
        # Needed because ComponentSearch doesnt have a verbose option
        # and output everything to stderr
        error_code += subprocess.call(cmd_line, shell=True, stderr=FNULL)

    # Output running time
    logger.debug('Graph compaction & Components identification terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Tag tmp files for removal
    to_rm_filepath_list.append(ovgraphbuild_nodes_csv_filepath)
    to_rm_filepath_list.append(ovgraphbuild_edges_csv_filepath)

    ################
    # LCA labelling

    logger.info('LCA labelling')

    # Set t0
    t0_wall = time.time()

    # Note: some of the manipulations here are needed because ComponentSearch
    # works with read ids rather than read names. The goal of this part is to
    # regroup all the infos in one file (component, metanode, read, ref taxo)
    # in order to compute LCA at metanode or component level

    # Convert CSV component file to TAB format and sort by read id
    cmd_line = 'tail -n +2 ' + componentsearch_basepath + '.components.csv'
    cmd_line += ' | sed "s/;/\\t/g" | sort -k2,2 > '
    cmd_line += componentsearch_basepath + '.components.tab'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Convert CSV node file to TAB format and sort by read id
    cmd_line = 'tail -n +2 ' + ovgraphbuild_nodes_csv_filepath
    cmd_line += ' | sed "s/;/\\t/g" | sort -k1,1 > '
    cmd_line += ovgraphbuild_basepath + '.nodes.tab'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Join component and node files on read id, and sort by read name
    cmd_line = 'join -a1 -e"NULL" -o "1.2,0,2.3,2.1" -11 -22 '
    cmd_line += ovgraphbuild_basepath + '.nodes.tab '
    cmd_line += componentsearch_basepath + '.components.tab '
    cmd_line += '| sed "s/ /\\t/g" | sort -k1,1  > '
    cmd_line += read_id_metanode_component_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Join sam file with ref taxo on ref name, and join it to the
    # component-node file on read name.
    # Output is (read, metanode, component, taxo) file
    cmd_line = 'cat ' + sam_filt_filepath + ' | cut -f1,3 | sort -k2,2'
    cmd_line += ' | join -12 -21 - ' + complete_ref_db_taxo_filepath
    cmd_line += ' | sort -k2,2 | awk "{print \$2,\$3}" | sed "s/ /\\t/g" '
    cmd_line += ' | join -11 -21 ' + read_id_metanode_component_filepath
    cmd_line += ' - | sed "s/ /\\t/g" | cut -f2-5 > '
    cmd_line += complete_taxo_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Compute LCA at component level using quorum threshold
    cmd_line = 'cat ' + complete_taxo_filepath + ' | sort -k3,3 -k1,1 | '
    cmd_line += compute_lca_bin + ' -t 4 -f 3 -g 1 -m ' + str(args.quorum)
    cmd_line += ' -o ' + components_lca_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Output running time
    logger.debug('LCA labelling terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Tag tmp files for removal
    to_rm_filepath_list.append(componentsearch_basepath + '.components.tab')
    to_rm_filepath_list.append(ovgraphbuild_basepath + '.nodes.tab')
    to_rm_filepath_list.append(complete_taxo_filepath)

    ###################################
    # Computing compressed graph stats

    logger.info('Computing compressed graph stats')

    cmd_line = compute_compressed_graph_stats_bin + ' --nodes_contracted '
    cmd_line += contracted_nodes_filepath + ' --edges_contracted '
    cmd_line += contracted_edges_filepath + ' --components_lca '
    cmd_line += components_lca_filepath
    if args.true_ref_taxo:
        cmd_line += ' --species_taxo ' + args.true_ref_taxo
        cmd_line += ' --read_node_component ' + read_id_metanode_component_filepath
    cmd_line += ' -o ' + stats_filepath

    # Set t0
    t0_wall = time.time()

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Output running time
    logger.debug('Computing compressed graph stats terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    ###################
    # Contigs assembly

    logger.info('Starting contigs assembly')

    # Set t0
    t0_wall = time.time()

    # Generate the read_id-->component_id file
    cmd_line = 'cat ' + read_id_metanode_component_filepath
    cmd_line += ' | awk \' {print $4"\\t"$1}\' | sort -k1,1n > '
    cmd_line += component_read_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Count the number of components to assemble
    cmd_line = 'cat ' + component_read_filepath
    cmd_line += ' | cut -f1 | grep -v "NULL" | sort | uniq | wc -l'

    logger.debug('CMD: {0}'.format(cmd_line))
    components_num = int(subprocess.check_output(cmd_line, shell=True))

    # TO DO, one day, maybe:
    # Convert input fastq to tab, join it to the component-read file
    # on the read name. Then generate a fastq file for each component.
    # This would take more space but be faster and less error prone
    # than using name filtering on the complete fastq each time

    # Reading components LCA and storing them in a dict
    logger.debug('Reading components LCA assignment from {0}'.format(components_lca_filepath))

    component_lca_dict = dict()
    with open(components_lca_filepath, 'r') as component_lca_fh:
        component_lca_dict = {t[0]: t[1] for t in (l.split() for l in component_lca_fh) if len(t) == 2}

    # Open output contigs file
    contigs_fh = open(contigs_filepath, 'w')

    # Begin reading read-->component file (sorted by component)
    contig_count = 0
    component_count = 0

    with open(component_read_filepath, 'r') as component_read_fh:
        # We will be working in the assembly directory because assembly tools
        # generate lots of files we dont want in our matam assembly directory
        os.chdir(contigs_assembly_wkdir)

        #
        logger.info('Assembling component #')

        # Assembling all reads of every component, one at a time
        for component_tab_list in read_tab_file_handle_sorted(component_read_fh, 0):

            component_id = component_tab_list[0][0]

            # To prevent assembly of singleton reads
            if component_id == 'NULL':
                continue

            component_count += 1

            # Starting component assembly
            if args.verbose:
                sys.stderr.write('\r{0} / {1}'.format(component_count, components_num))
                sys.stderr.flush()

            # Cleaning previous assembly
            if os.path.exists('contigs.fa'):
                os.remove('contigs.fa')
            if os.path.exists('tmp'):
                subprocess.call('rm -rf tmp', shell=True)

            # Write the component reads ids
            with open('reads_single_component.ids', 'w') as wfh:
                for tab in component_tab_list:
                    read_id = tab[1]
                    wfh.write('{0}\n'.format(read_id))

            # Generate a fastq file with this component reads
            cmd_line = fastq_name_filter_bin + ' -f reads_single_component.ids'
            cmd_line += ' -i ' + sortme_output_fastx_filepath + ' -o reads_single_component.fq'

            error_code += subprocess.call(cmd_line, shell=True)

            # Assemble those reads with SGA
            cmd_line = 'echo "component #' + component_id + '" >> '
            cmd_line += contigs_assembly_log_filepath + ' && '
            cmd_line += sga_assemble_bin + ' -i reads_single_component.fq'
            cmd_line += ' -o contigs.fa --sga_bin ' + assembler_bin
            cmd_line += ' --cpu ' + str(args.cpu)
            cmd_line += ' --tmp_dir tmp'
            cmd_line += ' >> ' + contigs_assembly_log_filepath + ' 2>&1'

            error_code += subprocess.call(cmd_line, shell=True)

            # Concatenate the component contigs in the output contigs file
            component_lca = 'NULL'
            if component_id in component_lca_dict:
                component_lca = component_lca_dict[component_id]
            with open('contigs.fa', 'r') as sga_contigs_fh:
                for header, seq in read_fasta_file_handle(sga_contigs_fh):
                    if len(seq):
                        contig_count += 1
                        contigs_fh.write('>{0} component={1} '.format(contig_count, component_id))
                        contigs_fh.write('lca={0}\n{1}\n'.format(component_lca, format_seq(seq)))

        #
        if args.verbose:
            sys.stderr.write('\n')

        # Return to matam assembly directory
        os.chdir(args.out_dir)

    # Close assembly contigs file
    contigs_fh.close()

    # Create symbolic link
    if os.path.exists(contigs_symlink_filepath):
        os.remove(contigs_symlink_filepath)
    os.symlink(os.path.basename(contigs_filepath), contigs_symlink_filepath)

    # TO DO, if it gets better results:
    # remove redundant sequences in contigs

    # Output running time
    logger.debug('Contigs assembly terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    ## Evaluate assembly if true ref are provided
    #if args.true_references:
        #cmd_line = evaluate_assembly_bin + ' -r ' + args.true_references
        #cmd_line += ' -i ' + contigs_symlink_filepath

        #logger.debug('CMD: {0}'.format(cmd_line))
        #error_code += subprocess.call(cmd_line, shell=True)

    # Tag tmp files for removal
    to_rm_filepath_list.append(read_id_metanode_component_filepath)

    # Delete assembly directory
    if not args.keep_tmp:
        subprocess.call('rm -rf {0}'.format(contigs_assembly_wkdir), shell=True)

    ##############
    # Scaffolding

    logger.info('Scaffolding')

    # Set t0
    t0_wall = time.time()

    # Contigs remapping on the complete database
    scaff_evalue = 1e-05

    cmd_line = sortmerna_bin
    cmd_line += ' --ref ' + complete_ref_db_filepath + ',' + complete_ref_db_basepath
    cmd_line += ' --reads ' + contigs_filepath
    cmd_line += ' --aligned ' + scaff_sortme_output_basepath
    cmd_line += ' --sam --blast "1"'
    cmd_line += ' --num_alignments 0 '
    #~ cmd_line += ' --best 0 --min_lis 10 '
    cmd_line += ' -e {0:.2e}'.format(scaff_evalue)
    cmd_line += ' -a ' + str(args.cpu)
    if args.verbose:
        cmd_line += ' -v '

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)
    if args.verbose:
        sys.stderr.write('\n')

    # Output running time
    logger.debug('[scaff] Contig mapping terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Set t0
    t0_wall = time.time()

    # Keep only quasi-equivalent best matches for each contig
    # -p 0.99 is used because of the tendency of SortMeRNA to soft-clip
    # a few nucleotides at 5’ and 3’ ends
    cmd_line = 'sort -k1,1V -k12,12nr ' + scaff_sortme_output_blast_filepath
    cmd_line += ' | ' + get_best_matches_bin + ' -p 0.99 -o ' + best_only_blast_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Select blast matches for scaffolding using a specific-first conserved-later approach
    cmd_line = gener_scaff_blast_bin + ' -i ' + best_only_blast_filepath
    cmd_line += ' -o ' + selected_best_only_blast_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Filter sam file based on blast scaffolding file
    cmd_line = filter_sam_blast_bin + ' -i ' + scaff_sortme_output_sam_filepath
    cmd_line += ' -b ' +  selected_best_only_blast_filepath
    cmd_line += ' | sort -k3,3 -k4,4n > ' + selected_sam_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Bin compatible contigs matching on the same reference
    cmd_line = compute_contigs_compatibility_bin
    cmd_line += ' -i ' + selected_sam_filepath
    cmd_line += ' | sort -k1,1n > ' + binned_sam_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Convert sam to bam
    cmd_line = 'samtools view -b -S ' + binned_sam_filepath
    cmd_line += ' -o ' + bam_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Sort bam
    cmd_line = 'samtools sort ' + bam_filepath + ' ' + sorted_bam_basepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Generate mpileup
    cmd_line = 'samtools mpileup -d 10000 ' + sorted_bam_filepath
    cmd_line += ' > ' + mpileup_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Scaffold contigs based on mpileup
    cmd_line = scaffold_contigs_bin + ' -i ' + mpileup_filepath
    cmd_line += ' -o ' + scaffolds_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Create symbolic link
    if os.path.exists(scaffolds_symlink_filepath):
        os.remove(scaffolds_symlink_filepath)
    os.symlink(os.path.basename(scaffolds_filepath), scaffolds_symlink_filepath)

    # Remove redundant scaffolds
    cmd_line = remove_redundant_bin + ' -i ' + scaffolds_symlink_filepath
    cmd_line += ' -o ' + scaffolds_NR_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Filter out small scaffolds
    cmd_line = fasta_length_filter_bin + ' -m ' + str(500)
    cmd_line += ' -i ' + scaffolds_NR_filepath
    cmd_line += ' -o ' + large_scaffolds_NR_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Output running time
    logger.debug('[scaff] Scaffolding terminated in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Evaluate assembly if true ref are provided
    if args.true_references:
        cmd_line = evaluate_assembly_bin + ' -r ' + args.true_references
        cmd_line += ' -i ' + scaffolds_symlink_filepath

        logger.debug('CMD: {0}'.format(cmd_line))
        error_code += subprocess.call(cmd_line, shell=True)

    # Tag tmp files for removal
    to_rm_filepath_list.append(scaff_sortme_output_blast_filepath)
    to_rm_filepath_list.append(scaff_sortme_output_sam_filepath)
    to_rm_filepath_list.append(best_only_blast_filepath)
    to_rm_filepath_list.append(selected_best_only_blast_filepath)
    to_rm_filepath_list.append(selected_sam_filepath)
    to_rm_filepath_list.append(binned_sam_filepath)
    to_rm_filepath_list.append(bam_filepath)
    to_rm_filepath_list.append(sorted_bam_filepath)
    to_rm_filepath_list.append(mpileup_filepath)

    ###############
    # Exit program

    exit_code = 0

    # Exit if everything went ok
    if not error_code:
        # Try to remove all tmp files
        # won't crash if it cannot
        if not args.keep_tmp:
            logger.info('Removing tmp files')
            rm_files(to_rm_filepath_list)
        #
        sys.stderr.write('\n{0} terminated with no error\n'.format(program_filename))
    # Deal with errors
    else:
        sys.stderr.write('\n{0} terminated with some errors. '.format(program_filename))
        if args.verbose:
            sys.stderr.write('Check the log for additional infos\n')
        else:
            sys.stderr.write('Rerun the program using --verbose or --debug option\n')
        exit_code = 1

    logger.info('Run terminated in {0:.4f} seconds wall time'.format(time.time() - global_t0_wall))

    exit(exit_code)
