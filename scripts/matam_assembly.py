#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import statistics
import subprocess
import time
import logging
import cProfile
from collections import defaultdict
import shutil

import runner
from compute_abundance import get_abundance_by_scaffold, complete_fasta_with_abundance, get_abundance_from_fasta
from rdp import run_rdp_classifier, filter_rdp_file
from krona import rdp_file_to_krona_text_file, make_krona_plot
from binary_utils import Binary
import components_assembly
from assembler_factory import AssemblerFactory

# Set LC_LANG to C for standard sort behaviour
os.environ["LC_ALL"] = "C"

# Create logger
logger = logging.getLogger()

logging_config = dict(
    version = 1,
    disable_existing_loggers = False,
    root = {
        'handlers': ['console'],
        'level': logging.INFO
    },
    formatters = {
        'default': {'format': '%(levelname)s - %(message)s'},
    },
    handlers = {
        'console': {'class': 'logging.StreamHandler',
              'formatter': 'default'},
    },


)

def update_logger_settings(logger_filepath, verbose, debug):
    """
    update logging level and logging format accordinlgy to verbose/debug flag.
    add a file logger for root/runner logger.
    """
    file_handler = logging.FileHandler(filename=logger_filepath, encoding='utf8')
    file_handler.formatter = logging.Formatter('%(levelname)s - %(message)s')
    logger.addHandler(file_handler)
    if debug:
        logger.setLevel(logging.DEBUG)
        # update formatter for debug level
        for handler in logger.handlers:
            handler.formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')


# Get program filename
program_filename = os.path.basename(sys.argv[0])

# Get MATAM root dir absolute path
matam_assembly_bin = os.path.realpath(sys.argv[0])
matam_bin_dir = os.path.dirname(matam_assembly_bin)
matam_root_dir = os.path.dirname(matam_bin_dir)
matam_db_dir = os.path.join(matam_root_dir, 'db')

# Set default ref db
default_ref_db = os.path.join(matam_db_dir, 'SILVA_128_SSURef_NR95')

# Get all dependencies bin
matam_script_dir = os.path.join(matam_root_dir, 'scripts')
clean_name_bin = os.path.join(matam_script_dir, 'fasta_clean_name.py')
sample_sam_cov_bin = Binary.assert_which('sample_sam_by_coverage.py')
filter_score_bin = os.path.join(matam_script_dir, 'filter_score_multialign.py')
compute_lca_bin = os.path.join(matam_script_dir, 'compute_lca_from_tab.py')
compute_compressed_graph_stats_bin = os.path.join(matam_script_dir, 'compute_compressed_graph_stats.py')
remove_redundant_bin = os.path.join(matam_script_dir, 'remove_redundant_sequences.py')
fastq_name_filter_bin = os.path.join(matam_script_dir, 'fastq_name_filter.py')
evaluate_assembly_bin = os.path.join(matam_script_dir, 'evaluate_assembly.py')
get_best_matches_bin = os.path.join(matam_script_dir, 'get_best_matches_from_blast.py')
gener_scaff_blast_bin = os.path.join(matam_script_dir, 'generate_scaffolding_blast.py')
filter_sam_blast_bin = os.path.join(matam_script_dir, 'filter_sam_based_on_blast.py')
compute_contigs_compatibility_bin = os.path.join(matam_script_dir, 'compute_contigs_compatibility.py')
scaffold_contigs_bin = os.path.join(matam_script_dir, 'scaffold_contigs.py')
fasta_length_filter_bin = os.path.join(matam_script_dir, 'fasta_length_filter.py')
sortmerna_bin = Binary.assert_which('sortmerna')
indexdb_bin = Binary.assert_which('indexdb_rna')
ovgraphbuild_bin = Binary.assert_which('ovgraphbuild')
componentsearch_bin = Binary.assert_which('componentsearch')
krona_bin = Binary.assert_which('ktImportText')

rdp_jar = Binary.which('classifier.jar')

# the rdp exe name is different between submodule installation and conda installation
if rdp_jar is not None:
    java = Binary.assert_which('java')
    rdp_exe = '{java} -Xmx1g -jar {jar}'.format(java=java, jar=rdp_jar)
else:
    rdp_exe = Binary.assert_which('classifier')


def force_symlink(target, link_name):
    if os.path.exists(link_name):
        os.remove(link_name)
    os.symlink(target, link_name)


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
    if header or seqlines:
        yield (header, ''.join(seqlines))
    # Close input file
    fasta_file_handle.close()


def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    # Variables initialization
    count = 0
    header = ''
    seq = ''
    qual = ''
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        count += 1
        if count % 4 == 1:
            if header:
                yield header, seq, qual
            header = line[1:].split()[0]
        elif count % 4 == 2:
            seq = line
        elif count % 4 == 0:
            qual = line
    # yield last fastq sequence
    yield header, seq, qual
    # Close input file
    fastq_file_handle.close()


class FastaStats():
    """
    """
    def __init__(self):
        self.seq_num = 0
        self.total_nt = int()
        self.seq_length_list = list()

    def add_sequence(self, seq):
        self.seq_length_list.append(len(seq))
        self.total_nt += len(seq)
        self.seq_num += 1

    def get_avg_length(self):
        if self.seq_length_list:
            return statistics.mean(self.seq_length_list)
        else:
            return -1

    def get_max_length(self):
        if self.seq_length_list:
            return max(self.seq_length_list)
        else:
            return -1

    def get_min_length(self):
        if self.seq_length_list:
            return min(self.seq_length_list)
        else:
            return -1


def compute_fasta_stats(fasta_filepath):
    """
    """
    fasta_stats = FastaStats()

    with open(fasta_filepath, 'r') as fh:
        for header, seq in read_fasta_file_handle(fh):
            fasta_stats.add_sequence(seq)

    return fasta_stats


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
        logger.fatal(message)
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
                                   'Default is $MATAM_DIR/db/SILVA_128_SSURef_NR95')
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
    # --coverage_threshold
    group_filt.add_argument('--coverage_threshold',
                            action = 'store',
                            metavar = 'INT',
                            type = int,
                            default = 0,
                            help = 'Ref coverage threshold. '
                                   'By default set to 0 to desactivate filtering')

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
                             default = 1,
                             help = 'Minimum number of read to keep a node. '
                                    'Default is %(default)s')
    # -E / --min_overlap_edge
    group_gcomp.add_argument('-E', '--min_overlap_edge',
                             action = 'store',
                             metavar = 'INT',
                             type = int,
                             default = 1,
                             help = 'Minimum number of overlap to keep an edge. '
                                    'Default is %(default)s')
    # --seed
    group_gcomp.add_argument('--seed',
                             action = 'store',
                             metavar = 'INT',
                             type = int,
                             help = 'Seed to initialize random generator. '
                                    'Default is picking seed from system time')
    # --optimize_components
    group_gcomp.add_argument('--optimize_components',
                             action = 'store_true',
                             help = argparse.SUPPRESS)

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
    group_contig = parser.add_argument_group('Contigs assembly')

    # -a/--assembler
    group_contig.add_argument('-a', '--assembler',
                        choices=[a.name() for a in AssemblerFactory.ASSEMBLER_ENGINES],
                        help="Select the assembler to be used. Default is %(default)s",
                        default="SGA")

    # --read_correction
    group_contig.add_argument('--read_correction',
                              action = 'store',
                              type = str,
                              choices = ['no', 'yes', 'auto'],
                              default = 'auto',
                              help = 'Set the assembler read correction step. '
                                     "'auto' means assemble the components with read correction enabled when the components coverage are sufficient (20X) and assemble the other components without read correction. "
                                     'Default is %(default)s')
    # --contig_coverage_threshold
    # this option is not shown to the user
    group_contig.add_argument( '--contig_coverage_threshold',
                              action = 'store',
                              type = int,
                              choices = [20, 50],
                              default = 20,
                              help = argparse.SUPPRESS)

    # Scaffolding
    group_scaff = parser.add_argument_group('Scaffolding')
    # --contigs_binning
    group_scaff.add_argument('--contigs_binning',
                             action = 'store_true',
                             help = 'Experimental. Perform contigs binning during scaffolding.')

    group_scaff.add_argument('--min_scaffold_length',
                             action = 'store',
                             type = int,
                             default = 500,
                             help = 'Filter out small scaffolds'
                             )

    # taxonomic assignment
    group_taxonomic_assign = parser.add_argument_group('Taxonomic assignment')

    # --perform_taxonomic_assignment
    group_taxonomic_assign.add_argument('--perform_taxonomic_assignment',
                           action = 'store_true',
                           help = 'Do the taxonomic assignment')

    # --training_model
    group_taxonomic_assign.add_argument('--training_model',
                           action = 'store',
                           type = str,
                           choices = ['16srrna', 'fungallsu', 'fungalits_warcup', 'fungalits_unite'],
                           default = '16srrna',
                           help = 'The training model used for taxonomic assignment. '
                                     'Default is %(default)s')


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
    # --resume_from
    group_adv.add_argument('--resume_from',
                           action = 'store',
                           metavar = 'STEP',
                           type = str,
                           choices = ['reads_mapping', 'alignments_filtering', 'overlap_graph_building',
                                      'graph_compaction', 'contigs_assembly', 'scaffolding',
                                      'abundance_calculation', 'taxonomic_assignment'],
                           help = 'Try to resume from given step. '
                                  'Steps are: %(choices)s')

    # --filter_only
    group_adv.add_argument('--filter_only',
                            action = 'store_true',
                            help = 'Perform only the first step of MATAM (i.e Reads mapping against ref db with sortmerna to filter the reads). '\
                                   'Relevant options for this step correspond to the "Read mapping" section.')

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

    # contig_coverage_threshold default value makes no sense when args.read_correction is not auto
    if args.read_correction != 'auto':
        args.contig_coverage_threshold = None

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

    logger.info('=== MATAM assembly ===')

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

    if args.resume_from:
        cmd_line += '--resume_from {} '.format(args.resume_from)

    if args.true_references:
        cmd_line += '--true_references {0} '.format(args.true_references)

    if args.true_ref_taxo:
        cmd_line += '--true_ref_taxo {0} '.format(args.true_ref_taxo)

    # Performance
    cmd_line += '--cpu {0} '.format(args.cpu)
    cmd_line += '--max_memory {0} '.format(args.max_memory)

    # Read mapping
    cmd_line += '--best {0} '.format(args.best)
    cmd_line += '--evalue {0:.2e} '.format(args.evalue)

    # Alignment filtering
    cmd_line += '--score_threshold {0:.2f} '.format(args.score_threshold)
    if args.straight_mode:
        cmd_line += '--straight_mode '
    cmd_line += '--coverage_threshold {0} '.format(args.coverage_threshold)

    # Overlap-graph building
    cmd_line += '--min_identity {0:.2f} '.format(args.min_identity)
    cmd_line += '--min_overlap_length {0} '.format(args.min_overlap_length)

    # Graph compaction & Components identification
    cmd_line += '--min_read_node {0} '.format(args.min_read_node)
    cmd_line += '--min_overlap_edge {0} '.format(args.min_overlap_edge)
    if args.optimize_components:
        cmd_line += '--optimize_components '

    # LCA labelling
    cmd_line += '--quorum {0:.2f} '.format(args.quorum)

    # Computing compressed graph stats

    # Contigs assembly
    cmd_line += '--read_correction {0} '.format(args.read_correction)
    if args.read_correction == 'auto':
        cmd_line += '--contig_coverage_threshold {0} '.format(args.contig_coverage_threshold)

    # Scaffolding
    if args.contigs_binning:
        cmd_line += '--contigs_binning '

    cmd_line += '--min_scaffold_length {0} '.format(args.min_scaffold_length)

    # Taxonomic assignment
    if args.perform_taxonomic_assignment:
        cmd_line += '--perform_taxonomic_assignment '
        cmd_line += '--training_model {} '.format(args.training_model)

    # Visualization

    # Main parameters
    cmd_line += '--out_dir {0} '.format(args.out_dir)
    cmd_line += '--ref_db {0} '.format(args.ref_db)
    cmd_line += '--input_fastx {0} '.format(args.input_fastx)

    # Print cmd line
    logger.info('CMD: {0}'.format(cmd_line))

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


def main():
    """
    """
    # Set global t0
    global_t0_wall = time.time()

    # Arguments parsing
    args = parse_arguments()

    # Init error code
    error_code = 0

    # Init list of tmp files to delete at the end
    to_rm_filepath_list = list()

    # Set resume tag
    run_step = True
    if args.resume_from:
        run_step = False

    ##############################################
    # Set all files and directories names + paths

    workdir = os.path.join(args.out_dir, 'workdir')

    try:
        if not os.path.exists(workdir):
            logger.debug('mkdir {0}'.format(workdir))
            os.makedirs(workdir)
    except OSError:
        logger.exception('Could not create output directory {0}'.format(workdir))
        raise

    logger_filepath = os.path.join(args.out_dir, 'matam.log')
    update_logger_settings(logger_filepath, args.verbose, args.debug)

    # Print intro infos
    if not args.filter_only:
        print_intro(args)

    sort_bin = 'sort -T ' + workdir + ' -S ' + str(args.max_memory)
    sort_bin += 'M --parallel ' + str(args.cpu)

    input_fastx_filepath = args.input_fastx
    input_fastx_filename = os.path.basename(input_fastx_filepath)
    input_fastx_basename, input_fastx_extension = os.path.splitext(input_fastx_filename)

    ref_db_basepath = args.ref_db
    ref_db_dir, ref_db_basename = os.path.split(ref_db_basepath)

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
    sortme_output_basepath = os.path.join(workdir, sortme_output_basename)

    sortme_output_fastx_filepath = sortme_output_basepath + input_fastx_extension
    sortme_output_sam_filepath = sortme_output_basepath + '.sam'

    # Alignments filtering
    score_threshold_int = int(args.score_threshold * 100)

    sam_filt_basename = sortme_output_basename + '.scr_filt_'
    if args.straight_mode:
        sam_filt_basename += 'str_'
    else:
        sam_filt_basename += 'geo_'
    sam_filt_basename += str(score_threshold_int) + 'pct'
    sam_filt_filename = sam_filt_basename + '.sam'
    sam_filt_filepath = os.path.join(workdir, sam_filt_filename)

    # Poorly covered references filtering
    sam_cov_filt_basename = sam_filt_basename
    if args.coverage_threshold:
        sam_cov_filt_basename += '.cov_filt_{0}'.format(args.coverage_threshold)
    sam_cov_filt_filename = sam_cov_filt_basename + '.sam'
    sam_cov_filt_filepath = os.path.join(workdir, sam_cov_filt_filename)

    # Overlap-graph building
    min_identity_int = int(args.min_identity * 100)

    ovgraphbuild_basename = sam_cov_filt_basename + '.ovgb_i' + str(min_identity_int)
    ovgraphbuild_basename += '_o' + str(args.min_overlap_length)
    ovgraphbuild_basepath = os.path.join(workdir, ovgraphbuild_basename)

    ovgraphbuild_asqg_filepath = ovgraphbuild_basepath + '.asqg'
    ovgraphbuild_nodes_csv_filepath = ovgraphbuild_basepath + '.nodes.csv'
    ovgraphbuild_edges_csv_filepath = ovgraphbuild_basepath + '.edges.csv'

    # Graph compaction & Components identification
    componentsearch_basename = ovgraphbuild_basename + '.cpts'
    componentsearch_basename += '_N' + str(args.min_read_node)
    componentsearch_basename += '_E' + str(args.min_overlap_edge)
    if args.optimize_components:
        componentsearch_basename += '_oc'
    componentsearch_basepath = os.path.join(workdir, componentsearch_basename)

    componentsearch_metanodes_csv_basename = componentsearch_basename + '.metaNodes'
    componentsearch_metanodes_csv_filename = componentsearch_metanodes_csv_basename + '.csv'
    componentsearch_metanodes_csv_filepath = os.path.join(workdir, componentsearch_metanodes_csv_filename)

    componentsearch_metaedges_csv_basename = componentsearch_basename + '.metaEdges'
    componentsearch_metaedges_csv_filename = componentsearch_metaedges_csv_basename + '.csv'
    componentsearch_metaedges_csv_filepath = os.path.join(workdir, componentsearch_metaedges_csv_filename)

    componentsearch_components_csv_basename = componentsearch_basename + '.components'
    componentsearch_components_csv_filename = componentsearch_components_csv_basename + '.csv'
    componentsearch_components_csv_filepath = os.path.join(workdir, componentsearch_components_csv_filename)

    # LCA labelling
    read_metanode_component_filepath = componentsearch_basepath + '.read_metanode_component.tab'
    complete_taxo_filepath = componentsearch_basepath + '.read_metanode_component_taxo.tab'

    quorum_int = int(args.quorum * 100)

    labelled_nodes_basename = componentsearch_metanodes_csv_basename
    labelled_nodes_basename += '.component_lca' + str(quorum_int) + 'pct'

    components_lca_filename = componentsearch_basename + '.component_lca' + str(quorum_int) + 'pct.tab'
    components_lca_filepath = os.path.join(workdir, components_lca_filename)

    # Computing compressed graph stats
    stats_filename = componentsearch_basename + '.graph.stats'
    stats_filepath = os.path.join(workdir, stats_filename)

    # Contigs assembly
    component_read_filepath = componentsearch_basepath + '.component_read.tab'

    contigs_assembly_wkdir = os.path.join(workdir, "components_assembly")

    contigs_basename = componentsearch_basename + '.'
    contigs_basename += args.assembler + '_by_component'
    contigs_basepath = os.path.join(workdir, contigs_basename)
    contigs_filename = contigs_basename + '.fasta'
    contigs_filepath = os.path.join(workdir, contigs_filename)

    contigs_assembly_log_filename = contigs_basename + '.log'
    contigs_assembly_log_filepath = os.path.join(workdir, contigs_assembly_log_filename)
    # Remove log file from previous assemblies
    # because we will only append to this file
    if os.path.exists(contigs_assembly_log_filepath):
        os.remove(contigs_assembly_log_filepath)

    contigs_symlink_basename = 'contigs'
    contigs_symlink_filename = contigs_symlink_basename + '.fasta'
    contigs_symlink_filepath = os.path.join(workdir, contigs_symlink_filename)

    contigs_NR_basename = contigs_symlink_basename + '.NR'
    contigs_NR_filename = contigs_NR_basename + '.fasta'
    contigs_NR_filepath = os.path.join(workdir, contigs_NR_filename)

    large_NR_contigs_basename = contigs_NR_basename + '.min_' + str(args.min_scaffold_length) + 'bp'
    large_NR_contigs_filename = large_NR_contigs_basename + '.fasta'
    large_NR_contigs_filepath = os.path.join(workdir, large_NR_contigs_filename)

    # Scaffolding
    scaff_sortme_output_basename = contigs_symlink_basename
    scaff_sortme_output_basename += '.sortmerna_vs_complete_' + ref_db_basename
    scaff_sortme_output_basename += '_num_align_0'
    scaff_sortme_output_basepath = os.path.join(workdir, scaff_sortme_output_basename)

    scaff_sortme_output_blast_filename = scaff_sortme_output_basename + '.blast'
    scaff_sortme_output_blast_filepath = os.path.join(workdir, scaff_sortme_output_blast_filename)

    scaff_sortme_output_sam_filename = scaff_sortme_output_basename + '.sam'
    scaff_sortme_output_sam_filepath = os.path.join(workdir, scaff_sortme_output_sam_filename)

    best_only_blast_basename = scaff_sortme_output_blast_filename + '.best_only'
    best_only_blast_filename = best_only_blast_basename + '.tab'
    best_only_blast_filepath = os.path.join(workdir, best_only_blast_filename)

    selected_best_only_blast_basename = best_only_blast_basename + '.selected'
    selected_best_only_blast_filename = selected_best_only_blast_basename + '.tab'
    selected_best_only_blast_filepath = os.path.join(workdir, selected_best_only_blast_filename)

    selected_sam_basename = selected_best_only_blast_basename
    selected_sam_filename = selected_sam_basename + '.sam'
    selected_sam_filepath = os.path.join(workdir, selected_sam_filename)

    binned_sam_basename = selected_best_only_blast_basename + '.binned'
    binned_sam_filename = binned_sam_basename + '.sam'
    binned_sam_filepath = os.path.join(workdir, binned_sam_filename)

    processed_sam_basename = selected_sam_basename
    processed_sam_filepath = selected_sam_filepath

    if args.contigs_binning:
        processed_sam_basename = binned_sam_basename
        processed_sam_filepath = binned_sam_filepath

    bam_basename = processed_sam_basename
    bam_filename = bam_basename + '.bam'
    bam_filepath = os.path.join(workdir, bam_filename)

    sorted_bam_basename = processed_sam_basename + '.sorted'
    sorted_bam_basepath = os.path.join(workdir, sorted_bam_basename)
    sorted_bam_filename = sorted_bam_basename + '.bam'
    sorted_bam_filepath = os.path.join(workdir, sorted_bam_filename)

    mpileup_filename = sorted_bam_basename + '.mpileup'
    mpileup_filepath = os.path.join(workdir, mpileup_filename)

    scaffolds_basename = sorted_bam_basename + '.scaffolds'
    scaffolds_filename = scaffolds_basename + '.fa'
    scaffolds_filepath = os.path.join(workdir, scaffolds_filename)

    scaffolds_symlink_basename = 'scaffolds'
    if args.contigs_binning:
        scaffolds_symlink_basename += '.contigs_binning'
    scaffolds_symlink_filename = scaffolds_symlink_basename + '.fa'
    scaffolds_symlink_filepath = os.path.join(workdir, scaffolds_symlink_filename)

    scaffolds_NR_basename = scaffolds_symlink_basename + '.NR'
    scaffolds_NR_filename = scaffolds_NR_basename + '.fa'
    scaffolds_NR_filepath = os.path.join(workdir, scaffolds_NR_filename)

    large_NR_scaffolds_basename = scaffolds_NR_basename + '.min_' + str(args.min_scaffold_length) + 'bp'
    large_NR_scaffolds_filename = large_NR_scaffolds_basename + '.fa'
    large_NR_scaffolds_filepath = os.path.join(workdir, large_NR_scaffolds_filename)

    final_assembly_symlink_basename = 'final_assembly'
    final_assembly_symlink_filename = final_assembly_symlink_basename + '.fa'
    final_assembly_symlink_filepath = os.path.join(args.out_dir, final_assembly_symlink_filename)

    final_krona_tab_symlink_filepath = os.path.join(args.out_dir, 'krona.tab')
    final_krona_html_symlink_filepath = os.path.join(args.out_dir, 'krona.html')

    #################################
    # Compute input reads statistics

    if args.resume_from == 'reads_mapping':
        run_step = True

    # Get input reads number
    input_reads_nb = -1

    if run_step:
        if input_fastx_extension in ('.fq', '.fastq'):
            input_fastx_line_nb = int(subprocess.check_output('wc -l {0}'.format(input_fastx_filepath), shell=True, bufsize=0).split()[0])
            if input_fastx_line_nb % 4 != 0:
                logger.warning('FastQ input file does not have a number of lines multiple of 4')
            input_reads_nb = input_fastx_line_nb // 4
        elif input_fastx_extension in ('.fa', '.fasta'):
            input_reads_nb = int(subprocess.check_output('grep -c "^>" {0}'.format(input_fastx_filepath), shell=True, bufsize=0))
        else:
            logger.warning('Input fastx file extension was not recognised ({0})'.format(input_fastx_extension))

        logger.info('=== Input ===')
        logger.info('Input file: {}'.format(input_fastx_filepath))
        logger.info('Input file reads nb: {} reads'.format(input_reads_nb))

    ###############################
    # Reads mapping against ref db

    if args.resume_from == 'reads_mapping':
        logger.info('Resuming from reads mapping')
        run_step = True

    if run_step:
        logger.info('=== Reads mapping against ref db ===')

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

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Reads mapping completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Get selected reads number
    selected_reads_nb = -1

    if input_fastx_extension in ('.fq', '.fastq'):
        selected_fastx_line_nb = int(subprocess.check_output('wc -l {0}'.format(sortme_output_fastx_filepath), shell=True, bufsize=0).split()[0])
        selected_reads_nb = selected_fastx_line_nb // 4
    elif input_fastx_extension in ('.fa', '.fasta'):
        selected_reads_nb = int(subprocess.check_output('grep -c "^>" {0}'.format(sortme_output_fastx_filepath), shell=True, bufsize=0))

    logger.info('Identified as marker: {} / {} reads ({:.2f}%)'.format(selected_reads_nb, input_reads_nb, selected_reads_nb*100.0/input_reads_nb))

    # Tag tmp files for removal
    to_rm_filepath_list.append(sortme_output_basepath + '.log')
    to_rm_filepath_list.append(sortme_output_basepath + '.blast')

    if args.filter_only:
        return error_code
    #############################
    # Alignment filtering

    if args.resume_from == 'alignments_filtering':
        logger.info('Resuming from alignments filtering')
        run_step = True

    if run_step:
        logger.info('=== Alignment filtering ===')

        cmd_line = 'cat ' + sortme_output_sam_filepath
        cmd_line += ' | grep -v "^@" | '
        cmd_line += sort_bin + ' -k 1,1V -k 12,12nr'
        cmd_line += ' | ' + filter_score_bin + ' -t ' + str(args.score_threshold)
        if not args.straight_mode:
            cmd_line += ' --geometric'
        cmd_line += ' > ' + sam_filt_filepath

        # Set t0
        t0_wall = time.time()

        logger.debug('CMD: {0}'.format(cmd_line))
        runner.logged_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Good alignments filtering completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        if args.coverage_threshold:
            cmd_line = 'cat ' + sam_filt_filepath
            cmd_line += ' | ' + sort_bin + ' -k 3,3 -k 4,4n'
            cmd_line += ' | ' + sample_sam_cov_bin + ' -c ' + str(args.coverage_threshold)
            cmd_line += ' -r ' + clustered_ref_db_filepath
            cmd_line += ' | ' + sort_bin + ' -k 1,1V -k 12,12nr'
            cmd_line += ' > ' + sam_cov_filt_filepath

            # Set t0
            t0_wall = time.time()

            runner.logged_check_call(cmd_line, verbose=args.verbose)

            # Output running time
            logger.info('Ref coverage filtering completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    # Tag tmp files for removal
    if args.coverage_threshold:
        to_rm_filepath_list.append(sam_filt_filepath)

    #########################
    # Overlap-graph building

    if args.resume_from == 'overlap_graph_building':
        logger.info('Resuming from overlap-graph building')
        run_step = True

    ovgraph_nodes_nb = -1
    ovgraph_edges_nb = -1

    if run_step:
        logger.info('=== Overlap-graph building ===')

        cmd_line = ovgraphbuild_bin
        cmd_line += ' -i ' + str(args.min_identity)
        cmd_line += ' -m ' + str(args.min_overlap_length)
        if args.debug:
            cmd_line += ' --asqg'
        cmd_line += ' --csv --output_basename '
        cmd_line += ovgraphbuild_basepath
        cmd_line += ' -r ' + clustered_ref_db_filepath
        cmd_line += ' -s ' + sam_cov_filt_filepath
        if args.verbose:
            cmd_line += ' -v'
        if args.debug:
            cmd_line += ' --debug'

        # Set t0
        t0_wall = time.time()

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Overlap-graph building completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # Tag tmp files for removal
        to_rm_filepath_list.append(ovgraphbuild_asqg_filepath)

    # Get overlap graph stats
    ovgraph_nodes_nb = int(subprocess.check_output('wc -l {0}'.format(ovgraphbuild_nodes_csv_filepath), shell=True, bufsize=0).split()[0]) - 1
    ovgraph_edges_nb = int(subprocess.check_output('wc -l {0}'.format(ovgraphbuild_edges_csv_filepath), shell=True, bufsize=0).split()[0]) - 1

    logger.info('Overlap graph stats: {} nodes, {} edges'.format(ovgraph_nodes_nb, ovgraph_edges_nb))

    ###############################################
    # Graph compaction & Components identification

    if args.resume_from == 'graph_compaction':
        logger.info('Resuming from graph compaction & components identification')
        run_step = True

    if run_step:
        logger.info('=== Graph compaction & Components identification ===')

        cmd_line = componentsearch_bin
        if args.optimize_components:
            cmd_line += ' -o'
        if args.seed:
            cmd_line += ' -r ' + str(args.seed)
        cmd_line += ' -N ' + str(args.min_read_node)
        cmd_line += ' -E ' + str(args.min_overlap_edge)
        cmd_line += ' -b ' + componentsearch_basepath
        cmd_line += ' -n ' + ovgraphbuild_nodes_csv_filepath
        cmd_line += ' -e ' + ovgraphbuild_edges_csv_filepath

        # Set t0
        t0_wall = time.time()

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Graph compaction & Components identification completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # Tag tmp files for removal
        #to_rm_filepath_list.append(ovgraphbuild_nodes_csv_filepath)
        #to_rm_filepath_list.append(ovgraphbuild_edges_csv_filepath)

    # Get compressed graph stats
    compressed_graph_nodes_nb = int(subprocess.check_output('wc -l {0}'.format(componentsearch_metanodes_csv_filepath), shell=True, bufsize=0).split()[0]) - 1
    compressed_graph_edges_nb = int(subprocess.check_output('wc -l {0}'.format(componentsearch_metaedges_csv_filepath), shell=True, bufsize=0).split()[0]) - 1
    compressed_graph_reads_nb = int(subprocess.check_output('cut -d ";" -f5 {0} | awk "\$1!=-1" | wc -l'.format(componentsearch_components_csv_filepath), shell=True, bufsize=0).split()[0]) - 1
    compressed_graph_excluded_reads_nb = ovgraph_nodes_nb - compressed_graph_reads_nb
    excluded_reads_percent = compressed_graph_excluded_reads_nb * 100.0 / ovgraph_nodes_nb
    components_nb = int(subprocess.check_output('cut -d ";" -f5 {0} | {1} | uniq | wc -l'.format(componentsearch_components_csv_filepath, sort_bin), shell=True).split()[0]) - 1

    logger.info('Compressed graph: {} components'.format(components_nb))

    ################
    # LCA labelling

    if run_step:
        logger.info('=== LCA labelling ===')

        # Set t0
        t0_wall = time.time()

        # Get read, metanode and component from the components file
        cmd_line = 'tail -n +2 ' + componentsearch_components_csv_filepath
        cmd_line += ' | sed "s/;/\\t/g" | cut -f2,4,5 | awk "\$3>=0" '
        cmd_line += ' | ' + sort_bin + ' -k1,1 > '
        cmd_line += read_metanode_component_filepath


        error_code += runner.logged_call(cmd_line, verbose=args.verbose)

        # Join sam file with ref taxo on ref name, and join it to the
        # component-node file on read name.
        # Output is (read, metanode, component, taxo) file
        cmd_line = 'cat ' + sam_filt_filepath + ' | cut -f1,3 | ' + sort_bin + ' -k2,2'
        cmd_line += ' | join -12 -21 - ' + complete_ref_db_taxo_filepath
        cmd_line += ' | ' + sort_bin + ' -k2,2 | awk "{print \$2,\$3}" | sed "s/ /\\t/g" '
        cmd_line += ' | join -11 -21 ' + read_metanode_component_filepath
        cmd_line += ' - | sed "s/ /\\t/g" > '
        cmd_line += complete_taxo_filepath


        error_code += runner.logged_call(cmd_line, verbose=args.verbose)

        # Compute LCA at component level using quorum threshold
        cmd_line = 'cat ' + complete_taxo_filepath + ' | ' + sort_bin + ' -k3,3 -k1,1 | '
        cmd_line += compute_lca_bin + ' -t 4 -f 3 -g 1 -m ' + str(args.quorum)
        cmd_line += ' -o ' + components_lca_filepath

        error_code += runner.logged_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('LCA labelling completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # Tag tmp files for removal
        to_rm_filepath_list.append(componentsearch_basepath + '.components.tab')
        to_rm_filepath_list.append(complete_taxo_filepath)

    # ###################################
    # # Computing compressed graph stats

    # if run_step:
    #     logger.info('=== Computing compressed graph stats ===')

    #     cmd_line = compute_compressed_graph_stats_bin + ' --nodes_contracted '
    #     cmd_line += contracted_nodes_filepath + ' --edges_contracted '
    #     cmd_line += contracted_edges_filepath + ' --components_lca '
    #     cmd_line += components_lca_filepath
    #     if args.true_ref_taxo:
    #         cmd_line += ' --species_taxo ' + args.true_ref_taxo
    #         cmd_line += ' --read_node_component ' + read_id_metanode_component_filepath
    #     cmd_line += ' -o ' + stats_filepath

    #     # Set t0
    #     t0_wall = time.time()

    #     logger.debug('CMD: {0}'.format(cmd_line))
    #     error_code += subprocess.call(cmd_line, shell=True, bufsize=0)

    #     # Output running time
    #     logger.info('Computing compressed graph stats completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))
    #     if args.verbose:
    #         sys.stderr.write('\n')

    ###################
    # Contigs assembly

    if args.resume_from == 'contigs_assembly':
        logger.info('Resuming from contigs assembly')
        run_step = True

    if run_step:
        logger.info('=== Contigs assembly ===')

        # Set t0
        t0_wall = time.time()

        components_assembly.assemble_all_components(args.assembler,
                                                    sortme_output_fastx_filepath, read_metanode_component_filepath, components_lca_filepath,
                                                    contigs_filepath, contigs_assembly_wkdir,
                                                    args.cpu, args.read_correction, args.contig_coverage_threshold)

        if not args.keep_tmp:
            shutil.rmtree(contigs_assembly_wkdir)

        # Create symbolic link
        if os.path.exists(contigs_symlink_filepath):
            os.remove(contigs_symlink_filepath)
        os.symlink(os.path.basename(contigs_filepath), contigs_symlink_filepath)

        # TO DO, if it gets better results:
        # remove redundant sequences in contigs

        # Remove redundant contigs
        cmd_line = remove_redundant_bin + ' -i ' + contigs_symlink_filepath
        cmd_line += ' -o ' + contigs_NR_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Filter out small contigs
        cmd_line = fasta_length_filter_bin + ' -m ' + str(args.min_scaffold_length)
        cmd_line += ' -i ' + contigs_NR_filepath
        cmd_line += ' -o ' + large_NR_contigs_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Contigs assembly completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # Evaluate assembly if true ref are provided
        if args.true_references:

            true_ref_filename = os.path.basename(args.true_references)
            true_ref_basename, true_ref_extension = os.path.splitext(true_ref_filename)

            # Evaluate all contigs
            cmd_line = evaluate_assembly_bin + ' -r ' + args.true_references
            cmd_line += ' -i ' + contigs_symlink_filepath

            runner.logged_check_call(cmd_line, verbose=args.verbose)

            contigs_assembly_stats_filename = contigs_symlink_basename + '.exonerate_vs_'
            contigs_assembly_stats_filename += true_ref_basename + '.assembly.stats'
            contigs_assembly_stats_filepath = os.path.join(workdir, contigs_assembly_stats_filename)

            contigs_error_rate = float(subprocess.check_output('grep "error rate   =" {0}'.format(contigs_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
            contigs_error_rate_2 = float(subprocess.check_output('grep "error rate 2 =" {0}'.format(contigs_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
            contigs_ref_coverage = float(subprocess.check_output('grep "ref coverage" {0}'.format(contigs_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])

            # Evaluate large NR contigs
            cmd_line = evaluate_assembly_bin + ' -r ' + args.true_references
            cmd_line += ' -i ' + large_NR_contigs_filepath

            runner.logged_check_call(cmd_line, verbose=args.verbose)

            large_NR_contigs_assembly_stats_filename = large_NR_contigs_basename + '.exonerate_vs_'
            large_NR_contigs_assembly_stats_filename += true_ref_basename + '.assembly.stats'
            large_NR_contigs_assembly_stats_filepath = os.path.join(workdir, large_NR_contigs_assembly_stats_filename)

            try:
                large_NR_contigs_error_rate = float(subprocess.check_output('grep "error rate   =" {0}'.format(large_NR_contigs_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
                large_NR_contigs_error_rate_2 = float(subprocess.check_output('grep "error rate 2 =" {0}'.format(large_NR_contigs_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
                large_NR_contigs_ref_coverage = float(subprocess.check_output('grep "ref coverage" {0}'.format(large_NR_contigs_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
            except subprocess.CalledProcessError:
                large_NR_contigs_error_rate = float()
                large_NR_contigs_error_rate_2 = float()
                large_NR_contigs_ref_coverage = float()

        # Tag tmp files for removal
        to_rm_filepath_list.append(read_metanode_component_filepath)
        to_rm_filepath_list.append(contigs_NR_filepath)

    # Compute contigs assembly stats
    contigs_stats = compute_fasta_stats(contigs_filepath)
    large_NR_contigs_stats = compute_fasta_stats(large_NR_contigs_filepath)

    ##############
    # Scaffolding

    if args.resume_from == 'scaffolding':
        logger.info('Resuming from scaffolding')
        run_step = True

    if run_step:
        logger.info('=== Scaffolding ===')

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
        #cmd_line += ' --num_seeds 3 '
        #cmd_line += ' --best 0 --min_lis 10 '
        cmd_line += ' -e {0:.2e}'.format(scaff_evalue)
        cmd_line += ' -a ' + str(args.cpu)
        if args.verbose:
            cmd_line += ' -v '

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Contig mapping completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # Set t0
        t0_wall = time.time()

        # Keep only quasi-equivalent best matches for each contig
        # -p 0.99 is used because of the tendency of SortMeRNA to soft-clip
        # a few nucleotides at 5’ and 3’ ends
        cmd_line = sort_bin + ' -k1,1V -k12,12nr ' + scaff_sortme_output_blast_filepath
        cmd_line += ' | ' + get_best_matches_bin + ' -p 0.99 -o ' + best_only_blast_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Select blast matches for scaffolding using a specific-first conserved-later approach
        cmd_line = gener_scaff_blast_bin + ' -i ' + best_only_blast_filepath
        cmd_line += ' -o ' + selected_best_only_blast_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Filter sam file based on blast scaffolding file
        cmd_line = filter_sam_blast_bin + ' -i ' + scaff_sortme_output_sam_filepath
        cmd_line += ' -b ' +  selected_best_only_blast_filepath
        cmd_line += ' | ' + sort_bin + ' -k3,3 -k4,4n > ' + selected_sam_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Bin compatible contigs matching on the same reference
        if args.contigs_binning:
            cmd_line = compute_contigs_compatibility_bin
            cmd_line += ' -i ' + selected_sam_filepath
            cmd_line += ' | ' + sort_bin + ' -k1,1n > ' + binned_sam_filepath

            runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Convert sam to bam
        cmd_line = 'samtools view -b -S ' + processed_sam_filepath
        if not args.contigs_binning:
            cmd_line += ' -T ' + complete_ref_db_filepath
        cmd_line += ' -o ' + bam_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Sort bam
        cmd_line = 'samtools sort -o ' + sorted_bam_filepath + ' ' + bam_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Generate mpileup
        cmd_line = 'samtools mpileup -d 10000 -o ' + mpileup_filepath
        cmd_line += ' ' + sorted_bam_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Scaffold contigs based on mpileup
        cmd_line = scaffold_contigs_bin + ' -i ' + mpileup_filepath
        cmd_line += ' -o ' + scaffolds_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Create symbolic link
        if os.path.exists(scaffolds_symlink_filepath):
            os.remove(scaffolds_symlink_filepath)
        os.symlink(os.path.basename(scaffolds_filepath), scaffolds_symlink_filepath)

        # Remove redundant scaffolds
        cmd_line = remove_redundant_bin + ' -i ' + scaffolds_symlink_filepath
        cmd_line += ' -o ' + scaffolds_NR_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Filter out small scaffolds
        cmd_line = fasta_length_filter_bin + ' -m ' + str(args.min_scaffold_length)
        cmd_line += ' -i ' + scaffolds_NR_filepath
        cmd_line += ' -o ' + large_NR_scaffolds_filepath

        runner.logged_check_call(cmd_line, verbose=args.verbose)

        # Output running time
        logger.info('Scaffolding completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # check that the the resulting file is not empty
        is_empty = True
        try:
            is_empty = os.path.getsize(large_NR_scaffolds_filepath) == 0
        except OSError:
            is_empty = True

        if is_empty:
            logger.fatal("Can't proceed all steps (abundance calculation ...) when no scaffolds are available: %s" % large_NR_scaffolds_filepath)
            sys.exit("Can't proceed further")

        # Evaluate assembly if true ref are provided
        if args.true_references:

            true_ref_filename = os.path.basename(args.true_references)
            true_ref_basename, true_ref_extension = os.path.splitext(true_ref_filename)

            # Evaluate all scaffolds
            cmd_line = evaluate_assembly_bin + ' -r ' + args.true_references
            cmd_line += ' -i ' + scaffolds_symlink_filepath

            runner.logged_check_call(cmd_line, verbose=args.verbose)

            scaffolds_assembly_stats_filename = scaffolds_symlink_basename + '.exonerate_vs_'
            scaffolds_assembly_stats_filename += true_ref_basename + '.assembly.stats'
            scaffolds_assembly_stats_filepath = os.path.join(workdir, scaffolds_assembly_stats_filename)

            scaffolds_error_rate = float(subprocess.check_output('grep "error rate   =" {0}'.format(scaffolds_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
            scaffolds_error_rate_2 = float(subprocess.check_output('grep "error rate 2 =" {0}'.format(scaffolds_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
            scaffolds_ref_coverage = float(subprocess.check_output('grep "ref coverage" {0}'.format(scaffolds_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])

            # Evaluate large NR scaffolds
            cmd_line = evaluate_assembly_bin + ' -r ' + args.true_references
            cmd_line += ' -i ' + large_NR_scaffolds_filepath

            runner.logged_check_call(cmd_line, verbose=args.verbose)

            large_NR_scaffolds_assembly_stats_filename = large_NR_scaffolds_basename + '.exonerate_vs_'
            large_NR_scaffolds_assembly_stats_filename += true_ref_basename + '.assembly.stats'
            large_NR_scaffolds_assembly_stats_filepath = os.path.join(workdir, large_NR_scaffolds_assembly_stats_filename)

            try:
                large_NR_scaffolds_error_rate = float(subprocess.check_output('grep "error rate   =" {0}'.format(large_NR_scaffolds_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
                large_NR_scaffolds_error_rate_2 = float(subprocess.check_output('grep "error rate 2 =" {0}'.format(large_NR_scaffolds_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
                large_NR_scaffolds_ref_coverage = float(subprocess.check_output('grep "ref coverage" {0}'.format(large_NR_scaffolds_assembly_stats_filepath), shell=True, bufsize=0).decode("utf-8").split('=')[1].strip()[:-1])
            except subprocess.CalledProcessError:
                large_NR_scaffolds_error_rate = float()
                large_NR_scaffolds_error_rate_2 = float()
                large_NR_scaffolds_ref_coverage = float()

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
        to_rm_filepath_list.append(scaffolds_NR_filepath)

    # Compute scaffolds assemblies stats
    scaffolds_stats = compute_fasta_stats(scaffolds_filepath)
    large_NR_scaffolds_stats = compute_fasta_stats(large_NR_scaffolds_filepath)

    ########################
    # Abundance calculation

    scaffolds_fasta = large_NR_scaffolds_filepath
    reads = sortme_output_fastx_filepath
    fasta_with_abundance_filepath =  '%s.abd%s' % os.path.splitext(scaffolds_fasta)
    abundance = None

    if args.resume_from == 'abundance_calculation':
        logger.info('Resuming from abundance calculation')
        run_step = True

    if run_step:

        logger.info('=== Abundance calculation ===')
        # Set t0
        t0_wall = time.time()

        abundance = get_abundance_by_scaffold(indexdb_bin, sortmerna_bin, get_best_matches_bin,
                                              scaffolds_fasta, reads,
                                              args.best, args.min_lis, args.evalue,
                                              args.max_memory, args.cpu,
                                              output_dir_basepath=workdir,
                                              verbose=args.verbose,
                                              keep_tmp=args.keep_tmp
        )


        complete_fasta_with_abundance(scaffolds_fasta, fasta_with_abundance_filepath, abundance)

        logger.debug('Write abundance informations to: %s' % fasta_with_abundance_filepath)
        # Output running time
        logger.info('Abundance calculation completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

    #################################
    # taxonomic assignment with rdp
    if args.resume_from == 'taxonomic_assignment':
        logger.info('Resuming from Taxonomic assignment ')
        run_step = True
        #If the user want to resume from this step, switch the args.perform_taxonomic_assignment
        #to true in case the user forget it.
        args.perform_taxonomic_assignment = True


    if run_step and args.perform_taxonomic_assignment:
        logger.info('=== Taxonomic assignment & Krona visualization ===')
        # Set t0
        t0_wall = time.time()
        rdp_cutoff = 0.8
        rdp_classification_filepath =  '%s.rdp.tab' % os.path.splitext(fasta_with_abundance_filepath)[0]
        run_rdp_classifier(rdp_exe, fasta_with_abundance_filepath,
                           rdp_classification_filepath, gene=args.training_model, cutoff=rdp_cutoff)
        logger.debug('Write taxonomic assignment to: %s' % rdp_classification_filepath)

        # tag results below the confidence cutoff as unclassified
        fltr_rdp_classification_filepath = '%s.fltr.tab' % os.path.splitext(rdp_classification_filepath)[0]
        filter_rdp_file(rdp_classification_filepath, fltr_rdp_classification_filepath, cutoff=rdp_cutoff)

        final_rdp_tab_symlink_filepath = os.path.join(args.out_dir, 'rdp.tab')
        force_symlink(
            os.path.relpath(fltr_rdp_classification_filepath, start=args.out_dir),
            final_rdp_tab_symlink_filepath,
        )
        #############################
        # build krona representation

        if not abundance:
            abundance = get_abundance_from_fasta(fasta_with_abundance_filepath)
        krona_text_filepath = '%s.krona.tab' % os.path.splitext(fltr_rdp_classification_filepath)[0]
        rdp_file_to_krona_text_file(fltr_rdp_classification_filepath, krona_text_filepath, abundance=abundance)

        krona_html_filepath =  '%s.html' % os.path.splitext(krona_text_filepath)[0]
        make_krona_plot(krona_bin, krona_text_filepath, krona_html_filepath)

        logger.debug('Write krona to: %s' % krona_html_filepath)
        logger.info('Taxonomic assignment & Krona visualization completed in {0:.4f} seconds wall time'.format(time.time() - t0_wall))

        # Expose final files
        force_symlink(
            os.path.relpath(krona_text_filepath, start=args.out_dir),
            final_krona_tab_symlink_filepath,
        )
        force_symlink(
            os.path.relpath(krona_html_filepath, start=args.out_dir),
            final_krona_html_symlink_filepath,
        )
    force_symlink(
        os.path.relpath(fasta_with_abundance_filepath, start=args.out_dir),
        final_assembly_symlink_filepath,
    )
    ###########################
    # Print Assembly Statistics

    if args.verbose:
        b = '=== MATAM Statistics ===\n\n'
        b += 'Input reads nb: {0}\n'.format(input_reads_nb)
        b += 'Selected reads nb: {0}\n\n'.format(selected_reads_nb)

        b += 'Ov. Graph nodes nb: {0}\n'.format(ovgraph_nodes_nb)
        b += 'Ov. Graph edges nb: {0}\n\n'.format(ovgraph_edges_nb)

        b += 'Compressed Graph nodes nb: {0}\n'.format(compressed_graph_nodes_nb)
        b += 'Compressed Graph edges nb: {0}\n'.format(compressed_graph_edges_nb)
        b += 'Compressed Graph reads nb: {0}\n'.format(compressed_graph_reads_nb)
        b += 'Compr. Graph excluded reads nb: {0} '.format(compressed_graph_excluded_reads_nb)
        b += '({0:.2f}% of total reads)\n'.format(excluded_reads_percent)
        b += 'Compr. Graph components nb: {0}\n\n'.format(components_nb)

        b += 'Contigs assembly:\n'
        b += '\tSeq nb: {0}\n'.format(contigs_stats.seq_num)
        b += '\tSeq min size: {0} nt\n'.format(contigs_stats.get_min_length())
        b += '\tSeq max size: {0} nt\n'.format(contigs_stats.get_max_length())
        b += '\tSeq avg. size: {0:.2f} nt\n'.format(contigs_stats.get_avg_length())
        b += '\tSeq total size: {0} nt\n'.format(contigs_stats.total_nt)
        if args.true_references:
            b += '\tError rate: {0:.2f}%\n'.format(contigs_error_rate)
            b += '\tError rate 2: {0:.2f}%\n'.format(contigs_error_rate_2)
            b += '\tRef coverage: {0:.2f}%\n'.format(contigs_ref_coverage)
        b += '\n'

        b += 'Large NR Contigs assembly:\n'
        b += '\tSeq nb: {0}\n'.format(large_NR_contigs_stats.seq_num)
        b += '\tSeq min size: {0} nt\n'.format(large_NR_contigs_stats.get_min_length())
        b += '\tSeq max size: {0} nt\n'.format(large_NR_contigs_stats.get_max_length())
        b += '\tSeq avg. size: {0:.2f} nt\n'.format(large_NR_contigs_stats.get_avg_length())
        b += '\tSeq total size: {0} nt\n'.format(large_NR_contigs_stats.total_nt)
        if args.true_references:
            b += '\tError rate: {0:.2f}%\n'.format(large_NR_contigs_error_rate)
            b += '\tError rate 2: {0:.2f}%\n'.format(large_NR_contigs_error_rate_2)
            b += '\tRef coverage: {0:.2f}%\n'.format(large_NR_contigs_ref_coverage)
        b += '\n'

        b += 'Scaffolds assembly:\n'
        b += '\tSeq nb: {0}\n'.format(scaffolds_stats.seq_num)
        b += '\tSeq min size: {0} nt\n'.format(scaffolds_stats.get_min_length())
        b += '\tSeq max size: {0} nt\n'.format(scaffolds_stats.get_max_length())
        b += '\tSeq avg. size: {0:.2f} nt\n'.format(scaffolds_stats.get_avg_length())
        b += '\tSeq total size: {0} nt\n'.format(scaffolds_stats.total_nt)
        if args.true_references:
            b += '\tError rate: {0:.2f}%\n'.format(scaffolds_error_rate)
            b += '\tError rate 2: {0:.2f}%\n'.format(scaffolds_error_rate_2)
            b += '\tRef coverage: {0:.2f}%\n'.format(scaffolds_ref_coverage)
        b += '\n'

        b += 'Large NR Scaffolds assembly:\n'
        b += '\tSeq nb: {0}\n'.format(large_NR_scaffolds_stats.seq_num)
        b += '\tSeq min size: {0} nt\n'.format(large_NR_scaffolds_stats.get_min_length())
        b += '\tSeq max size: {0} nt\n'.format(large_NR_scaffolds_stats.get_max_length())
        b += '\tSeq avg. size: {0:.2f} nt\n'.format(large_NR_scaffolds_stats.get_avg_length())
        b += '\tSeq total size: {0} nt\n'.format(large_NR_scaffolds_stats.total_nt)
        if args.true_references:
            b += '\tError rate: {0:.2f}%\n'.format(large_NR_scaffolds_error_rate)
            b += '\tError rate 2: {0:.2f}%\n'.format(large_NR_scaffolds_error_rate_2)
            b += '\tRef coverage: {0:.2f}%\n'.format(large_NR_scaffolds_ref_coverage)
        b += '\n'

        if args.debug:
            b += 'One-line stats\n'
            b += '{}\t{}\t{}\t{}\t'.format(input_reads_nb, selected_reads_nb,
                                           ovgraph_nodes_nb, ovgraph_edges_nb)
            b += '{}\t{}\t{}\t{}\t{:.2f}%\t{}\t'.format(compressed_graph_nodes_nb, compressed_graph_edges_nb,
                                                       compressed_graph_reads_nb, compressed_graph_excluded_reads_nb,
                                                       excluded_reads_percent, components_nb)
            b += '{}\t{}\t{}\t{:.2f}\t{}\t'.format(contigs_stats.seq_num, contigs_stats.get_min_length(), contigs_stats.get_max_length(),
                                                  contigs_stats.get_avg_length(), contigs_stats.total_nt)
            if args.true_references:
                b += '{:.2f}%\t{:.2f}%\t{:.2f}%\t'.format(contigs_error_rate, contigs_error_rate_2, contigs_ref_coverage)
            if large_NR_contigs_stats.seq_num:
                b += '{}\t{}\t{}\t{:.2f}\t{}\t'.format(large_NR_contigs_stats.seq_num, large_NR_contigs_stats.get_min_length(), large_NR_contigs_stats.get_max_length(),
                                                      large_NR_contigs_stats.get_avg_length(), large_NR_contigs_stats.total_nt)
                if args.true_references:
                    b += '{:.2f}%\t{:.2f}%\t{:.2f}%\t'.format(large_NR_contigs_error_rate, large_NR_contigs_error_rate_2, large_NR_contigs_ref_coverage)
            else:
                b += '0\tNA\tNA\tNA\t0\t'
                if args.true_references:
                    b += 'NA\tNA\t0%\t'
            b += '{}\t{}\t{}\t{:.2f}\t{}\t'.format(scaffolds_stats.seq_num, scaffolds_stats.get_min_length(), scaffolds_stats.get_max_length(),
                                                  scaffolds_stats.get_avg_length(), scaffolds_stats.total_nt)
            if args.true_references:
                b += '{:.2f}%\t{:.2f}%\t{:.2f}%\t'.format(scaffolds_error_rate, scaffolds_error_rate_2, scaffolds_ref_coverage)
            if large_NR_scaffolds_stats.seq_num:
                b += '{}\t{}\t{}\t{:.2f}\t{}\t'.format(large_NR_scaffolds_stats.seq_num, large_NR_scaffolds_stats.get_min_length(), large_NR_scaffolds_stats.get_max_length(),
                                                      large_NR_scaffolds_stats.get_avg_length(), large_NR_scaffolds_stats.total_nt)
                if args.true_references:
                    b += '{:.2f}%\t{:.2f}%\t{:.2f}%\t'.format(large_NR_scaffolds_error_rate, large_NR_scaffolds_error_rate_2, large_NR_scaffolds_ref_coverage)
            else:
                b += '0\tNA\tNA\tNA\t0'
                if args.true_references:
                    b += '\tNA\tNA\t0%'
            b += '\n'

        logger.info(b)

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
        logger.info('{0} terminated with no error'.format(program_filename))
    # Deal with errors
    else:
        logger.error('{0} terminated with some errors. '.format(program_filename))
        if args.verbose:
            logger.info('Check the log for additional infos')
        else:
            logger.info('Rerun the program using --verbose or --debug option')
        exit_code = 1

    logger.info('Run completed in {0:.4f} seconds wall time'.format(time.time() - global_t0_wall))

    return exit_code


if __name__ == '__main__':

    # Set logging

    import logging.config
    logging.config.dictConfig(logging_config)

    exit_code = 0

    # cProfile.run('main()')
    exit_code = main()

    exit(exit_code)
