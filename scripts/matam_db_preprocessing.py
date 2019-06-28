#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import time
import logging
from binary_utils import Binary

# Create logger
logger = logging.getLogger(__name__)

# Get program filename
program_filename = os.path.basename(sys.argv[0])

# Get MATAM root dir absolute path
matam_db_prepro_bin = os.path.realpath(sys.argv[0])
matam_bin_dir = os.path.dirname(matam_db_prepro_bin)
matam_root_dir = os.path.dirname(matam_bin_dir)

# Get all dependencies bin
matam_script_dir = os.path.join(matam_root_dir, 'scripts')
extract_taxo_bin = os.path.join(matam_script_dir, 'extract_taxo_from_fasta.py')
replace_Ns_bin = os.path.join(matam_script_dir, 'replace_Ns_by_As.py')
sort_fasta_bin = os.path.join(matam_script_dir, 'sort_fasta_by_length.py')
fasta_length_filter_bin = os.path.join(matam_script_dir, 'fasta_length_filter.py')
fasta_name_filter_bin = os.path.join(matam_script_dir, 'fasta_name_filter.py')
clean_name_bin = os.path.join(matam_script_dir, 'fasta_clean_name.py')
indexdb_bin = Binary.assert_which('indexdb_rna')
vsearch_bin = Binary.assert_which('vsearch')

# Define a null file handle
FNULL = open(os.devnull, 'w')


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
    parser = DefaultHelpParser(description='MATAM db preprocessing',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))

    # Main parameters
    group_main = parser.add_argument_group('Main parameters')
    # -i / --input_ref
    group_main.add_argument('-i', '--input_ref',
                            action = 'store',
                            metavar = 'INREF',
                            type = str,
                            required = True,
                            help = 'Input reference file (fasta format). '
                                   'Silva-formated taxonomies will be used if available')
    # -d / --db_dir
    group_main.add_argument('-d', '--db_dir',
                            action = 'store',
                            metavar = 'DBDIR',
                            type = str,
                            default = os.getcwd(),
                            help = 'Database output directory. '
                                   'Default is cwd')
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

    # Advanced parameters
    group_adv = parser.add_argument_group('Advanced parameters')
    # -m / --min_length
    group_adv.add_argument('-m', '--min_length',
                           action = 'store',
                           metavar = 'MINLGTH',
                           type = int,
                           default = None,
                           help = 'Ref sequences minimum length')
    # -M / --max_length
    group_adv.add_argument('-M', '--max_length',
                           action = 'store',
                           metavar = 'MAXLGTH',
                           type = int,
                           default = None,
                           help = 'Ref sequences maximum length')
    # -n / --max_consecutive_n
    group_adv.add_argument('-n', '--max_consecutive_n',
                           action = 'store',
                           metavar = 'MAXN',
                           type = int,
                           default = 5,
                           help = 'Maximum nb of consecutive Ns a sequence is allowed to have. '
                                  'Default is %(default)s. Setting it to 0 will remove all '
                                  'sequences with Ns. Ns in accepted sequences will be replaced '
                                  'by As')
    # --clustering_id_threshold
    group_adv.add_argument('--clustering_id_threshold',
                           action = 'store',
                           metavar = 'REAL',
                           type = float,
                           default = 0.95,
                           help = 'Identity threshold for clustering. '
                                  'Default is %(default)s')
    # --by_kingdoms
    group_adv.add_argument('--by_kingdom',
                           action = 'store_true',
                           help = 'Perform clustering by kingdom')
    # --kingdoms
    group_adv.add_argument('--kingdoms',
                           action = 'store',
                           metavar = 'STR',
                           type = str,
                           nargs = '+',
                           default = ['archaea', 'bacteria', 'eukaryota'],
                           help = 'Kingdoms to clusterize the DB for. '
                                  'Default is %(default)s')
    # -o / --out_db_name
    group_adv.add_argument('--out_db_name',
                           action = 'store',
                           metavar = 'OUTNAME',
                           type = str,
                           help = 'Output MATAM db name. '
                                  'Default is composed from parameters')
    # --keep_tmp
    group_adv.add_argument('--keep_tmp',
                            action = 'store_true',
                            help = 'Do not remove tmp files')
    # --debug
    group_adv.add_argument('--debug',
                            action = 'store_true',
                            help = 'Output debug infos')

    #
    args = parser.parse_args()

    # Arguments checking
    if args.clustering_id_threshold < 0 or args.clustering_id_threshold > 1:
        parser.print_help()
        raise Exception("clustering id threshold not in range [0,1]")

    # Set debug parameters
    if args.debug:
        args.verbose = True
        args.keep_tmp = True

    # Get absolute path for all arguments
    args.input_ref = os.path.abspath(args.input_ref)
    args.db_dir = os.path.abspath(args.db_dir)

    #
    return args


def print_intro(args):
    """
    Print the introduction
    """

    sys.stdout.write("""
#################################
    MATAM db pre-processing
#################################\n\n""")

    # Retrieve complete cmd line
    cmd_line = '{binpath} '.format(binpath=matam_db_prepro_bin)

    # Verbose
    if args.verbose:
        cmd_line += '--verbose '

    # Debug
    if args.debug:
        cmd_line += '--debug '

    # Performance
    cmd_line += """--cpu {cpu} --max_memory {memory} \
""".format(cpu=args.cpu,
           memory=args.max_memory)

    # Advanced parameters
    if args.min_length:
        cmd_line += '--min_length {} '.format(args.min_length)
    if args.max_length:
        cmd_line += '--max_length {} '.format(args.max_length)

    cmd_line += '--max_consecutive_n {0} '.format(args.max_consecutive_n)

    cmd_line += '--clustering_id_threshold {0} '.format(args.clustering_id_threshold)

    if args.by_kingdom:
        cmd_line += '--by_kingdom --kingdoms '
        for kingdom in args.kingdoms:
            cmd_line += '{0} '.format(kingdom)

    if args.keep_tmp:
        cmd_line += '--keep_tmp '

    if args.out_db_name:
        cmd_line += '--out_db_name {0} '.format(args.out_db_name)

    # Main parameters
    cmd_line += '--db_dir {0} '.format(args.db_dir)

    cmd_line += '--input_ref {0}'.format(args.input_ref)

    # Print cmd line
    sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))

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

    complete_ref_db_filepath = args.input_ref
    complete_ref_db_filename = os.path.basename(complete_ref_db_filepath)
    complete_ref_db_basename = os.path.splitext(complete_ref_db_filename)[0]

    try:
        if not os.path.exists(args.db_dir):
            logger.debug('mkdir {0}'.format(args.db_dir))
            os.makedirs(args.db_dir)
    except OSError:
        logger.exception('Could not create output directory {0}'.format(args.db_dir))
        raise

    complete_ref_db_taxo_filename = complete_ref_db_basename + '.taxo.tab'
    complete_ref_db_taxo_filepath = os.path.join(args.db_dir, complete_ref_db_taxo_filename)

    cleaned_complete_ref_db_basename = complete_ref_db_basename
    cleaned_complete_ref_db_filename = cleaned_complete_ref_db_basename + '.cleaned.fasta'
    cleaned_complete_ref_db_filepath = os.path.join(args.db_dir, cleaned_complete_ref_db_filename)

    clustering_id_threshold_int = int(args.clustering_id_threshold * 100)

    vsearch_output_basename = complete_ref_db_basename + '.vsearch_'
    vsearch_output_basename += '{0}'.format(clustering_id_threshold_int)
    if args.by_kingdom:
        vsearch_output_basename += '_by_kingdom'
    vsearch_output_filename = vsearch_output_basename + '.fasta'
    vsearch_output_filepath = os.path.join(args.db_dir, vsearch_output_filename)

    vsearch_centroids_basename = vsearch_output_basename + '.centroids'
    vsearch_centroids_basepath = os.path.join(args.db_dir, vsearch_output_basename)
    vsearch_centroids_filename = vsearch_centroids_basename + '.fasta'
    vsearch_centroids_filepath = os.path.join(args.db_dir, vsearch_centroids_filename)

    clustered_ref_db_basename = cleaned_complete_ref_db_basename + '_NR{0}'.format(clustering_id_threshold_int)
    if args.by_kingdom:
        clustered_ref_db_basename += '_bk'
    clustered_ref_db_basepath = os.path.join(args.db_dir, clustered_ref_db_basename)
    clustered_ref_db_filename = clustered_ref_db_basename + '.fasta'
    clustered_ref_db_filepath = os.path.join(args.db_dir, clustered_ref_db_filename)

    # This is the output MATAM db basepath to pass to matam_assembly.py
    output_ref_db_basename = clustered_ref_db_basename
    if args.out_db_name:
        output_ref_db_basename = os.path.basename(args.out_db_name)
    output_ref_db_basepath = os.path.join(args.db_dir, output_ref_db_basename)
    # This is the output MATAM db file names
    # For the complete db fasta file
    output_complete_ref_db_basename = output_ref_db_basename + '.complete'
    output_complete_ref_db_basepath = os.path.join(args.db_dir, output_complete_ref_db_basename)
    output_complete_ref_db_filename = output_complete_ref_db_basename + '.fasta'
    output_complete_ref_db_filepath = os.path.join(args.db_dir, output_complete_ref_db_filename)
    # For the complete db taxo file
    output_complete_ref_db_taxo_filename = output_complete_ref_db_basename + '.taxo.tab'
    output_complete_ref_db_taxo_filepath = os.path.join(args.db_dir, output_complete_ref_db_taxo_filename)
    # For the clustered db fasta file
    output_clustered_ref_db_basename = output_ref_db_basename + '.clustered'
    output_clustered_ref_db_basepath = os.path.join(args.db_dir, output_clustered_ref_db_basename)
    output_clustered_ref_db_filename = output_clustered_ref_db_basename + '.fasta'
    output_clustered_ref_db_filepath = os.path.join(args.db_dir, output_clustered_ref_db_filename)

    ##############################################################
    # Ref DB pre-processing (cleaning, extracting taxo, indexing)

    logger.info('Starting ref db pre-processing')

    # Extract taxo from ref DB and sort by ref id
    logger.info('Extracting taxonomies from reference DB')

    cmd_line = extract_taxo_bin + ' -i ' + complete_ref_db_filepath + ' | sort -k1,1 > '
    cmd_line += complete_ref_db_taxo_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # TO DO, maybe, one day:
    # Trim Ns from both sides

    # Convert Us in Ts
    # Option: Either filter out seq with Ns or replace Ns with random nucl
    # Option: Filter too small or too long sequences
    # Sort sequences by decreasing length
    logger.info('Cleaning reference db')

    cmd_line = 'cat ' + complete_ref_db_filepath
    cmd_line += ' | sed "/^>/!s/U/T/g" | sed "/^>/!s/u/t/g" | sed "/^>/!s/ //g"'
    cmd_line += ' | ' + replace_Ns_bin + ' -n {0} '.format(args.max_consecutive_n)
    if args.min_length or args.max_length:
        cmd_line += ' | ' + fasta_length_filter_bin
        if args.min_length:
            cmd_line += ' -m ' + str(args.min_length)
        if args.max_length:
            cmd_line += ' -M ' + str(args.max_length)
    cmd_line += ' | ' + sort_fasta_bin + ' --reverse > '
    cmd_line += cleaned_complete_ref_db_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    ####################
    # Ref DB clustering

    logger.info('Starting ref db clustering')

    # Perform clustering, either by kingdom or at once
    if args.by_kingdom:

        # Perform clustering by kingdom
        for kingdom in args.kingdoms:

            # Define by-kingdom files names + paths
            cleaned_complete_ref_db_kingdom_basename = cleaned_complete_ref_db_basename + '.' + kingdom
            cleaned_complete_ref_db_kingdom_filename = cleaned_complete_ref_db_kingdom_basename + '.fasta'
            cleaned_complete_ref_db_kingdom_filepath = os.path.join(args.db_dir, cleaned_complete_ref_db_kingdom_filename)

            vsearch_output_kingdom_basename = cleaned_complete_ref_db_kingdom_basename + '.vsearch_'
            vsearch_output_kingdom_basename += '{0}'.format(clustering_id_threshold_int)
            vsearch_output_kingdom_filename = vsearch_output_kingdom_basename + '.fasta'
            vsearch_output_kingdom_filepath = os.path.join(args.db_dir, vsearch_output_kingdom_filename)

            # Extracting kingdoms fasta files
            logger.info('Extracting sequences from {0} kingdom'.format(kingdom))

            cmd_line = fasta_name_filter_bin + ' -i ' + cleaned_complete_ref_db_filepath
            cmd_line += ' -s \' ' + kingdom + '\' > '  # !! need to be a space before the kingdom
            cmd_line += cleaned_complete_ref_db_kingdom_filepath

            logger.debug('CMD: {0}'.format(cmd_line))
            error_code += subprocess.call(cmd_line, shell=True)

            # Clutering with vsearch
            # Aplly a null-penalty value to left and right gaps (--gaopoen "20I/0E") to use vsearch as a semi-global aligner
            logger.info('Clustering {0} sequences @ {1} pct id'.format(kingdom, clustering_id_threshold_int))

            cmd_line = vsearch_bin + ' --cluster_fast'
            cmd_line += ' ' + cleaned_complete_ref_db_kingdom_filepath
            cmd_line += ' --centroids ' + vsearch_output_kingdom_filepath
            cmd_line += ' --id {0:.2f}'.format(args.clustering_id_threshold)
            cmd_line += ' --gapopen "20I/0E"'
            cmd_line += ' --threads ' + str(args.cpu)

            logger.debug('CMD: {0}'.format(cmd_line))
            if args.verbose:
                error_code += subprocess.call(cmd_line, shell=True)
            else:
                # Needed because VSearch doesnt have a verbose option
                # and output everything to stderr
                error_code += subprocess.call(cmd_line, shell=True, stdout=FNULL, stderr=FNULL)

            # Concatenate vsearch outputs
            cmd_line = 'cat ' + vsearch_output_kingdom_filepath + ' >> '
            cmd_line += vsearch_output_filepath

            logger.debug('CMD: {0}'.format(cmd_line))
            error_code += subprocess.call(cmd_line, shell=True)

            # Tag tmp files for removal
            to_rm_filepath_list.append(cleaned_complete_ref_db_kingdom_filepath)
            to_rm_filepath_list.append(vsearch_output_kingdom_filepath)

    else:
        # Clutering with vsearch
        # Aplly a null-penalty value to left and right gaps (--gaopoen "20I/0E") to use vsearch as a semi-global aligner
        logger.info('Clustering sequences @ {0} pct id'.format(clustering_id_threshold_int))
        cmd_line = vsearch_bin + ' --cluster_fast'
        cmd_line += ' ' + cleaned_complete_ref_db_filepath
        cmd_line += ' --centroids ' + vsearch_output_filepath
        cmd_line += ' --id {0:.2f}'.format(args.clustering_id_threshold)
        cmd_line += ' --gapopen "20I/0E"'
        cmd_line += ' --threads ' + str(args.cpu)

        logger.debug('CMD: {0}'.format(cmd_line))
        if args.verbose:
            error_code += subprocess.call(cmd_line, shell=True)
        else:
            # Needed because VSearch doesnt have a verbose option
            # and output everything to stderr
            error_code += subprocess.call(cmd_line, shell=True, stdout=FNULL, stderr=FNULL)
    vsearch_centroids_filepath = vsearch_output_filepath

    # Clean fasta headers
    cmd_line = clean_name_bin + ' -i ' + vsearch_centroids_filepath
    cmd_line += ' -o ' + clustered_ref_db_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Tag tmp files for removal
    to_rm_filepath_list.append(vsearch_output_filepath)
    to_rm_filepath_list.append(vsearch_centroids_filepath)

    ##########################################
    # Renaming output files as MATAM db files

    logger.info('Renaming output files as MATAM db files')

    try:
        os.rename(cleaned_complete_ref_db_filepath, output_complete_ref_db_filepath)
        os.rename(complete_ref_db_taxo_filepath, output_complete_ref_db_taxo_filepath)
        os.rename(clustered_ref_db_filepath, output_clustered_ref_db_filepath)
    except OSError:
        logger.exception('Could not rename tmp files into MATAM db files')
        raise

    ######################################################
    # SortMeRNA indexing of complete and clustered ref db

    # SortMeRNA complete ref db indexing
    logger.info('Indexing complete ref db')

    cmd_line = indexdb_bin + ' --ref ' + output_complete_ref_db_filepath
    cmd_line += ',' + output_complete_ref_db_basepath
    cmd_line += ' -m {0}'.format(args.max_memory)
    if args.verbose:
        cmd_line += ' -v '

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)
    if args.verbose:
        sys.stdout.write('\n')

    # SortMeRNA clustered ref db indexing
    logger.info('Indexing clustered ref db')

    cmd_line = indexdb_bin + ' --ref ' + output_clustered_ref_db_filepath
    cmd_line += ',' + output_clustered_ref_db_basepath
    cmd_line += ' -m {0}'.format(args.max_memory)
    if args.verbose:
        cmd_line += ' -v '

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)
    if args.verbose:
        sys.stdout.write('\n')

    ###############
    # Exit program

    sys.stdout.write('Output MATAM db: {0}\n'.format(output_ref_db_basepath))

    # Exit if everything went ok
    if not error_code:
        # Try to remove all tmp files
        # won't crash if it cannot
        if not args.keep_tmp:
            sys.stdout.write('\n')
            logger.info('Removing tmp files')
            rm_files(to_rm_filepath_list)
        #
        sys.stdout.write('\n{0} terminated with no error\n'.format(program_filename))
        exit(0)
    # Deal with errors
    else:
        sys.stdout.write('\n{0} terminated with some errors. '.format(program_filename))
        if args.verbose:
            sys.stdout.write('Check the log for additional infos\n')
        else:
            sys.stdout.write('Rerun the program using --verbose or --debug option\n')
        exit(1)
