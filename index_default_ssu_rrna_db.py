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

# Get program filename, filepath and directory
program_filepath = os.path.realpath(sys.argv[0])
program_dirpath, program_filename = os.path.split(program_filepath)
program_name, program_extension = os.path.splitext(program_filename)

# Set dependencies directories path
matam_root_dirpath = program_dirpath
matam_db_dirpath = os.path.join(matam_root_dirpath, 'db')
matam_scripts_dirpath = os.path.join(matam_root_dirpath, 'scripts')
index_ref_db_bin = os.path.join(matam_scripts_dirpath, 'index_ref_db.py')

# Set default ref db name
default_ref_db_basename = 'SILVA_128_SSURef_rdNs_NR95'
default_ref_db_basepath = os.path.join(matam_db_dirpath, default_ref_db_basename)
default_ref_db_archive_filename = default_ref_db_basename + '.tar.bz2'
default_ref_db_archive_filepath = os.path.join(matam_db_dirpath, default_ref_db_archive_filename)


def parse_arguments():
    """
    Parse the command line, and check if arguments are correct
    """
    # Initiate argument parser
    parser = argparse.ArgumentParser(description='Index default SSU rRNA DB')

    # -m / --max_memory
    parser.add_argument('-m', '--max_memory',
                        action = 'store',
                        metavar = 'MAXMEM',
                        type = int,
                        default = 10000,
                        help = 'Maximum memory to use (in MBi). '
                               'Default is %(default)s MBi')

    args = parser.parse_args()

    #
    return args


if __name__ == '__main__':

    # Set logging
    # create console handler
    ch = logging.StreamHandler()
    # set logging level
    logger.setLevel(logging.DEBUG)
    # create formatter for debug level
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # add the formatter to the console handler
    ch.setFormatter(formatter)
    # add the handler to logger
    logger.addHandler(ch)

    # Set global t0
    global_t0_wall = time.time()

    # Init global error code
    global_error_code = 0

    # Arguments parsing
    args = parse_arguments()

    #
    sys.stderr.write('\n')

    #########################
    # Get compressed archive

    logger.info('-- Get compressed archive --')

    os.chdir(matam_root_dirpath)
    logger.debug('PWD: {0}'.format(matam_root_dirpath))

    cmd_line = 'git lfs pull'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while getting the compressed archive. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    ############################
    # Extracting default ref db

    logger.info('-- Extracting default ref db --')

    os.chdir(matam_db_dirpath)
    logger.debug('PWD: {0}'.format(matam_db_dirpath))

    cmd_line = 'tar jxvf ' + default_ref_db_archive_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while extracting default ref db. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    ##########################
    # Indexing default ref db

    logger.info('-- Indexing default ref db --')

    logger.debug('PWD: {0}'.format(matam_db_dirpath))

    cmd_line = index_ref_db_bin + ' -v -i ' + default_ref_db_basepath
    cmd_line += ' --max_memory ' + str(args.max_memory)

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while indexing default ref db. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    #######
    # Exit

    logger.info('-- Completed default SSU rRNA DB indexing --')

    if global_error_code > 0:
        logger.warning('Problems might have happened during indexing. Please check log above')
    else:
        logger.debug('Indexing completed in {0:.2f} seconds'.format(time.time() - global_t0_wall))
        logger.info('Indexing went well. '
                    'Default SSU rRNA DB and its indexes can be found in '
                    'MATAM db directory: {0}*'.format(default_ref_db_basepath))

    sys.stderr.write('\n')































