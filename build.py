#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import subprocess
import sys
import time

# Create logger
logger = logging.getLogger(__name__)

# Get program filename, filepath and directory
program_filepath = os.path.realpath(sys.argv[0])
program_dirpath, program_filename = os.path.split(program_filepath)
program_name, program_extension = os.path.splitext(program_filename)

# Set dependencies directories path
matam_root_dirpath = program_dirpath
matam_scripts_dirpath = os.path.join(matam_root_dirpath, 'scripts')
componentsearch_dirpath = os.path.join(matam_root_dirpath, 'componentsearch')
ovgraphbuild_dirpath = os.path.join(matam_root_dirpath, 'ovgraphbuild')
sortmerna_dirpath = os.path.join(matam_root_dirpath, 'sortmerna')
vsearch_dirpath = os.path.join(matam_root_dirpath, 'vsearch')
bamtools_lib_dirpath = os.path.join(matam_root_dirpath, 'lib', 'bamtools')
sga_dirpath = os.path.join(matam_root_dirpath, 'sga')
rdptools_dirpath = os.path.join(matam_root_dirpath, 'RDPTools')
kronatools_dirpath = os.path.join(matam_root_dirpath, 'Krona', 'KronaTools')


def makedir(dirpath):
    """
    Make a directory if it doesnt exist
    """
    try:
        if not os.path.exists(dirpath):
            logger.debug('mkdir {0}'.format(dirpath))
            os.makedirs(dirpath)
    except OSError:
        logger.exception('{0} directory cannot be created'.format(dirpath))
        raise


def execute_cmd(cmd_line, directory, info, warning, createdir=False):
    logger.info(info)

    if createdir:
        makedir(directory)

    os.chdir(directory)
    logger.debug('PWD: {0}'.format(directory))

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning(warning)

    sys.stderr.write('\n')

    return error_code


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

    #
    sys.stderr.write('\n')

    parser = argparse.ArgumentParser()
    valid_targets = ['build', 'clean']
    parser.add_argument(
        "target",
        nargs='?',
        default='build',
        help="Delete all files in the current directory that are \
normally created by building the program. Default is %(default)s",
        choices=valid_targets)

    args = parser.parse_args()
    print(args)

    ########################
    # Update git submodules

    if args.target == 'build':
        info = '-- Updating git submodules --'
        cmd_line = 'git submodule update --init --recursive'
        warning = 'A problem might have happened while updating git submodules. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         matam_root_dirpath, info, warning)

    ############################
    # Compiling ComponentSearch

    if args.target == 'clean':
        info = '-- Cleaning ComponentSearch --'
        cmd_line = 'make clean'
        warning = 'A problem might have happened while cleaning ComponentSearch. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         componentsearch_dirpath, info, warning)

    elif args.target == 'build':
        info = '-- Compiling ComponentSearch --'
        cmd_line = 'make'
        warning = 'A problem might have happened while compiling ComponentSearch. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         componentsearch_dirpath, info, warning)

    #########################
    # Compiling ovgraphbuild

    if args.target == 'clean':
        info = '-- Cleaning ovgraphbuild --'
        cmd_line = 'rm -rf build bin'
        warning = 'A problem might have happened while cleaning ovgraphbuild. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         ovgraphbuild_dirpath, info, warning)

    elif args.target == 'build':
        info = '-- Compiling ovgraphbuild --'
        ovgraphbuild_build_dirpath = os.path.join(
            ovgraphbuild_dirpath, 'build')
        cmd_line = 'cmake .. -G"CodeBlocks - Unix Makefiles"'
        cmd_line += '&& make'
        warning = 'A problem might have happened while compiling ovgraphbuild. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         ovgraphbuild_build_dirpath,
                                         info,
                                         warning,
                                         createdir=True)

    #####################
    # Building SortMeRNA

    if args.target == 'clean':
        info = '-- Cleaning SortMeRNA --'
        cmd_line = 'make distclean'
        warning = 'A problem might have happened while cleaning SortMeRNA. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         sortmerna_dirpath, info, warning)

    elif args.target == 'build':
        info = '-- Building SortMeRNA --'
        cmd_line = './build.sh'
        warning = 'A problem might have happened while building SortMeRNA. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         sortmerna_dirpath, info, warning)

    ######################
    # Compiling VSearch

    if args.target == 'clean':
        info = '-- Cleaning VSearch --'
        cmd_line = 'make clean'
        warning = 'A problem might have happened while cleaning VSearch. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         vsearch_dirpath, info, warning)

    elif args.target == 'build':
        info = '-- Compiling VSearch --'

        cmd_line = './autogen.sh && ./configure && make'
        warning = 'A problem might have happened while compiling VSearch. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         vsearch_dirpath, info, warning)

    #########################
    # Compiling Bamtools lib

    if args.target == 'clean':
        info = '-- Cleaning Bamtools lib (for SGA) --'
        warning = 'A problem might have happened while cleaning Bamtools lib. Check log above'
        cmd_line = 'rm -rf build bin include lib src/toolkit/bamtools_version.h'
        global_error_code += execute_cmd(cmd_line,
                                         bamtools_lib_dirpath, info, warning,)
    elif args.target == 'build':
        info = '-- Compiling Bamtools lib (for SGA) --'
        bamtools_lib_build_dirpath = os.path.join(
            bamtools_lib_dirpath, 'build')
        cmd_line = 'CC=gcc CXX=g++ cmake .. && make'
        warning = 'A problem might have happened while compiling Bamtools lib. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         bamtools_lib_build_dirpath,
                                         info,
                                         warning,
                                         createdir=True)

    ################
    # Compiling SGA

    if args.target == 'clean':
        info = '-- Cleaning SGA --'
        # "make distclean" is not enough, use git clean instead.
        # Be aware that all local changes to SGA submodule
        # will be lost
        cmd_line = 'git clean -xfd'
        warning = 'A problem might have happened while cleaning SGA. Check log above'
        global_error_code += execute_cmd(cmd_line, sga_dirpath, info, warning)
    elif args.target == 'build':
        info = '-- Compiling SGA --'
        sga_src_dirpath = os.path.join(sga_dirpath, 'src')
        cmd_line = './autogen.sh && '
        cmd_line += 'CC=gcc CXX=g++ ./configure --with-bamtools=' + \
            bamtools_lib_dirpath + ' && make'
        warning = 'A problem might have happened while compiling SGA. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         sga_src_dirpath, info, warning)

    #####################
    # Compiling RDPTools

    if args.target == 'clean':
        info = '-- Cleaning RDPTools --'
        cmd_line = 'make clean'
        warning = 'A problem might have happened while cleaning RDPTools. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         rdptools_dirpath, info, warning)

    elif args.target == 'build':
        info = '-- Compiling RDPTools --'
        cmd_line = 'make && chmod +x classifier.jar'
        warning = 'A problem might have happened while compiling RDPTools. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         rdptools_dirpath,
                                         info,
                                         warning,
                                         createdir=True)

    ########################
    # Installing KronaTools
    krona_install_dirpath = os.path.join(
        kronatools_dirpath, 'bin')

    if args.target == 'clean':
        info = '-- Cleaning Krona --'
        warning = 'A problem might have happened while cleaning Krona. Check log above'
        cmd_line = 'rm -rf %s' % krona_install_dirpath
        global_error_code += execute_cmd(cmd_line,
                                         krona_install_dirpath, info, warning)
    elif args.target == 'build':
        info = '-- Installing Krona --'
        cmd_line = './install.pl --prefix %s' % kronatools_dirpath
        warning = 'A problem might have happened while installing Krona. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         kronatools_dirpath,
                                         info,
                                         warning)

    ##############################
    # Creating links into bin dir

    if args.target == 'clean':
        pass
    elif args.target == 'build':
        info = '-- Creating links into bin dir --'
        matam_bin_dirpath = os.path.join(matam_root_dirpath, 'bin')
        cmd_line = 'ln -sf ' + \
            os.path.join(matam_scripts_dirpath, 'matam_*.py')
        cmd_line += ' ' + os.path.join(matam_root_dirpath, "index_default_ssu_rrna_db.py")
        cmd_line += ' ' + os.path.join(matam_bin_dirpath, '.')
        warning = 'A problem might have happened while creating links into bin dir. Check log above'
        global_error_code += execute_cmd(cmd_line,
                                         matam_bin_dirpath,
                                         info,
                                         warning,
                                         createdir=True)

    #######
    # Exit

    if args.target == 'clean':
        logger.info('-- MATAM cleaning complete --')
        if global_error_code > 0:
            logger.warning(
                'Problems might have happened during MATAM cleaning. Please check log above')
        else:
            logger.debug(
                'Cleaning completed in {0:.2f} seconds'.format(
                    time.time() - global_t0_wall))
            logger.info('MATAM cleaning went well. ')

    elif args.target == 'build':
        logger.info('-- MATAM building complete --')
        if global_error_code > 0:
            logger.warning(
                'Problems might have happened during MATAM building. Please check log above')
            sys.exit(global_error_code)
        else:
            logger.debug(
                'Building completed in {0:.2f} seconds'.format(
                    time.time() - global_t0_wall))
            logger.info('MATAM building went well. '
                        'Program executables can be found in '
                        'MATAM bin directory: {0}'.format(matam_bin_dirpath))

    sys.stderr.write('\n')
