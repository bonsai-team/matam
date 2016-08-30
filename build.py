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
matam_scripts_dirpath = os.path.join(matam_root_dirpath, 'scripts')
componentsearch_dirpath = os.path.join(matam_root_dirpath, 'componentsearch')
ovgraphbuild_dirpath = os.path.join(matam_root_dirpath, 'ovgraphbuild')
sortmerna_dirpath = os.path.join(matam_root_dirpath, 'sortmerna')
sumaclust_dirpath = os.path.join(matam_root_dirpath, 'sumaclust')
bamtools_lib_dirpath = os.path.join(matam_root_dirpath, 'lib', 'bamtools')
sga_dirpath = os.path.join(matam_root_dirpath, 'sga')


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

    ########################
    # Update git submodules

    logger.info('-- Updating git submodules --')

    os.chdir(matam_root_dirpath)
    logger.debug('PWD: {0}'.format(matam_root_dirpath))

    cmd_line = 'git submodule update --init --recursive'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while updating git submodules. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    ############################
    # Compiling ComponentSearch

    logger.info('-- Compiling ComponentSearch --')

    os.chdir(componentsearch_dirpath)
    logger.debug('PWD: {0}'.format(componentsearch_dirpath))

    cmd_line = './compile.sh'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while compiling ComponentSearch. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    #########################
    # Compiling ovgraphbuild

    logger.info('-- Compiling ovgraphbuild --')

    ovgraphbuild_build_dirpath = os.path.join(ovgraphbuild_dirpath, 'build')
    makedir(ovgraphbuild_build_dirpath)
    os.chdir(ovgraphbuild_build_dirpath)
    logger.debug('PWD: {0}'.format(ovgraphbuild_build_dirpath))

    cmd_line = 'cmake .. && make'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while compiling ovgraphbuild. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    #####################
    # Building SortMeRNA

    logger.info('-- Building SortMeRNA --')

    os.chdir(sortmerna_dirpath)
    logger.debug('PWD: {0}'.format(sortmerna_dirpath))

    cmd_line = './build.sh'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while building SortMeRNA. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    ######################
    # Compiling Sumaclust

    logger.info('-- Compiling Sumaclust --')

    # Update Sumalibs to the last version
    # !! Tmp fix to get the last patches
    sumalibs_dirpath = os.path.join(sumaclust_dirpath, 'sumalibs')

    os.chdir(sumalibs_dirpath)
    logger.debug('PWD: {0}'.format(sumalibs_dirpath))

    cmd_line = 'git checkout master'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while updating Sumalibs. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    #
    os.chdir(sumaclust_dirpath)
    logger.debug('PWD: {0}'.format(sumaclust_dirpath))

    cmd_line = 'make'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while compiling Sumaclust. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    #########################
    # Compiling Bamtools lib

    logger.info('-- Compiling Bamtools lib (for SGA) --')

    bamtools_lib_build_dirpath = os.path.join(bamtools_lib_dirpath, 'build')
    makedir(bamtools_lib_build_dirpath)
    os.chdir(bamtools_lib_build_dirpath)
    logger.debug('PWD: {0}'.format(bamtools_lib_build_dirpath))

    cmd_line = 'cmake .. && make'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while compiling Bamtools lib. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    ################
    # Compiling SGA

    logger.info('-- Compiling SGA --')

    sga_src_dirpath = os.path.join(sga_dirpath, 'src')
    os.chdir(sga_src_dirpath)
    logger.debug('PWD: {0}'.format(sga_src_dirpath))

    cmd_line = './autogen.sh && '
    cmd_line += './configure --with-bamtools=' + bamtools_lib_dirpath + ' && make'

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while compiling SGA. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    ##############################
    # Creating links into bin dir

    logger.info('-- Creating links into bin dir --')

    matam_bin_dirpath = os.path.join(matam_root_dirpath, 'bin')
    makedir(matam_bin_dirpath)
    os.chdir(matam_bin_dirpath)
    logger.debug('PWD: {0}'.format(matam_bin_dirpath))

    cmd_line = 'ln -sf ' + os.path.join(matam_scripts_dirpath, 'matam_*.py')
    cmd_line += ' ' + os.path.join(matam_bin_dirpath, '.')

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code = subprocess.call(cmd_line, shell=True)

    if error_code > 0:
        logger.warning('A problem might have happened while creating links into bin dir. Check log above')

    global_error_code += error_code

    sys.stderr.write('\n')

    #######
    # Exit

    logger.info('-- MATAM building complete --')

    if global_error_code > 0:
        logger.warning('Problems might have happened during MATAM building. Please check log above')
    else:
        logger.debug('Building completed in {0:.2f} seconds'.format(time.time() - global_t0_wall))
        logger.info('MATAM building went well. '
                    'Program executables can be found in '
                    'MATAM bin directory: {0}'.format(matam_bin_dirpath))

    sys.stderr.write('\n')
