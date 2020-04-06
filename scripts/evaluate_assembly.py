#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
from binary_utils import Binary

# Get script dir absolute path
program_filename = os.path.basename(sys.argv[0])
program_filepath = os.path.realpath(sys.argv[0])
script_dir = os.path.dirname(program_filepath)

# Get dependencies bin
exonerate_bin = Binary.assert_which('exonerate')
exonerate_to_sam_bin = os.path.join(script_dir, 'exonerate_to_sam.py')
compute_assembly_stats_bin = os.path.join(script_dir, 'compute_assembly_stats.py')

if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_contigs
    parser.add_argument('-i', '--input_contigs',
                        metavar='CONTIGS',
                        type=str,
                        required=True,
                        help='Assembly contigs fasta file')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar='REF',
                        type=str,
                        required=True,
                        help='References fasta file')

    args = parser.parse_args()

    # Get absolute path for all arguments
    args.input_contigs = os.path.abspath(args.input_contigs)
    args.references = os.path.abspath(args.references)

    #
    input_contigs_basepath = '.'.join(args.input_contigs.split('.')[:-1])
    references_filename = args.references.split('/')[-1]
    references_basename = '.'.join(references_filename.split('.')[:-1])

    # Exit if the fasta file is empty
    if os.path.getsize(args.input_contigs) == 0:
        sys.stderr.write('\n{0} fasta file is empty\n\n'.format(args.input_contigs))
        exit(0)

    # Align contigs against references with Exonerate
    exonerate_output_basepath = input_contigs_basepath + '.exonerate_vs_'
    exonerate_output_basepath += references_basename
    exonerate_output_filepath = exonerate_output_basepath + '.tab'

    cmd_line = exonerate_bin + ' --model affine:overlap --exhaustive yes --bestn 1'
    cmd_line += ' --subopt no --verbose 0 --showalignment no --showvulgar no'
    cmd_line += ' --ryo "%qi\\t%ti\\t%qab\\t%qae\\t%tab\\t%tae\\t%C\\n"'
    cmd_line += ' ' + args.input_contigs + ' ' + args.references
    cmd_line += ' | sort -k1,1V > ' + exonerate_output_filepath

    sys.stderr.write('\nCMD: {0}\n'.format(cmd_line))
    subprocess.check_call(cmd_line, shell=True)

    # Convert Exonerate tab output to sam file
    exonerate_sam_filepath = exonerate_output_basepath + '.sam'

    cmd_line = exonerate_to_sam_bin + ' -i ' + exonerate_output_filepath
    cmd_line += ' -o ' + exonerate_sam_filepath
    cmd_line += ' -q ' + args.input_contigs

    sys.stderr.write('\nCMD: {0}\n'.format(cmd_line))
    subprocess.check_call(cmd_line, shell=True)
    sys.stderr.write('\n')

    # Compute assembly stats
    assembly_stats_filepath = exonerate_output_basepath + '.assembly.stats'

    cmd_line = compute_assembly_stats_bin + ' -r ' + args.references
    cmd_line += ' -i ' + exonerate_sam_filepath
    cmd_line += ' > ' + assembly_stats_filepath

    sys.stderr.write('CMD: {0}\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    sys.stderr.write('\n')
