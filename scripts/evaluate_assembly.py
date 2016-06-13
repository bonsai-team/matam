#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess


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
    
    
    # Align contigs against references with Exonerate
    
    exonerate_output_basepath = input_contigs_basepath + '.exonerate_vs_'
    exonerate_output_basepath += references_basename
    exonerate_output_filepath = exonerate_output_basepath + '.tab'
    
    cmd_line = 'exonerate --model affine:overlap --exhaustive yes --bestn 1'
    cmd_line += ' --subopt no --verbose 0 --showalignment no --showvulgar no'
    cmd_line += ' --ryo "%qi\\t%ti\\t%qab\\t%qae\\t%tab\\t%tae\\t%C\\t%em\\n"'
    cmd_line += ' ' + args.input_contigs + ' ' + args.references
    cmd_line += ' | sort -k1,1V > ' + exonerate_output_filepath
    
    sys.stdout.write('\nCMD: {0}\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    sys.stdout.write('\n')
    
    
    # Compute assembly stats
    
    assembly_stats_filepath = exonerate_output_basepath + '.assembly.stats'
    output_sam_filepath = exonerate_output_basepath + '.sam'
    
    cmd_line = 'compute_assembly_stats.py -c ' + args.input_contigs
    cmd_line += ' -r ' + args.references
    cmd_line += ' -i ' + exonerate_output_filepath
    cmd_line += ' -s ' + output_sam_filepath
    cmd_line += ' > ' + assembly_stats_filepath
    
    sys.stdout.write('CMD: {0}\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    sys.stdout.write('\n')
    
    
    
    
    
    
    
    
