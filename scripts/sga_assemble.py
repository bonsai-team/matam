#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
import subprocess

if __name__ == '__main__':
    
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_fastq
    parser.add_argument('-i', '--input_fastq', 
                        action='store',
                        metavar='FASTQ', 
                        type=str, 
                        required=True,
                        help='Input fastq file')
    # --paired_end
    parser.add_argument('--paired_end',
                        action='store_true',
                        help='Input fastq is made of '
                             'interleaved paired-end reads')
    # -o / --output_contigs
    parser.add_argument('-o', '--output_contigs', 
                        action='store',
                        metavar='FASTA', 
                        type=str, 
                        default='contigs.fa',
                        help='Ouput contigs fasta file. '
                             'Default is ${default}')
    # --sga_bin
    parser.add_argument('--sga_bin', 
                        action='store',
                        metavar='PATH', 
                        type=str, 
                        default='sga',
                        help='SGA bin path (by default sga is searched in $PATH)')
    # --tmp_dir
    parser.add_argument('--tmp_dir', 
                        action='store',
                        metavar='DIR',
                        type=str, 
                        default='/tmp',
                        help='Tmp directory. '
                             'Default is ${default}')
    # --cpu
    parser.add_argument('--cpu', 
                        action='store',
                        metavar='INT',
                        type=int, 
                        default=3,
                        help='Max number of CPU to use. '
                             'Default is ${default}')
    #
    args = parser.parse_args()
    
    assembly_output_basename = 'assemble'
    
    # Get input and output files absolute paths
    input_filepath = os.path.realpath(args.input_fastq)
    current_working_dir = os.getcwd()
    
    if (args.output_contigs[0] == '/' or args.output_contigs[0] == '~'):
        # output file was given using an absolute path
        output_filepath = args.output_contigs
    else:
        # output file was given using a relative path
        output_filepath = current_working_dir + '/' + args.output_contigs
    
    # Mkdir tmp dir if need be
    try:
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir)
    except OSError:
        sys.stderr.write("\nERROR: {0} tmp dir cannot be created\n\n".format(args.tmp_dir))
        raise
    
    # Change cwd to tmp dir
    os.chdir(args.tmp_dir)
    
    # Cleaning last assembly contigs
    if os.path.exists(assembly_output_basename + '-contigs.fa'):
        os.remove(assembly_output_basename + '-contigs.fa')
    
    # Preprocessing
    preprocess_output = 'preprocess_output.fq'
    
    cmd_line = args.sga_bin + ' preprocess -v ' + input_filepath
    cmd_line += ' -o ' + preprocess_output
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    ## Error correction
    # Build the index that will be used for error correction
    cmd_line = args.sga_bin + ' index -a ropebwt' # ropebwt algo will only work for sequences < 200bp
    cmd_line += ' -t ' + str(args.cpu) + ' --no-reverse '
    cmd_line += preprocess_output
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Perform error correction
    kmer_cutoff = 41
    error_corrected_output_basename = 'error_corrected'
    
    cmd_line = args.sga_bin + ' correct -k ' + str(kmer_cutoff)
    cmd_line += ' --discard -x 2 -t ' + str(args.cpu)
    #~ cmd_line += ' --discard --learn -t ' + str(args.cpu)
    cmd_line += ' -o ' + error_corrected_output_basename + '.fq'
    cmd_line += ' ' + preprocess_output
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    ## Contig assembly
    # Index the corrected data
    cmd_line = args.sga_bin + ' index -a ropebwt' # ropebwt algo will only work for sequences < 200bp
    cmd_line += ' -t ' + str(args.cpu)
    cmd_line += ' ' + error_corrected_output_basename + '.fq'
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Remove exact-match duplicates and reads with low-frequency k-mers
    filtered_output = error_corrected_output_basename + '.filter.pass.fa'
    min_kmer_coverage = 2
    
    cmd_line = args.sga_bin + ' filter -x ' + str(min_kmer_coverage)
    cmd_line += ' -t ' + str(args.cpu) + ' --homopolymer-check '
    cmd_line += '--low-complexity-check -o ' + filtered_output
    cmd_line += ' ' + error_corrected_output_basename + '.fq'
    
    sys.stdout.write('\nCMD: {0}\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Merge simple, unbranched chains of vertices
    fm_merge_overlap = 55
    merged_output_basename = 'merged_output'
    
    cmd_line = args.sga_bin + ' fm-merge -m ' + str(fm_merge_overlap)
    cmd_line += ' -t ' + str(args.cpu) + ' -o ' + merged_output_basename + '.fa'
    cmd_line += ' ' + filtered_output
    
    sys.stdout.write('\nCMD: {0}\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Build an index of the merged sequences
    cmd_line = args.sga_bin + ' index -d 1000000'
    cmd_line += ' -t ' + str(args.cpu)
    cmd_line += ' ' + merged_output_basename + '.fa'
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Remove any substrings that were generated from the merge process
    cmd_line = args.sga_bin + ' rmdup'
    cmd_line += ' -t ' + str(args.cpu)
    cmd_line += ' ' + merged_output_basename + '.fa'
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Compute the structure of the string graph
    min_overlap = fm_merge_overlap
    
    cmd_line = args.sga_bin + ' overlap -m ' + str(min_overlap)
    cmd_line += ' -t ' + str(args.cpu)
    cmd_line += ' ' + merged_output_basename + '.rmdup.fa'
    
    sys.stdout.write('\nCMD: {0}\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Perform the contig assembly without bubble popping
    assembly_output_basename = 'assemble'
    assembly_contigs_filename = assembly_output_basename + '-contigs.fa'
    
    cmd_line = args.sga_bin + ' assemble -m ' + str(min_overlap)
    #~ cmd_line += ' -b 3 -d 0.03 -g 0.01 '
    cmd_line += ' -b 0'
    #~ cmd_line += ' -r 10'
    #~ cmd_line += ' --max-edges 10000 -x 10 -l 100 '
    cmd_line += ' -x 10 -l 100 '
    cmd_line += ' -o ' + assembly_output_basename
    cmd_line += ' ' + merged_output_basename + '.rmdup.asqg.gz'
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    # Scaffolding
    assembly_scaffolds_filename = assembly_output_basename + '-scaffolds.fa'
    
    if args.paired_end:
        cmd_line = 'bwa index ' + assembly_contigs_filename
        #~ cmd_line
    
    ## Final post-processing
    assembly_output_filename = assembly_output_basename
    if args.paired_end:
        assembly_output_filename += '-scaffolds.fa'
    else:
        assembly_output_filename += '-contigs.fa'
    
    #~ cmd_line = 'cp ' + assembly_output_filename + ' '
    cmd_line = 'cp ' + assembly_contigs_filename + ' '
    cmd_line += output_filepath
    
    sys.stdout.write('\nCMD: {0}\n\n'.format(cmd_line))
    subprocess.call(cmd_line, shell=True)
    
    exit(0)
