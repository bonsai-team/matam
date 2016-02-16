#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re

if __name__ == '__main__':
    
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_file', metavar='IN', 
                        type=argparse.FileType('r', 0), default='-',
                        help='input file')
    parser.add_argument('-l', '--predicted_lca', metavar='INT',
                        type=int, required=True,
                        help='Column number containing the predicted LCA')
    parser.add_argument('-t', '--true_taxo', metavar='INT', 
                        type=int, required=True,
                        help='Column number containing the true taxonomy')
    parser.add_argument('--header', action='store_true', 
                        help='Indicate there is an header line in the input file')
    parser.add_argument('-s', '--separator', metavar='CHAR',
                        type=str, default=';',
                        help='Input file separator character (";" by default for csv)')
    parser.add_argument('-o', '--output_file', metavar='OUT', 
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput file')
    args = parser.parse_args()
    
    # to get indices from column numbers
    args.predicted_lca -= 1
    args.true_taxo -= 1
    
    # Deals with headers
    if args.header:
        args.input_file.readline()
    
    level_count_list = [0 for i in xrange(7)]
    stats_level_list = [[0,0] for i in xrange(7)] # [PredictedLCA==TrueTaxo, PredictedLCA!=TrueTaxo]
    
    # Count 
    for line in args.input_file:
        tab = line.strip().split(args.separator)
        predicted_lca = tab[args.predicted_lca].split(',')
        true_taxo = tab[args.true_taxo].split(',')
        
        # Predicted LCA level count
        level_count_list[len(predicted_lca)-1] += 1
        
        # Is LCA prediction compatible with true taxonomy 
        if true_taxo[0] != 'NULL':
            is_same = True
            for i in xrange(len(predicted_lca)):
                if is_same:
                    is_same = (predicted_lca[i]==true_taxo[i])
                if is_same:
                    stats_level_list[i][0] += 1
                else:
                    stats_level_list[i][1] += 1
    
    # Output stats
    args.output_file.write('#Nodes={0}\n\n'.format(sum(level_count_list)))
    
    args.output_file.write('Level Stats on predicted LCA:\nLevel\tNumber\n')
    for i in xrange(7):
        args.output_file.write('{0}\t{1}\n'.format(i+1, level_count_list[i]))
        
    args.output_file.write('\nLCA prediction compatible with true taxo:\n')
    args.output_file.write('Level\t#Pos\t#Neg\t#Total\t%Pos\n')
    
    for i in xrange(7):
        args.output_file.write('{}\t{}\t{}\t{}'.format(i+1, stats_level_list[i][0], stats_level_list[i][1], sum(stats_level_list[i])))
        if (sum(stats_level_list[i])):
            args.output_file.write('\t{0:.2f}%'.format(stats_level_list[i][0] * 100.0 / sum(stats_level_list[i])))
        else:
            args.output_file.write('\tNA')
        args.output_file.write('\n')
    
    return 0
