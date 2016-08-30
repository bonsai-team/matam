#!/usr/bin/env python3
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
                        help='Input tabulated file')
    parser.add_argument('-l', '--predicted_lca', metavar='INT',
                        type=int, required=True,
                        help='Predicted LCA column number')
    parser.add_argument('-t', '--true_taxo', metavar='INT',
                        type=int, required=True,
                        help='True taxonomy column number')
    parser.add_argument('-s', '--size', metavar='INT',
                        type=int, required=True,
                        help='Node size column number')
    parser.add_argument('--count_size', action='store_true',
                        help='Compute stats based on the node size rather than the node count')
    parser.add_argument('--header', action='store_true',
                        help='Indicate there is an header line in the input file')
    parser.add_argument('--separator', metavar='CHAR',
                        type=str, default=';',
                        help='Input file separator character (";" by default for csv)')
    parser.add_argument('-o', '--output_file', metavar='OUT',
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput file')
    args = parser.parse_args()

    # to get indices from column numbers
    args.predicted_lca -= 1
    args.true_taxo -= 1
    args.size -= 1

    # Deals with headers
    if args.header:
        args.input_file.readline()

    level_count_list = [0 for i in range(7)]
    stats_level_list = [[0,0] for i in range(7)] # [PredictedLCA==TrueTaxo, PredictedLCA!=TrueTaxo]

    # Count
    for line in args.input_file:
        tab = line.strip().split(args.separator)
        predicted_lca = tab[args.predicted_lca].split(',')
        true_taxo = tab[args.true_taxo].split(',')
        size = int(tab[args.size])

        # Predicted LCA level count
        if args.count_size:
            level_count_list[len(predicted_lca)-1] += size
        else:
            level_count_list[len(predicted_lca)-1] += 1

        # Is LCA prediction compatible with true taxonomy
        if true_taxo[0] != 'NULL':
            is_same = True
            for i in range(len(predicted_lca)):
                if is_same:
                    is_same = (predicted_lca[i]==true_taxo[i])

                if is_same:
                    column = 0
                else:
                    column = 1

                if args.count_size:
                    stats_level_list[i][column] += size
                else:
                    stats_level_list[i][column] += 1

    # Output stats
    total = sum(level_count_list)
    monospecific_total = sum(stats_level_list[0])
    monospecific_percent = monospecific_total * 100.0 / total

    if args.count_size:
        args.output_file.write('# Total size = {0}\n'.format(total))
    else:
        args.output_file.write('# Nodes num = {0}\n'.format(total))
    args.output_file.write('# {0:.2f}% in mono-specific nodes\n\n'.format(monospecific_percent))

    args.output_file.write('Level Stats on predicted LCA:\nLevel\tNumber\n')
    for i in range(7):
        args.output_file.write('{0}\t{1}\n'.format(i+1, level_count_list[i]))

    args.output_file.write('\nLCA prediction compatible with true taxo:\n')
    args.output_file.write('Level\t#Pos\t#Neg\t#Total\t%Pos\n')

    for i in range(7):
        args.output_file.write('{}\t{}\t{}\t{}'.format(i+1, stats_level_list[i][0], stats_level_list[i][1], sum(stats_level_list[i])))
        if (sum(stats_level_list[i])):
            args.output_file.write('\t{0:.2f}%'.format(stats_level_list[i][0] * 100.0 / sum(stats_level_list[i])))
        else:
            args.output_file.write('\tNA')
        args.output_file.write('\n')

    exit(0)
