#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
from collections import defaultdict


def load_nodes_arity (edges_contracted_fh):
    # Init node arity dict
    nodes_arity_dict = defaultdict(int)
    # Remove header
    args.edges_contracted.readline()
    # Count node arity --> Unitigs categorization
    edges_tabs = (l.strip().split(';') for l in edges_contracted_fh if l.strip())
    for edge_t in edges_tabs:
        nodes_arity_dict[edge_t[0]] += 1 # Source
        nodes_arity_dict[edge_t[1]] += 1 # Target
    # Close contracted edges file
    args.edges_contracted.close()
    # Debug
    #~ print sorted(nodes_arity_dict.items(), key=lambda x: x[1], reverse=True)
    #~ print
    #
    return nodes_arity_dict


def load_unitigs_lca (unitigs_lca_fh):
    # Init unitigs lca dict
    unitigs_lca_dict = dict()
    # Load unitigs LCA
    lca_tabs = (l.split() for l in args.unitigs_lca if l.strip())
    for lca_t in lca_tabs:
        unitig_id = int(lca_t[0])
        lca = lca_t[1]
        unitigs_lca_dict[unitig_id] = lca
    # Close unitigs lca file
    args.unitigs_lca.close()
    # Debug
    #~ print sorted(unitigs_lca_dict.items(), key=lambda x: x[0])
    #~ print
    #
    return unitigs_lca_dict



if __name__ == '__main__':
    
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    #~ parser.add_argument('-i', '--input_file', metavar='IN', 
                        #~ type=argparse.FileType('r', 0), default='-',
                        #~ help='Input tabulated file')
    parser.add_argument('--edges_contracted', metavar='IN_FILE',
                        type=argparse.FileType('r'), required=True,
                        help='Contracted graph edges file')
    parser.add_argument('--nodes_contracted', metavar='IN_FILE',
                        type=argparse.FileType('r'), required=True,
                        help='Contracted graph nodes file')
    parser.add_argument('--unitigs_lca', metavar='IN_FILE',
                        type=argparse.FileType('r'), required=True,
                        help='Unitigs LCA file')
    #~ parser.add_argument('-l', '--predicted_lca', metavar='INT',
                        #~ type=int, required=True,
                        #~ help='Predicted LCA column number')
    #~ parser.add_argument('-t', '--true_taxo', metavar='INT', 
                        #~ type=int, required=True,
                        #~ help='True taxonomy column number')
    #~ parser.add_argument('-s', '--size', metavar='INT', 
                        #~ type=int, required=True,
                        #~ help='Node size column number')
    #~ parser.add_argument('--count_size', action='store_true', 
                        #~ help='Compute stats based on the node size rather than the node count')
    #~ parser.add_argument('--header', action='store_true', 
                        #~ help='Indicate there is an header line in the input file')
    #~ parser.add_argument('--separator', metavar='CHAR',
                        #~ type=str, default=';',
                        #~ help='Input file separator character (";" by default for csv)')
    parser.add_argument('-o', '--output_file', metavar='OUT', 
                        type=argparse.FileType('w'), default='-',
                        help='ouput file')
    args = parser.parse_args()
    
    
    # Load nodes arity dict
    nodes_arity_dict = load_nodes_arity(args.edges_contracted)
    
    # Load unitigs LCA
    unitigs_lca_dict = load_unitigs_lca(args.unitigs_lca)
    unitigs_num = len(unitigs_lca_dict)
    
    
    # Init stats metrics
    reads_num = 0
    nodes_num = 0
    singleton_nodes_num = 0
    hubs_num = 0
    node_strings_num = 0
    lca_level_count_list = [0 for i in xrange(7)]
    lca_level_count_by_category_list = [[0 for j in xrange(7)] for i in xrange(3)]
        
    # Remove header
    args.nodes_contracted.readline()
    # Scan contracted nodes
    nodes_tabs = (l.strip().split(';') for l in args.nodes_contracted if l.strip())
    for node_t in nodes_tabs:
        node_id = node_t[0]
        node_size = int(node_t[1])
        unitig_id = int(node_t[3])
        # Get LCA
        predicted_lca_tab = unitigs_lca_dict[unitig_id].split(';')
        lca_level_count_list[len(predicted_lca_tab)-1] += node_size
        # Basic stats
        reads_num += node_size
        nodes_num += 1
        # Categorize
        node_arity = int(nodes_arity_dict[node_id])
        category = 2
        if node_arity < 1:
            singleton_nodes_num += 1
            category = 0
        elif node_arity >= 3:
            hubs_num += 1
            category = 1
        lca_level_count_by_category_list[category][len(predicted_lca_tab)-1] += node_size
    
    node_strings_num = unitigs_num - singleton_nodes_num - hubs_num
    
    
    # Compute final stats
    nodes_average_size = float(reads_num) / float(nodes_num)
    unitigs_average_size = float(reads_num) / float(unitigs_num)
    
    
    # Start outputing infos 
    args.output_file.write('Graph: {0}\n\n'.format(args.edges_contracted.name))
    
    # Print descriptive statistics
    args.output_file.write('\nDescriptive Stats\n\n')
    args.output_file.write('#Reads:   {0}\n'.format(reads_num))
    args.output_file.write('#Nodes:   {0} '.format(nodes_num))
    args.output_file.write('(av. size = {0:.2f} reads)\n'.format(nodes_average_size))
    args.output_file.write('#Unitigs: {0} '.format(unitigs_num))
    args.output_file.write('(av. size = {0:.2f} reads)\n'.format(unitigs_average_size))
    args.output_file.write('\t#Singletons   = {0}\n'.format(singleton_nodes_num))
    args.output_file.write('\t#Hubs         = {0}\n'.format(hubs_num))
    args.output_file.write('\t#Node strings = {0}\n\n'.format(node_strings_num))
    
    # Print taxonomical assignment statistics
    args.output_file.write('\nTaxo assignment Stats\n\n')
    #~ args.output_file.write('Unitigs LCA level distribution:\n')
    #~ args.output_file.write('Level\t#Reads\n')
    #~ for i in xrange(7):
        #~ args.output_file.write('{0}\t{1}\n'.format(i+1, lca_level_count_list[i]))
    #~ args.output_file.write('\n')
    
    # And now by categories
    for j in xrange(3):
        if j == 0:
            args.output_file.write('Singletons ')
        elif j == 1:
            args.output_file.write('Hubs ')
        else:
            args.output_file.write('Node strings ')
        args.output_file.write('LCA level distribution:\n')
        args.output_file.write('Level\t#Reads\n')
        for i in xrange(7):
            args.output_file.write('{0}\t{1}\n'.format(i+1, lca_level_count_by_category_list[j][i]))
        args.output_file.write('\n')

    
    #~ # to get indices from column numbers
    #~ args.predicted_lca -= 1
    #~ args.true_taxo -= 1
    #~ args.size -= 1
    #~ 
    #~ # Deals with headers
    #~ if args.header:
        #~ args.input_file.readline()
    #~ 
    #~ level_count_list = [0 for i in xrange(7)]
    #~ stats_level_list = [[0,0] for i in xrange(7)] # [PredictedLCA==TrueTaxo, PredictedLCA!=TrueTaxo]
    #~ 
    #~ # Count 
    #~ for line in args.input_file:
        #~ tab = line.strip().split(args.separator)
        #~ predicted_lca = tab[args.predicted_lca].split(',')
        #~ true_taxo = tab[args.true_taxo].split(',')
        #~ size = int(tab[args.size])
        #~ 
        #~ # Predicted LCA level count
        #~ if args.count_size:
            #~ level_count_list[len(predicted_lca)-1] += size
        #~ else:
            #~ level_count_list[len(predicted_lca)-1] += 1
        #~ 
        #~ # Is LCA prediction compatible with true taxonomy 
        #~ if true_taxo[0] != 'NULL':
            #~ is_same = True
            #~ for i in xrange(len(predicted_lca)):
                #~ if is_same:
                    #~ is_same = (predicted_lca[i]==true_taxo[i])
                #~ 
                #~ if is_same:
                    #~ column = 0
                #~ else:
                    #~ column = 1
                #~ 
                #~ if args.count_size:
                    #~ stats_level_list[i][column] += size
                #~ else:
                    #~ stats_level_list[i][column] += 1
    #~ 
    #~ # Output stats
    #~ total = sum(level_count_list)
    #~ monospecific_total = sum(stats_level_list[0])
    #~ monospecific_percent = monospecific_total * 100.0 / total
    #~ 
    #~ if args.count_size:
        #~ args.output_file.write('# Total size = {0}\n'.format(total))
    #~ else:
        #~ args.output_file.write('# Nodes num = {0}\n'.format(total))
    #~ args.output_file.write('# {0:.2f}% in mono-specific nodes\n\n'.format(monospecific_percent))
    #~ 
    #~ args.output_file.write('Level Stats on predicted LCA:\nLevel\tNumber\n')
    #~ for i in xrange(7):
        #~ args.output_file.write('{0}\t{1}\n'.format(i+1, level_count_list[i]))
        #~ 
    #~ args.output_file.write('\nLCA prediction compatible with true taxo:\n')
    #~ args.output_file.write('Level\t#Pos\t#Neg\t#Total\t%Pos\n')
    #~ 
    #~ for i in xrange(7):
        #~ args.output_file.write('{}\t{}\t{}\t{}'.format(i+1, stats_level_list[i][0], stats_level_list[i][1], sum(stats_level_list[i])))
        #~ if (sum(stats_level_list[i])):
            #~ args.output_file.write('\t{0:.2f}%'.format(stats_level_list[i][0] * 100.0 / sum(stats_level_list[i])))
        #~ else:
            #~ args.output_file.write('\tNA')
        #~ args.output_file.write('\n')
    #~ 
    exit(0)
