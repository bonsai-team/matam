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
    edges_contracted_fh.close()
    # Debug
    #~ print sorted(nodes_arity_dict.items(), key=lambda x: x[1], reverse=True)
    #~ print
    #
    return nodes_arity_dict


def load_unitigs_lca (unitigs_lca_fh):
    # Init unitigs lca dict
    unitigs_lca_dict = dict()
    # Load unitigs LCA
    lca_tabs = (l.split() for l in unitigs_lca_fh if l.strip())
    for lca_t in lca_tabs:
        unitig_id = lca_t[0]
        lca = lca_t[1]
        unitigs_lca_dict[unitig_id] = lca
    # Close unitigs lca file
    unitigs_lca_fh.close()
    # Debug
    #~ print sorted(unitigs_lca_dict.items(), key=lambda x: x[0])
    #~ print
    #
    return unitigs_lca_dict


def load_species_taxo (species_taxo_fh):
    # Init species taxo dict
    species_taxo_dict = dict()
    # Load species taxo
    species_taxo_tabs = (l.split() for l in species_taxo_fh if l.strip())
    for specie_taxo_tab in species_taxo_tabs:
        specie_id = specie_taxo_tab[0]
        specie_taxo = specie_taxo_tab[1]
        species_taxo_dict[specie_id] = specie_taxo
    # Close species taxo file
    species_taxo_fh.close()
    # Debug
    #~ print sorted(species_taxo_dict.items(), key=lambda x: x[0])
    #~ print
    #
    return species_taxo_dict


def get_node_category (node_id, nodes_arity_dict):
    #
    category = 2
    #
    node_arity = int(nodes_arity_dict[node_id])
    #
    if node_arity < 1:
        category = 0
    elif node_arity >= 3:
        category = 1
    #
    return category


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
    parser.add_argument('--species_taxo', metavar='IN_FILE',
                        type=argparse.FileType('r'),
                        help='Test species taxonomies')
    parser.add_argument('--read_node_unitig', metavar='IN_FILE',
                        type=argparse.FileType('r'),
                        help='')
    #~ parser.add_argument('--count_size', action='store_true', 
                        #~ help='Compute stats based on the node size rather than the node count')
    #~ parser.add_argument('--header', action='store_true', 
                        #~ help='Indicate there is an header line in the input file')
    #~ parser.add_argument('--separator', metavar='CHAR',
                        #~ type=str, default=';',
                        #~ help='Input file separator character (";" by default for csv)')
    parser.add_argument('--test_dataset', action='store_true',
                        help='Ouput additional stats when using a test dataset')
    parser.add_argument('-o', '--output_file', metavar='OUT', 
                        type=argparse.FileType('w'), default='-',
                        help='ouput file')
    args = parser.parse_args()
    
    if args.test_dataset:
        if args.species_taxo == None:
            sys.stderr.write('Test species taxonomies file needed when using a test dataset\n')
            exit(1)
    
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
    reads_in_monospecific_nodes_num = 0
    lca_level_count_list = [0 for i in xrange(7)]
    lca_level_count_by_category_list = [[0 for j in xrange(7)] for i in xrange(3)]
    
    # Remove header
    args.nodes_contracted.readline()
    # Scan contracted nodes
    nodes_tabs = (l.strip().split(';') for l in args.nodes_contracted if l.strip())
    for node_t in nodes_tabs:
        node_id = node_t[0]
        node_size = int(node_t[1])
        node_specie = node_t[2]
        unitig_id = node_t[3]
        # Get LCA
        predicted_lca_tab = unitigs_lca_dict[unitig_id].split(';')
        lca_level_count_list[len(predicted_lca_tab)-1] += node_size
        # Basic stats
        reads_num += node_size
        nodes_num += 1
        if node_specie != '?':
            reads_in_monospecific_nodes_num += node_size
        # Categorize
        node_category = get_node_category(node_id, nodes_arity_dict)
        if node_category == 0:
            singleton_nodes_num += 1
        elif node_category == 1:
            hubs_num += 1
        lca_level_count_by_category_list[node_category][len(predicted_lca_tab)-1] += node_size
    
    # Compute final stats
    node_strings_num = unitigs_num - singleton_nodes_num - hubs_num
    nodes_average_size = float(reads_num) / float(nodes_num)
    unitigs_average_size = float(reads_num) / float(unitigs_num)
    
    #
    if args.test_dataset:
        # Load test species taxonomies
        species_taxo_dict = load_species_taxo(args.species_taxo)
        
        #
        true_taxo_level_count_list = [[0,0] for i in xrange(7)] # [PredictedLCA==TrueTaxo, PredictedLCA!=TrueTaxo]
        true_taxo_level_count_by_category_list = [[[0,0] for i in xrange(7)] for j in xrange(3)]
        total_reads_num = 0
        
        #
        read_node_unitig_tabs = (l.split() for l in args.read_node_unitig if l.strip())
        for read_node_unitig_tab in read_node_unitig_tabs:
            total_reads_num += 1
            read_name = read_node_unitig_tab[0]
            node_id = read_node_unitig_tab[2]
            unitig_id = read_node_unitig_tab[3]
            #
            if unitig_id != 'NULL':
                specie_id = read_name[:3]
                node_category = get_node_category(node_id, nodes_arity_dict)
                read_taxo = species_taxo_dict[specie_id].split(';')
                predicted_lca = unitigs_lca_dict[unitig_id].split(';')
                # Compare true taxo vs. predicted LCA
                is_same = True
                for i in xrange(len(predicted_lca)):
                    if is_same:
                        is_same = (predicted_lca[i]==read_taxo[i])
                    if is_same:
                        column = 0
                    else:
                        column = 1
                    true_taxo_level_count_list[i][column] += 1
                    true_taxo_level_count_by_category_list[node_category][i][column] += 1
    
    # Start outputing infos 
    args.output_file.write('Graph: {0}\n\n'.format(args.edges_contracted.name))
    
    # Print descriptive statistics
    args.output_file.write('\nDescriptive Stats\n\n')
    args.output_file.write('#Reads:   {0} mapped reads'.format(reads_num))
    if args.test_dataset:
        args.output_file.write(' / {0} total reads '.format(total_reads_num))
        mapped_reads_percent = reads_num * 100.0 / total_reads_num
        args.output_file.write('({0:.2f}%)\n'.format(mapped_reads_percent))
        reads_in_monospecific_nodes_percent = reads_in_monospecific_nodes_num * 100.0 / reads_num
        args.output_file.write('#Reads in monospecific nodes: {0}'.format(reads_in_monospecific_nodes_num))
        args.output_file.write(' ({0:.2f}% of all reads in the graph)'.format(reads_in_monospecific_nodes_percent))
    args.output_file.write('\n')
    args.output_file.write('#Nodes:   {0} '.format(nodes_num))
    args.output_file.write('(av. size = {0:.2f} reads / node)\n'.format(nodes_average_size))
    args.output_file.write('#Unitigs: {0} '.format(unitigs_num))
    args.output_file.write('(av. size = {0:.2f} reads / unitig)\n'.format(unitigs_average_size))
    args.output_file.write('\t#Singletons   = {0}\n'.format(singleton_nodes_num))
    args.output_file.write('\t#Hubs         = {0}\n'.format(hubs_num))
    args.output_file.write('\t#Node strings = {0}\n\n'.format(node_strings_num))
    
    # Print taxonomical assignment statistics
    args.output_file.write('\nTaxo assignment Stats\n\n')
    args.output_file.write('Global unitigs LCA level distribution:\n')
    args.output_file.write('Level\t#Reads\n')
    for i in xrange(7):
        args.output_file.write('{0}\t{1}\n'.format(i+1, lca_level_count_list[i]))
    args.output_file.write('\n')
    
    # And now by categories
    for category in xrange(3):
        if category == 0:
            args.output_file.write('Singletons ')
        elif category == 1:
            args.output_file.write('Hubs ')
        else:
            args.output_file.write('Node strings ')
        args.output_file.write('LCA level distribution:\n')
        args.output_file.write('Level\t#Reads\n')
        for i in xrange(7):
            args.output_file.write('{0}\t{1}\n'.format(i+1, lca_level_count_by_category_list[category][i]))
        args.output_file.write('\n')
    
    #
    if args.test_dataset:
        # Print test dataset statistics
        args.output_file.write('\nTest dataset Stats\n\n')
        args.output_file.write('Global LCA prediction compatible with true taxo:\n')
        args.output_file.write('Level\t#Pos\t#Neg\t#Total\t%Pos\n')
        for i in xrange(7):
            args.output_file.write('{}\t{}\t{}\t{}'.format(i+1, 
                                                           true_taxo_level_count_list[i][0], 
                                                           true_taxo_level_count_list[i][1], 
                                                           sum(true_taxo_level_count_list[i])))
            if (sum(true_taxo_level_count_list[i])):
                true_assignment_percent = true_taxo_level_count_list[i][0] * 100.0 / sum(true_taxo_level_count_list[i])
                args.output_file.write('\t{0:.2f}%'.format(true_assignment_percent))
            else:
                args.output_file.write('\tNA')
            args.output_file.write('\n')
        
        # And now by categories
        for category in xrange(3):
            if category == 0:
                args.output_file.write('Singletons ')
            elif category == 1:
                args.output_file.write('Hubs ')
            else:
                args.output_file.write('Node strings ')
            args.output_file.write('LCA prediction compatible with true taxo:\n')
            args.output_file.write('Level\t#Pos\t#Neg\t#Total\t%Pos\n')
            for i in xrange(7):
                args.output_file.write('{}\t{}\t{}\t{}'.format(i+1, 
                                                               true_taxo_level_count_by_category_list[category][i][0], 
                                                               true_taxo_level_count_by_category_list[category][i][1], 
                                                               sum(true_taxo_level_count_by_category_list[category][i])))
                if (sum(true_taxo_level_count_by_category_list[category][i])):
                    true_assignment_percent = true_taxo_level_count_by_category_list[category][i][0] * 100.0 / sum(true_taxo_level_count_by_category_list[category][i])
                    args.output_file.write('\t{0:.2f}%'.format(true_assignment_percent))
                else:
                    args.output_file.write('\tNA')
                args.output_file.write('\n')
            args.output_file.write('\n')
    
    exit(0)
