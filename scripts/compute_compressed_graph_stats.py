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
    # Count node arity --> Components categorization
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


def load_components_lca (components_lca_fh):
    # Init components lca dict
    components_lca_dict = dict()
    # Load components LCA
    lca_tabs = (l.split() for l in components_lca_fh if l.strip())
    for lca_t in lca_tabs:
        unitig_id = lca_t[0]
        lca = lca_t[1]
        components_lca_dict[unitig_id] = lca
    # Close components lca file
    components_lca_fh.close()
    # Debug
    #~ print sorted(components_lca_dict.items(), key=lambda x: x[0])
    #~ print
    #
    return components_lca_dict


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
    parser.add_argument('--edges_contracted', metavar='IN_FILE',
                        type=argparse.FileType('r'), required=True,
                        help='Contracted graph edges file')
    parser.add_argument('--nodes_contracted', metavar='IN_FILE',
                        type=argparse.FileType('r'), required=True,
                        help='Contracted graph nodes file')
    parser.add_argument('--components_lca', metavar='IN_FILE',
                        type=argparse.FileType('r'), required=True,
                        help='components LCA file')
    parser.add_argument('--species_taxo', metavar='IN_FILE',
                        type=argparse.FileType('r'),
                        help='Test species taxonomies')
    parser.add_argument('--read_node_component', metavar='IN_FILE',
                        type=argparse.FileType('r'),
                        help='')
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
    
    # Load components LCA
    components_lca_dict = load_components_lca(args.components_lca)
    components_num = len(components_lca_dict)
    if 'NULL' in components_lca_dict:
        # Null is not a component, it is the LCA from all the singletons
        components_num -= 1
    
    # Init stats metrics
    reads_num = 0
    nodes_num = 0
    singleton_nodes_num = 0
    hubs_num = 0
    node_strings_num = 0
    reads_in_monospecific_nodes_num = 0
    max_lca_length = 0
    lca_level_count_list = [0 for i in xrange(100)]
    lca_level_count_by_category_list = [[0 for j in xrange(100)] for i in xrange(3)]
    
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
        predicted_lca_tab = components_lca_dict[unitig_id].split(';')
        if len(predicted_lca_tab) > max_lca_length:
            max_lca_length = len(predicted_lca_tab)
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
    node_strings_num = components_num - singleton_nodes_num - hubs_num
    nodes_average_size = float(reads_num) / float(nodes_num)
    components_average_size = float(reads_num) / float(components_num)
    
    #
    if args.test_dataset:
        # Load test species taxonomies
        species_taxo_dict = load_species_taxo(args.species_taxo)
        
        #
        true_taxo_level_count_list = [[0,0] for i in xrange(max_lca_length)] # [PredictedLCA==TrueTaxo, PredictedLCA!=TrueTaxo]
        true_taxo_level_count_by_category_list = [[[0,0] for i in xrange(max_lca_length)] for j in xrange(3)]
        total_reads_num = 0
        
        #
        read_node_component_tabs = (l.split() for l in args.read_node_component if l.strip())
        for read_node_component_tab in read_node_component_tabs:
            total_reads_num += 1
            read_name = read_node_component_tab[0]
            node_id = read_node_component_tab[2]
            unitig_id = read_node_component_tab[3]
            #
            if unitig_id != 'NULL':
                specie_id = read_name[:3]
                node_category = get_node_category(node_id, nodes_arity_dict)
                read_taxo = species_taxo_dict[specie_id].split(';')
                predicted_lca = components_lca_dict[unitig_id].split(';')
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
    args.output_file.write('#Components: {0} '.format(components_num))
    args.output_file.write('(av. size = {0:.2f} reads / component)\n'.format(components_average_size))
    args.output_file.write('\t#Singletons   = {0}\n'.format(singleton_nodes_num))
    args.output_file.write('\t#Hubs         = {0}\n'.format(hubs_num))
    args.output_file.write('\t#Node strings = {0}\n\n'.format(node_strings_num))
    
    # Print taxonomical assignment statistics
    args.output_file.write('\nTaxo assignment Stats\n\n')
    args.output_file.write('Global components LCA level distribution:\n')
    args.output_file.write('Level\t#Reads\n')
    for i in xrange(max_lca_length):
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
        for i in xrange(max_lca_length):
            args.output_file.write('{0}\t{1}\n'.format(i+1, lca_level_count_by_category_list[category][i]))
        args.output_file.write('\n')
    
    #
    if args.test_dataset:
        # Print test dataset statistics
        args.output_file.write('\nTest dataset Stats\n\n')
        args.output_file.write('Global LCA prediction compatible with true taxo:\n')
        args.output_file.write('Level\t#Pos\t#Neg\t#Total\t%Pos\n')
        for i in xrange(max_lca_length):
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
            for i in xrange(max_lca_length):
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
