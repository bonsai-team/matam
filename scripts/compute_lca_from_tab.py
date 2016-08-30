#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
import operator
from collections import Counter, defaultdict

def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in fasta_file_handle:
        if line[0] == '>':
            # Yield the last read header and sequence
            if sequence_nb:
                yield (header, ''.join(seqlines))
                del seqlines[:]
            # Get header
            header = line[1:].rstrip()
            sequence_nb += 1
        else:
            # Concatenate sequence
            seqlines.append(line.strip())
    # Yield the input file last sequence
    yield (header, ''.join(seqlines))
    # Close input file
    fasta_file_handle.close()


def read_tab_file_handle_sorted(tab_file_handle, factor_index, group_by_index, sep):
    """
    Parse a tab file (sorted by a column) and return a generator
    """
    previous_factor_id = ''
    previous_group_by_id = ''
    factor_tab_list = list()
    group_by_tab_list = list()
    # Reading tab file
    for tab in (l.strip().split(sep) for l in tab_file_handle if l.strip()):
        current_factor = tab[factor_index]
        current_group_by = tab[group_by_index]
        # Group by
        if current_group_by != previous_group_by_id:
            if previous_group_by_id:
                factor_tab_list.append(group_by_tab_list)
                group_by_tab_list = list()
        # Yield the previous factor tab list
        if current_factor != previous_factor_id:
            if previous_factor_id:
                yield factor_tab_list
                factor_tab_list = list()
        group_by_tab_list.append(tab)
        previous_factor_id = current_factor
        previous_group_by_id = current_group_by
    # Yield the last tab list
    factor_tab_list.append(group_by_tab_list)
    yield factor_tab_list
    # Close tab file
    tab_file_handle.close()


def compute_lca(factor_taxo_list, min_proportion):
    """
    Given a list of taxonomies, compute the LCA using a minimal proportion
    """
    lca_tab = list()

    # Deals with the canonic LCA right away
    if min_proportion >= 1.0:
        # List comprehension for 2 embriqued loops into a flat list
        total_taxo_tab_list = [t.split(';') for g in factor_taxo_list for t in g]
        lca_tab = os.path.commonprefix(total_taxo_tab_list)
    # Try to find a more recent LCA using at least min_proportion of the population
    else:
        lca_count = defaultdict(float)
        for group_by_taxo_list in factor_taxo_list:
            for taxo in group_by_taxo_list:
                taxo_tab = taxo.split(';')
                #~ for i in range(1, len(taxo_tab)): # Here, there was a bug, but in fact we learned that we can have better results by setting the LCA level we want in advance
                for i in range(len(taxo_tab)):
                    lca = ';'.join(taxo_tab[:i+1])
                    lca_count[lca] += 1.0/len(group_by_taxo_list)
        #~ print(len(factor_taxo_list), lca_count['Bacteria'])
        lca_count_list = [ (lca.split(';'), count) for (lca, count) in lca_count.items()]
        lca_count_list.sort(key = lambda x: (x[1], -len(x[0])), reverse=True)
        #~ print(lca_count_list)

        threshold = float(len(factor_taxo_list)) * float(min_proportion)

        for cur_lca_tab, count in lca_count_list:
            if float(count) < float(threshold):
                break
            lca_tab = cur_lca_tab

    lca = ';'.join(lca_tab)
    if not lca:
        lca = 'Root'

    return lca


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_tab', metavar='TAB',
                        type=argparse.FileType('r'), default='-',
                        help='input tab file')
    parser.add_argument('-t', '--taxo', metavar='INT',
                        type=int, required=True,
                        help='Column number containing the taxo')
    parser.add_argument('-f', '--factor', metavar='INT',
                        type=int, required=True,
                        help='Column number of the factor to compute the LCA for (eg. node or contig)')
    parser.add_argument('-g', '--group_by', metavar='INT',
                        type=int, required=True,
                        help='Column number on wich to group by (eg. reads). Taxo are weighted')
    parser.add_argument('-m', '--min_proportion', metavar='FLOAT',
                        type=float, default=1.0,
                        help='Find the most recent LCA covered by at least min_proportion '
                             'of the total population (between 0 and 1, canonic LCA by default)')
    parser.add_argument('--header', action='store_true',
                        help='Indicate there is an header line in the input file')
    parser.add_argument('-s', '--separator', metavar='CHAR',
                        type=str, default=None,
                        help='Input file separator character (whitespaces by default)')
    parser.add_argument('-o', '--output_tab', metavar='TAB',
                        type=argparse.FileType('w'), default='-',
                        help='ouput tab file')
    args = parser.parse_args()

    if args.min_proportion > 1 or args.min_proportion < 0:
        sys.stderr.write('ERROR: Min proportion has to be a real number between 0 and 1\n')
        exit(1)

    # to get indices from column numbers
    args.taxo -= 1
    args.factor -= 1
    args.group_by -= 1

    # Deals with headers
    if args.header:
        header_tab = args.input_tab.readline().strip().split(args.separator)
        factor_title = header_tab[args.factor]
        args.output_tab.write('{0};LCA\n'.format(factor_title))

    #
    for factor_tab_list in read_tab_file_handle_sorted(args.input_tab, args.factor, args.group_by, args.separator):
        factor_id = factor_tab_list[0][0][args.factor]
        #~ print(factor_id)
        factor_taxo_list = list()
        #
        for group_by_tab_list in factor_tab_list:
            group_by_taxo_list = list()
            for group_by_tab in group_by_tab_list:
                group_by_taxo_list.append(group_by_tab[args.taxo])
            factor_taxo_list.append(group_by_taxo_list)
            #~ print("\t", group_by_taxo_list)
        lca = compute_lca(factor_taxo_list, args.min_proportion)
        args.output_tab.write('{0}'.format(factor_id))
        if args.separator:
            args.output_tab.write(args.separator)
        else:
            args.output_tab.write("\t")
        args.output_tab.write('{0}\n'.format(lca))

    exit(0)

