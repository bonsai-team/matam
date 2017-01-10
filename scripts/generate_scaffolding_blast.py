#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict


def read_tab_file_handle_sorted(tab_file_handle, factor_index=0):
    """
    Parse a tab file (sorted by a column) and return a generator
    """
    previous_factor_id = ''
    factor_tab_list = list()
    # Reading tab file
    for line in tab_file_handle:
        l = line.strip()
        if l:
            tab = l.split()
            current_factor = tab[factor_index]
            # Yield the previous factor tab list
            if current_factor != previous_factor_id:
                if previous_factor_id:
                    yield factor_tab_list
                    factor_tab_list = list()
            factor_tab_list.append(tab)
            previous_factor_id = current_factor
    # Yield the last tab list
    yield factor_tab_list
    # Close tab file
    tab_file_handle.close()


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_tab
    parser.add_argument('-i', '--input_tab',
                        metavar='INBLAST',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input blast tab file. '
                             'Only best matches should be kept')
    # -o / --output_tab
    parser.add_argument('-o', '--output_tab',
                        metavar='OUTBLAST',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Ouput filtered blast tab file')

    args = parser.parse_args()

    # Variables initialization
    tab_list_list = list()

    # Reading blast tab file
    for tab_list in read_tab_file_handle_sorted(args.input_tab, 0):
        query_id = tab_list[0][0]
        tab_list_list.append(sorted(tab_list, key=lambda x: x[1]))

    # Sort by specificity (decreasing match number) and then
    # by decreasing alignment length. So the first seen contigs
    # are the long specific ones
    tab_list_list.sort(key=lambda x: (len(x), -int(x[0][3])))

    #
    kept_references_ids_set = set()
    fillers_tab_list_list = list()
    dealt_contigs_num = 0
    dealt_contigs_id = set()

    # Count refs
    ref_count_dict = defaultdict(int)
    for tab_list in tab_list_list:
        for tab in tab_list:
            ref_count_dict[tab[1]] += 1

    # Init a buffer for the output file
    output_tab_buffer = list()

    # Start iterating
    while (len(tab_list_list)):

        # Deal with specific contigs
        tab_list_buffer = list()

        # Purely specific contigs
        current_conservation_count = len(tab_list_list[0])
        if current_conservation_count == 1:
            for tab_list in tab_list_list:
                if len(tab_list) == 1: # Purely specific contig
                    tab = tab_list[0]
                    kept_references_ids_set.add(tab[1]) # Add the matching ref id to the set
                    # Write the alignment
                    output_tab_buffer.add('\t'.join(tab))
                    dealt_contigs_id.add(tab[0])
                    dealt_contigs_num += 1
                else:
                    tab_list_buffer.append(tab_list)
        # Multiple-choices specific contigs
        else:
            # Get most specific contig so far
            selected_tab_list = tab_list_list[0]
            #~ print(selected_tab_list)
            tab_list_buffer = tab_list_list[1:]

            selected_tab_list.sort(key=lambda t: (-ref_count_dict[t[1]], t[1]))
            selected_tab = selected_tab_list[0]

            kept_references_ids_set.add(selected_tab[1]) # Add the matching ref id to the set
            # Write the alignment
            output_tab_buffer.add('\t'.join(selected_tab))
            dealt_contigs_id.add(selected_tab[0])
            dealt_contigs_num += 1

        tab_list_list = tab_list_buffer
        tab_list_buffer = list()

        # Uses already known fillers
        fillers_tab_list_buffer = list()
        for tab_list in fillers_tab_list_list:
            #
            references_ids_set = frozenset(x[1] for x in tab_list)
            references_intersection = references_ids_set & kept_references_ids_set
            #
            tab_buffer = list()
            for tab in tab_list:
                reference_id = tab[1]
                if reference_id in references_intersection:
                    output_tab_buffer.add('\t'.join(tab))
                    dealt_contigs_id.add(tab[0])
                else:
                    # Store the alignments not already used
                    tab_buffer.append(tab)
            if tab_buffer:
                fillers_tab_list_buffer.append(tab_buffer)
            else:
                dealt_contigs_num += 1
        fillers_tab_list_list = fillers_tab_list_buffer
        fillers_tab_list_buffer = list()

        # Deal with new fillers
        for tab_list in tab_list_list:
            #
            references_ids_set = frozenset(x[1] for x in tab_list)
            references_intersection = references_ids_set & kept_references_ids_set
            #
            if len(references_intersection):  # So fillers
                tab_buffer = list()
                for tab in tab_list:
                    reference_id = tab[1]
                    if reference_id in references_intersection:
                        output_tab_buffer.add('\t'.join(tab))
                        dealt_contigs_id.add(tab[0])
                    else:
                        # Store the alignments not already used
                        tab_buffer.append(tab)
                if tab_buffer:
                    fillers_tab_list_list.append(tab_buffer)
                else:
                    dealt_contigs_num += 1
            else: # Still specific
                tab_list_buffer.append(tab_list)
        tab_list_list = tab_list_buffer
        tab_list_buffer = list()

    # Write output file
    for line in output_tab_buffer:
        args.output_tab.write('{}\n'.format(line))
