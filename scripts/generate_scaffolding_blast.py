#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


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
    
    tab_list_list.sort(key=lambda x: (len(x), -int(x[0][3])))
    
    #~ print tab_list_list
    
    kept_references_ids_set = set()
    
    for tab_list in tab_list_list:
        query_id = tab_list[0][0]
        references_ids_set = frozenset(x[1] for x in tab_list)
        #~ print references_ids_set
        
        references_intersection = references_ids_set & kept_references_ids_set
        
        if not len(references_intersection):
            selected_tab = tab_list[0]
            kept_references_ids_set.add(selected_tab[1])
            args.output_tab.write('{0}\n'.format('\t'.join(selected_tab)))
        else:
            for tab in tab_list:
                reference_id = tab[1]
                if reference_id in references_intersection:
                    args.output_tab.write('{0}\n'.format('\t'.join(tab)))
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
