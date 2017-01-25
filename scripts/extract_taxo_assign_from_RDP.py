#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --rdp_file
    parser.add_argument('-i', '--rdp_file',
                        metavar = 'PATH',
                        action = 'store',
                        type = argparse.FileType('r'),
                        default = '-',
                        help = 'RDP classification file')
    # -o / --output_tab
    parser.add_argument('-o', '--output_tab',
                        metavar = 'PATH',
                        action = 'store',
                        type = argparse.FileType('w'),
                        default = '-',
                        help = 'Taxonomies output tab file')
    # -t / --confidence_threshold
    parser.add_argument('-t', '--confidence_threshold',
                        metavar = 'INT',
                        action = 'store',
                        type = int,
                        default = 50,
                        help = 'Minimal confidence threshold. '
                               'Default is %(default)s')

    args = parser.parse_args()

    # Deal with header
    for i in range(7):
        args.rdp_file.readline()
    
    for tab in (l.strip().split(';') for l in args.rdp_file if l.strip()):
        taxo_tuple_list = list(zip(tab[0::2], tab[1::2]))
        taxon_list = list()
        for taxon_str, confidence_str in taxo_tuple_list[1:]:
            taxon = taxon_str.replace('"', '').replace(' ', '_')
            confidence = int(confidence_str.split()[0])
            if confidence >= args.confidence_threshold:
                taxon_list.append(taxon)
            else:
                break
        if taxon_list:
            args.output_tab.write('{}\n'.format('\t'.join(taxon_list)))
            
            
    
