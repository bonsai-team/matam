#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from collections import defaultdict

def group_by_species_name(ref_db_taxonomies_list):
    ref_db_taxonomies_buffer = list()
    previous_species_name = str()
    for ref_db_taxo_tuple in ref_db_taxonomies_list:
        species_name = ref_db_taxo_tuple[0][-1]
        if previous_species_name and previous_species_name != species_name:
            yield(previous_species_name, ref_db_taxonomies_buffer)
            ref_db_taxonomies_buffer.clear()
        ref_db_taxonomies_buffer.append(ref_db_taxo_tuple)
        previous_species_name = species_name
    yield(previous_species_name, ref_db_taxonomies_buffer)


def get_binominal_name(full_name):
    binominal_name = str()
    if ' sp. ' in full_name:
        binominal_name = ' '.join(full_name.split()[0:3])
    else:
        binominal_name = ' '.join(full_name.split()[0:2])
    return binominal_name


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_species
    parser.add_argument('-i', '--input_species',
                        metavar='FILE',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input file. One species name by line.')
    # -r / --ref_db_taxonomies
    parser.add_argument('-r', '--ref_db_taxonomies',
                        metavar='FILE',
                        type=argparse.FileType('r'),
                        required=True,
                        help='Taxonomy reference database. Format is "count[TAB]taxonomy"')
    # -o / --output_taxonomies
    parser.add_argument('-o', '--output_taxonomies',
                        metavar='FILE',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Ouput full taxonomies for each species.')

    args = parser.parse_args()

    #
    ref_db_taxonomies_list = list()

    for ref_db_taxo in args.ref_db_taxonomies:
        tab = ref_db_taxo.split()
        count = int(tab[0])
        taxo_tab = ' '.join(tab[1:]).split(';')
        # print('{}\t{}'.format(taxo_tab, count))
        ref_db_taxonomies_list.append((taxo_tab, count))

    ref_db_taxonomies_list.sort(key=lambda x: (x[0][-1],-x[1]))
    full_name_taxo_dict = dict()

    # for ref_db_taxo_tuple in ref_db_taxonomies_list:
        # print('{}\t{}'.format(ref_db_taxo_tuple[0], ref_db_taxo_tuple[1]))

    for species_name, ref_db_taxonomies_buffer in group_by_species_name(ref_db_taxonomies_list):
        # print('{}\t{}'.format(species_name, ref_db_taxonomies_buffer))
        full_name_taxo_dict[species_name] = list(ref_db_taxonomies_buffer[0][0])

    # print(full_name_taxo_dict['Escherichia coli'])

    binominal_name_taxo_dict = dict()

    for ref_db_taxo_tuple in ref_db_taxonomies_list:
        ref_db_taxo_tuple[0][-1] = get_binominal_name(ref_db_taxo_tuple[0][-1])

    # for ref_db_taxo_tuple in ref_db_taxonomies_list:
    #     print('{}\t{}'.format(ref_db_taxo_tuple[0], ref_db_taxo_tuple[1]))

    ref_db_taxonomies_list.sort(key=lambda x: (x[0][-1],-x[1]))

    for species_name, ref_db_taxonomies_buffer in group_by_species_name(ref_db_taxonomies_list):
        count_dict = defaultdict(int)
        # print(ref_db_taxonomies_buffer)
        for ref_db_taxo_tuple in ref_db_taxonomies_buffer:
            count_dict[';'.join(ref_db_taxo_tuple[0])] += ref_db_taxo_tuple[1]
        # print(count_dict)
        # print(sorted(count_dict.items(), key=lambda x: x[1], reverse=True))
        binominal_name_taxo_dict[species_name] = list(sorted(count_dict.items(), key=lambda x: x[1], reverse=True)[0][0].split(';'))

    # print(binominal_name_taxo_dict['Escherichia coli'])
    # print(binominal_name_taxo_dict)

    genus_name_taxo_dict = dict()

    for ref_db_taxo_tuple in ref_db_taxonomies_list:
        ref_db_taxo_tuple[0][-1] = ref_db_taxo_tuple[0][-1].split()[0]

    ref_db_taxonomies_list.sort(key=lambda x: (x[0][-1],-x[1]))

    for species_name, ref_db_taxonomies_buffer in group_by_species_name(ref_db_taxonomies_list):
        count_dict = defaultdict(int)
        # print(ref_db_taxonomies_buffer)
        for ref_db_taxo_tuple in ref_db_taxonomies_buffer:
            count_dict[';'.join(ref_db_taxo_tuple[0])] += ref_db_taxo_tuple[1]
        # print(count_dict)
        # print(sorted(count_dict.items(), key=lambda x: x[1], reverse=True))
        genus_name_taxo_dict[species_name] = list(sorted(count_dict.items(), key=lambda x: x[1], reverse=True)[0][0].split(';'))

    # print(full_name_taxo_dict['Escherichia coli'])

    # Get taxonomies

    for full_name in args.input_species:
        taxo = list()
        full_name = full_name.strip()
        if full_name in full_name_taxo_dict:
            taxo = full_name_taxo_dict[full_name]
        else:
            binominal_name = get_binominal_name(full_name)
            if binominal_name in binominal_name_taxo_dict:
                taxo = binominal_name_taxo_dict[binominal_name]
            else:
                genus_name = full_name.split()[0]
                if genus_name in genus_name_taxo_dict:
                    taxo = genus_name_taxo_dict[genus_name]
        if taxo:
            args.output_taxonomies.write('{}\t{}\n'.format(full_name, ';'.join(taxo)))
        else:
            print('WARNING: {} could not be found in the taxonomy'.format(full_name), file=sys.stderr)