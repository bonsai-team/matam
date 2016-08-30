#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import defaultdict


# Compile read bases regular expression
read_bases_re = re.compile('\.|,|>|<|[ACGTNacgtn]|\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+|\*|\^.|\$')


def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


def iter_read_bases(read_bases):
    """
    """

    read_base = ''
    i = 0
    while i < len(read_bases):

        character = read_bases[i]

        # Identify indel
        if character in frozenset('-+'):
            count_str = ''
            seq = ''
            while (read_bases[i+1].isdigit()):
                i += 1
                count_str += read_bases[i]
            count = int(count_str)
            for x in range(count):
                i += 1
                seq += read_bases[i]
            #~ read_base = character + count_str + seq
            read_base = character + seq
        elif character == '^':
            i += 1
            mapping_qual = read_bases[i]
            read_base = character + mapping_qual
        elif character in frozenset('.,$ACGTNacgtn*><'):
            read_base = character
        else:
            print(read_bases)
            raise ParsingError('Character is not recognized')

        yield read_base
        i += 1


def find_called_base(ref_base, read_bases, coverage):
    """
    """

    called_base = ''

    base_dict = defaultdict(int)
    insert_dict = defaultdict(int)

    for read_base in (r.upper() for r in iter_read_bases(read_bases)):
        if read_base in frozenset('ACGTN*'):
            base_dict[read_base] += 1
        elif read_base[0] == '+':
            insert_dict[read_base[1:]] += 1

    base_sorted_list = sorted(base_dict.items(), reverse=True, key=lambda x: (x[1], x[0]))

    most_common_base = base_sorted_list[0][0]

    if most_common_base != '*':
        called_base += most_common_base

    #~ print(base_sorted_list)

    if len(insert_dict):
        insert_sorted_list = sorted(insert_dict.items(), reverse=True, key=lambda x: (x[1], x[0]))
        most_common_insert_count = insert_sorted_list[0][1]
        if most_common_insert_count >= (coverage / 2.0):
            called_base += insert_sorted_list[0][0]

    #~ print(called_base)

    return called_base


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_mpileup
    parser.add_argument('-i', '--input_mpileup',
                        metavar='INMPILEUP',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input mpileup tab file. ')
    # -o / --output_scaffolds
    parser.add_argument('-o', '--output_scaffolds',
                        metavar='OUTSCAF',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Ouput scaffolds fasta file')

    args = parser.parse_args()

    max_N_string_length = 1


    #
    previous_ref_id = ''
    previous_pos = -1
    scaffolds_list = list()
    scaffold_seq = ''

    #
    for tab in (l.split() for l in args.input_mpileup if l.strip()):
        #~ print(tab)
        ref_id = tab[0]
        position = int(tab[1])
        ref_base = tab[2].upper()
        coverage = int(tab[3])

        read_bases = tab[4]
        base_qualities = tab[5]

        gap_length = position - previous_pos

        if (ref_id != previous_ref_id) or (gap_length > max_N_string_length):
            if scaffold_seq:
                scaffolds_list.append(scaffold_seq)
                scaffold_seq = ''
        elif position > previous_pos + 1:
            for i in range(gap_length):
                scaffold_seq += 'N'

        called_base = find_called_base(ref_base, read_bases, coverage)

        scaffold_seq += called_base

        previous_ref_id = ref_id
        previous_pos = position

    if scaffold_seq:
        scaffolds_list.append(scaffold_seq)


    scaffolds_list.sort(key=lambda x: len(x))

    #~ print(scaffolds_list)

    scaffolds_to_keep_list = list()

    for i in range(len(scaffolds_list)-1):
        short_scaff = scaffolds_list[i]
        to_keep = True
        for j in range(i+1, len(scaffolds_list)):
            long_scaff = scaffolds_list[j]
            if short_scaff in long_scaff:
                to_keep = False
                break
        if to_keep:
            scaffolds_to_keep_list.append(short_scaff)
    scaffolds_to_keep_list.append(scaffolds_list[-1])


    # Write scaffolds
    scaffold_num = 0
    for scaffold_seq in scaffolds_to_keep_list:
        scaffold_num += 1
        args.output_scaffolds.write('>{0}\n{1}\n'.format(scaffold_num,
                                                         format_seq(scaffold_seq)))
