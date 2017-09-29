#!/usr/bin/env python3

import sys
import os
import subprocess
import logging
import re
import argparse
import collections
import operator
import itertools

from fasta_clean_name import read_fasta_file_handle
from compute_abundance import get_abundance_from_fasta

from rdp import read_rpd_file
from rdp import get_lineage

logger = logging.getLogger(__name__)


class SampleCollection():

    def __init__(self, samples_path):
        self.float_precision = 4
        self.samples = samples_path
        self.comparaison_table = [('Taxonomy', 'SampleID', 'SequenceID', 'Abundance', 'Normalized_Abundance(%)')]

        if not samples_path:
            logger.fatal('The sample collection is empty')
            sys.exit('Empty collection')
        
        self.check_path_validity()
        self.build_tables()


    def write_comparaison_table(self, out_handler):
        for row in self.comparaison_table:
            str_row = [str(v) for v in row]
            print('\t'.join(str_row), file=out_handler)


    def build_tables(self):
        for sample_id, (fasta_path, rdp_path) in self.samples.items():
            sample_comp_table = []
            abundance_by_sequence = get_abundance_from_fasta(fasta_path)
            total_abd = sum(abundance_by_sequence.values())
            for rdp_line in read_rpd_file(rdp_path):
                sequence_id = rdp_line[0]
                abundance = abundance_by_sequence[sequence_id]
                normalized_abundance = round(abundance / total_abd * 100, self.float_precision)
                taxonomy = ';'.join(get_lineage(rdp_line))
                row = [taxonomy, sample_id, sequence_id, abundance, normalized_abundance]
                sample_comp_table.append(row)

            sample_comp_table = sorted(sample_comp_table,key=operator.itemgetter(0)) #sort sample by taxo
            self.comparaison_table.extend(sample_comp_table)


    def check_path_validity(self):
        for sample_id, (fasta_path, rdp_path) in self.samples.items():
            if not os.path.isfile(fasta_path):
                logger.fatal('Invalid fasta path (sample: %s):%s' % (sample_id, fasta_path))
                sys.exit('Invalid file')
            if not os.path.isfile(rdp_path):
                logger.fatal('Invalid RDP path (sample: %s):%s' % (sample_id, rdp_path))
                sys.exit('Invalid file')



def retrieve_samples_path(listing_file):
    """
    From a tabulated file, return an ordered dict with fasta and rdp path foreach sample
    as: { sample_id1:(fasta_path, rdp_path) ...}
    First col: sampleid
    Second col: fasta_path
    Third col: rdp_path
    """
    samples = collections.OrderedDict()
    with open(listing_file, 'r') as lst_handler:
        for line_number, line in enumerate(lst_handler):
            line = line.strip()
            if not line: break
            arr_line = line.split('\t')
            arr_line = [v.strip() for v in arr_line]
            if len(arr_line) != 3:
                logger.fatal("Wrong number of fields (line number:%s, file:%s)" % (line_number, listing_file))
                sys.exit("Wrong number of fields")
            sample_id, fasta_path, rdp_path = arr_line
            if sample_id in samples:
                logger.fatal("Duplicated sample_id (id:%s, file:%s)" % (sample_id, listing_file))
                sys.exit("Duplicated id")
            samples[sample_id] = (fasta_path,rdp_path)
        return samples


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--samples_file',
                        type=argparse.FileType('r'),
                        help="A tabulated file with one sample by row." \
                         "The first column contains the sample id (must be uniq)" \
                         "The second column contains the fasta path" \
                         "The third, the rdp path",
                         required=True)

    parser.add_argument('-o', '--ouput_table',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output table file')

    args = parser.parse_args()

    lst_path = os.path.abspath(args.samples_file.name)
    samples_path = retrieve_samples_path(lst_path)
    sample_collection = SampleCollection(samples_path)
    sample_collection.write_comparaison_table(args.ouput_table)
    
