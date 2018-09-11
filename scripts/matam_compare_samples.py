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
        self.samples_path = samples_path
        #keep the same order to generate the tables
        self.samples_id = [k for (k,v) in samples_path.items()] #no keys method for ordereddict
        self.comparaison_table = [('Taxonomy', 'SampleID', 'SequenceID', 'Abundance', 'Normalized_Abundance(%)')]
        self.contingency_table = [['Taxonomy/Samples', *self.samples_id]]

        if not samples_path:
            logger.fatal('The sample collection is empty')
            sys.exit('Empty collection')

        self._check_path_validity()
        #contingency table is build from the comparaison table
        self._build_comparaison_table()
        self._build_contingency_table()


    def write_comparaison_table(self, out_handler):
        self._write_table(self.comparaison_table, out_handler)


    def write_contingency_table(self, out_handler):
        self._write_table(self.contingency_table, out_handler)

    def _write_table(self, table, out_handler):
        for row in table:
            str_row = [str(v) for v in row]
            print('\t'.join(str_row), file=out_handler)

    def _groupby(self, table, operator):
        return itertools.groupby(sorted(table,key=operator), operator)


    def _build_comparaison_table(self):
        taxonomic_col = operator.itemgetter(0)
        sample_col = operator.itemgetter(1)

        for sample_id, (fasta_path, rdp_path) in self.samples_path.items():
            sample_comp_table = []
            abundance_by_sequence = get_abundance_from_fasta(fasta_path)
            total_abundance = sum(abundance_by_sequence.values())
            for rdp_line in read_rpd_file(rdp_path):
                sequence_id = rdp_line[0]
                abundance = abundance_by_sequence[sequence_id]
                normalized_abundance = round(abundance / total_abundance * 100, self.float_precision)
                taxonomy = ';'.join(get_lineage(rdp_line))
                row = [taxonomy, sample_id, sequence_id, abundance, normalized_abundance]
                sample_comp_table.append(row)

            sample_comp_table = sorted(sample_comp_table, key=taxonomic_col)
            self.comparaison_table.extend(sample_comp_table)


    def _build_contingency_table(self):
        taxonomic_col = operator.itemgetter(0)
        sample_col = operator.itemgetter(1)
        unpivoted_contingence_table = []

        # Build long format contingency table
        for sample, sample_group in self._groupby(self.comparaison_table[1:],sample_col):
            sample_group_list = list(sample_group)
            total_abundance = sum([ row[3] for row in sample_group_list ])

            for taxonomy, taxonomic_group in self._groupby(sample_group_list, taxonomic_col):
                taxonomic_group_list = list(taxonomic_group)
                taxonomic_abundance = sum([ row[3] for row in taxonomic_group_list ])
                normalized_taxo_abd = round(taxonomic_abundance / total_abundance * 100, self.float_precision)
                row = [taxonomy, sample, normalized_taxo_abd]
                unpivoted_contingence_table.append(row)

        # Pivot the table to wide format
        for taxonomy, taxonomic_group in self._groupby(unpivoted_contingence_table, taxonomic_col):
            row = [taxonomy]
            abd_by_sample = dict([(row[1], row[2]) for row in taxonomic_group])
            row.extend( [abd_by_sample.get(sample, None) for sample in self.samples_id])
            self.contingency_table.append(row)


    def _check_path_validity(self):
        for sample_id, (fasta_path, rdp_path) in self.samples_path.items():
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
            if not line: continue
            arr_line = line.split('\t')
            arr_line = [v.strip() for v in arr_line]
            if len(arr_line) != 3:
                logger.fatal("Wrong number of fields (line number:%s, file:%s)" % (line_number, listing_file))
                sys.exit("Wrong number of fields")
            sample_id, fasta_path, rdp_path = arr_line
            if sample_id in samples:
                logger.fatal("Duplicated sample_id (id:%s, file:%s)" % (sample_id, listing_file))
                sys.exit("Duplicated id")
            samples[sample_id] = (os.path.expanduser(fasta_path),
                                  os.path.expanduser(rdp_path))
        return samples


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    # Arguments parsing
    parser = argparse.ArgumentParser(description='This script let you compare two or more samples coming from MATAM -- v1.5.1 or superior')

    parser.add_argument('-s', '--samples_file',
                        type=argparse.FileType('r'),
                        help="A tabulated file with one sample by row. "
                        "The first column contains the sample id (must be unique) "
                        "The second column contains the fasta path. The abundances must be present into this file. "
                        "The third, the rdp path. "
                        "Paths can be absolute or relative to the current working directory.",
                         required=True)

    parser.add_argument('-t', '--ouput_comparaison_table',
                        type=argparse.FileType('w'),
                        help='Output a table with the abundance for each sequence',
                        required=True)

    parser.add_argument('-c', '--ouput_contingency_table',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output a contingency table (taxonomy vs samples)',
                        required=True)

    args = parser.parse_args()


    logger.info("Parsing the sample list")
    lst_path = os.path.abspath(args.samples_file.name)
    samples_path = retrieve_samples_path(lst_path)

    logger.info("Build the tables")
    sample_collection = SampleCollection(samples_path)

    logger.info("Write the tables")
    sample_collection.write_comparaison_table(args.ouput_comparaison_table)
    sample_collection.write_contingency_table(args.ouput_contingency_table)

    logger.info("Done")
