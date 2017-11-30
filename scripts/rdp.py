#!/usr/bin/env python3

import runner
import logging
import re
import sys

logger = logging.getLogger(__name__)

def run_rdp_classifier(rdp_exe, in_fasta, out_classification_file, cutoff=0.8, gene='16srrna'):

    parameters = { 'fa': in_fasta, 'out': out_classification_file, 'cutoff': cutoff, 'gene': gene }
    cmd_line = '{rdp_exe} classify -c {cutoff} -f fixrank -g {gene} -o {out} {fa}'.format(rdp_exe=rdp_exe, **parameters)
    runner.logged_check_call(cmd_line)


def filter_rdp_file(rdp_file, out_fltr_rdp_file, cutoff=0.8):
    """
    For a given line, if any level has a confidence score under cutoff:
        tag all nodes as "unclassified"
    """

    with open(out_fltr_rdp_file, 'w') as rdp_out:
        for tab in read_rpd_file(rdp_file):
            seqid = tab.pop(0)
            rank_names, rank_levels, rank_scores = zip(*[tab[i:i+3] for i in range(0,len(tab),3)])
            if any([float(score) < cutoff for score in rank_scores]):
                rank_names = ['unclassified'] * 6

            rdp_line = [seqid] + ['\t'.join(triplet) for triplet in zip(rank_names, rank_levels, rank_scores)]
            print('\t'.join(rdp_line), file=rdp_out)


def read_rpd_file(rdp_path):
    """
    parse a rdp file and return a generator
    """
    with open(rdp_path, 'r') as in_rdp_handler:
        for l in in_rdp_handler:
            l = l.strip()
            if not l or l.startswith('#'): continue
            rdp_line = re.split('"?\t+"?', l)
            rdp_line = [field.strip() for field in rdp_line]
            if len(rdp_line) != 19: # seqid + 6 taxonomic levels * 3
                logger.fatal('RDP: wrong number of fields -- %s, expected 19, line: %s' % (len(rdp_line), l))
                sys.exit('Failed to parse RDP file:%s' % rdp_path)

            yield rdp_line


def get_lineage(splitted_rdp_line):
    """
    Extract lineage from splitted rdp line
    """

    splitted_rdp_line.pop(0) # remove id
    return [ splitted_rdp_line[i] for i in range(0, len(splitted_rdp_line), 3) ]