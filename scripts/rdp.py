#!/usr/bin/env python3

import subprocess
import logging
import re
import sys

logger = logging.getLogger(__name__)

def run_rdp_classifier(rdp_exe, in_fasta, out_classification_file, cutoff=0.9, gene='16srrna'):

    parameters = { 'fa': in_fasta, 'out': out_classification_file, 'cutoff': cutoff, 'gene': gene }
    cmd_line = '{rdp_exe} classify -c {cutoff} -f allrank -g {gene} -o {out} {fa}'.format(rdp_exe=rdp_exe, **parameters)
    logger.debug('CMD: {}'.format(cmd_line))
    rc = subprocess.call(cmd_line, shell=True, bufsize=0)
    return rc


def read_rpd_file(rdp_path):
    """
    parse a rdp file and return a generator
    """
    with open(rdp_path, 'r') as in_rdp_handler:
        for l in in_rdp_handler:
            l = l.strip()
            if l == '' or l.startswith('#'): continue
            rdp_line = re.split('"?\t+"?', l)
            rdp_line = [field.strip() for field in rdp_line]
            if rdp_line[1] != 'Root' or (len(rdp_line) - 1) % 3 != 0:
                logger.fatal('The RDP file is no well formatted, line: %s' % l)
                sys.exit('Failed to parse RDP file:%s' % rdp_path)
            yield rdp_line


def get_lineage(splitted_rdp_line):
    """
    Extract lineage from splitted rdp line
    """

    lineage = []
    splitted_rdp_line.pop(0) # remove id
    for i in range(0, len(splitted_rdp_line), 3):
        if splitted_rdp_line[i] == 'Root': continue
        lineage.append(splitted_rdp_line[i])
    return lineage