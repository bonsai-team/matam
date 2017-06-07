#!/usr/bin/env python3

import sys
import subprocess
import logging
import re

from fasta_clean_name import read_fasta_file_handle
from compute_abundance import get_abundance_from_fasta

logger = logging.getLogger(__name__)


def rdp_file_to_krona_text_file(rdp_file, krona_text_file, abundance=None):

    in_rdp_handler = open(rdp_file, 'r')
    out_krona_handler = open(krona_text_file, 'w')

    for l in in_rdp_handler:
        l = l.strip()
        if l == '' or l.startswith('#'): continue
        rdp_line = re.split('"?\t+"?', l)
        if rdp_line[1] != 'Root' or (len(rdp_line) - 1) % 3 != 0:
            logger.fatal('The RDP file is no well formatted, line: %s' % l)
            sys.exit('Failed to parse RDP file:%s' % rdp_file)
        id = rdp_line.pop(0)
        count = 1
        if abundance is not None:
            if id not in abundance:
                logger.fatal("No abundance for this id:%s" % id)
                sys.exit("No abundance for this id:%s" % id)
            count = abundance[id]

        lineage = []
        for i in range(0, len(rdp_line), 3):
            if rdp_line[i] == 'Root': continue
            lineage.append(rdp_line[i])
            #rank = rdp_line[i+1]
            #confidence = rdp_line[i+2]
        out_line = '{abundance}\t{lineage}\n'.format(abundance=count,
                                                     lineage = '\t'.join(lineage))

        out_krona_handler.write(out_line)

    in_rdp_handler.close()
    out_krona_handler.close()

def make_krona_plot(krona_bin, krona_text_file, krona_html_file):
    logger.info('Make krona plot with:%s' % krona_text_file)
    cmd_line = '{bin} {txt} -o {html}'.format(bin=krona_bin, txt=krona_text_file, html=krona_html_file)

    logger.debug('CMD: {}'.format(cmd_line))
    rc = subprocess.call(cmd_line, shell=True)
    return rc
