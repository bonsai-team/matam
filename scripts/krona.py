#!/usr/bin/env python3

import sys
import subprocess
import logging

from fasta_clean_name import read_fasta_file_handle
from compute_abundance import get_abundance_from_fasta
from rdp import read_rpd_file, get_lineage

logger = logging.getLogger(__name__)


def rdp_file_to_krona_text_file(rdp_file, krona_text_file, abundance=None):

    out_krona_handler = open(krona_text_file, 'w')
    for rdp_line in read_rpd_file(rdp_file):
        seq_id = rdp_line[0]
        count = 1
        if abundance is not None:
            if seq_id not in abundance:
                logger.fatal("No abundance for this seq_id:%s" % seq_id)
                sys.exit("No abundance for this seq_id:%s" % seq_id)
            count = abundance[seq_id]

        lineage = get_lineage(rdp_line)
        out_line = '{abundance}\t{lineage}\n'.format(abundance=count,
                                                     lineage = '\t'.join(lineage))

        out_krona_handler.write(out_line)
    out_krona_handler.close()


def make_krona_plot(krona_bin, krona_text_file, krona_html_file):
    logger.info('Make krona plot with:%s' % krona_text_file)
    cmd_line = '{bin} {txt} -o {html}'.format(bin=krona_bin, txt=krona_text_file, html=krona_html_file)

    logger.debug('CMD: {}'.format(cmd_line))
    rc = subprocess.call(cmd_line, shell=True)
    return rc
