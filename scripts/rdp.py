#!/usr/bin/env python3

import subprocess
import logging


logger = logging.getLogger(__name__)

def run_rdp_classifier(rdp_exe, in_fasta, out_classification_file, cutoff=0.9, gene='16srrna'):

    parameters = { 'fa': in_fasta, 'out': out_classification_file, 'cutoff': cutoff, 'gene': gene }
    cmd_line = '{rdp_exe} classify -c {cutoff} -f allrank -g {gene} -o {out} {fa}'.format(rdp_exe=rdp_exe, **parameters)
    logger.debug('CMD: {}'.format(cmd_line))
    rc = subprocess.call(cmd_line, shell=True, bufsize=0)
    return rc
