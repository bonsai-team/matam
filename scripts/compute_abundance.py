#!/usr/bin/env python3

import os
import sys
import logging
import runner
import collections
import tempfile
import shutil
import re

from fasta_clean_name import read_fasta_file_handle, format_seq

logger = logging.getLogger(__name__)

def index_ref(indexdb_bin_path, input_fasta_ref_path, output_basepath, max_mem, verbose=False):
    logger.info('Indexing %s' % input_fasta_ref_path)
    parameters = { 'input': input_fasta_ref_path, 'output': output_basepath, 'max_mem': max_mem }
    cmd_line = '{bin} --ref {input},{output} -m {max_mem}'.format(bin=indexdb_bin_path, **parameters)

    if verbose:
        cmd_line += ' -v'

    runner.logged_check_call(cmd_line, verbose=verbose)


def reads_mapping(sortmerna_bin, fasta_ref_path, index_ref_basepath, reads_path, output_basepath, best, min_lis, evalue, cpu, verbose=False):

    logger.info('=== Reads mapping against %s ===' % index_ref_basepath)
    parameters = { 'fasta_ref': fasta_ref_path, 'index_ref': index_ref_basepath, 'reads': reads_path, 'output':output_basepath,
                   'best': best, 'min': min_lis, 'evalue': evalue, 'cpu': cpu }

    cmd_line = '{bin} --ref {fasta_ref},{index_ref} --reads {reads} --aligned {output}' \
               ' --fastx --sam --blast "1" --log --best {best} --min_lis {min}' \
               ' -e {evalue:.2e} -a {cpu}'.format(bin=sortmerna_bin, **parameters)
    if verbose:
        cmd_line += ' -v'

    runner.logged_check_call(cmd_line, verbose=verbose)


def get_best_matches(best_matches_bin, input_blast_path, out_blast_path, max_mem, cpu):

    outdir = os.path.dirname(out_blast_path)

    sort_bin = 'sort -T ' + outdir + ' -S ' + str(max_mem)
    sort_bin += 'M --parallel ' + str(cpu)

    cmd_line = sort_bin + ' -k1,1V -k12,12nr ' + input_blast_path
    cmd_line += ' | ' + best_matches_bin + ' -p 0.99 -o ' + out_blast_path

    runner.logged_check_call(cmd_line)


def abundance_calculation(blast_path):
    scaffolds_by_read = collections.defaultdict(list)
    reads_by_scaffold = collections.defaultdict(list)

    with open(blast_path, 'r') as f:
        for line in f:
            if line.strip() == '': continue
            read, scaffold = [ s.strip() for s in line.split('\t')[:2] ]
            scaffolds_by_read[read].append(scaffold)
            reads_by_scaffold[scaffold].append(read)

    #asign a weight for each read
    read_weight = {}
    for read, scaffolds in scaffolds_by_read.items():
        uniq_scaffolds = set(scaffolds)
        if len(uniq_scaffolds) == 0:
            weight = 0
        elif len(uniq_scaffolds) == 1:
            weight = 1
        else:
            weight = 1/len(uniq_scaffolds)
        read_weight[read] = weight

        if len(scaffolds) != len(uniq_scaffolds):
            # keep scaffolds names where this read map more than once
            more_than_once = {k:v for k,v in collections.Counter(scaffolds).items() if v > 1}
            logger.warning('%s is mapped more than once on the same scaffold (%s) but it will contribute \
to the abundance of this scaffold only as 1 weight where weight=1/uniq_scaffolds_nb=%s' % (read, more_than_once, weight))

    # compute abundance for each scaffold depending on read weight
    abundance_by_scaffold = collections.defaultdict(lambda: 0)
    for scaffold, reads in reads_by_scaffold.items():
        for read in set(reads):
            abundance_by_scaffold[scaffold]+= read_weight[read]

        abundance_by_scaffold[scaffold] = round(abundance_by_scaffold[scaffold], 2)
    return abundance_by_scaffold


def get_abundance_by_scaffold(idx_bin, map_bin, best_bin,
                              input_fasta_ref, input_fasta_reads,
                              best=10, min_lis=10, evalue=1e-05,
                              max_mem=10000, cpu=4,
                              output_dir_basepath="/tmp/",
                              verbose=False,
                              keep_tmp=False):

    output_dir_basepath = os.path.join(output_dir_basepath , '') # add a trailing slash
    outdir = tempfile.mkdtemp(dir=output_dir_basepath, prefix='abundance_')

    #index ref
    idx_ref_basepath = os.path.join(outdir, 'idx_prefix')

    index_ref(idx_bin, input_fasta_ref, idx_ref_basepath, max_mem, verbose=verbose)

    #reads mapping
    filtered_basepath = os.path.join(outdir, 'filt_prefix')
    reads_mapping(map_bin, input_fasta_ref, idx_ref_basepath, input_fasta_reads, filtered_basepath, best, min_lis, evalue, cpu, verbose=verbose)

    #best matches
    blast_path = '%s.blast' % filtered_basepath
    best_path = '%s.best' % blast_path
    get_best_matches(best_bin, blast_path, best_path, max_mem, cpu)

    #abundance calculation
    abundance = abundance_calculation(best_path)

    #delete tempfiles
    if not keep_tmp:
        shutil.rmtree(outdir)

    return abundance


def complete_fasta_with_abundance(input_fasta, output_fasta, abundance):
    in_fasta_handler = open(input_fasta, 'r')
    out_fasta_handler = open(output_fasta, 'w')

    for header, seq in read_fasta_file_handle(in_fasta_handler):
        id = header.split()[0].strip()
        if id not in abundance:
            logger.fatal("Can't find the abundance for:%s" % id)
            sys.exit("Can't find abundance")
        else:
            ab = abundance[id]
            header = '{header} count={abundance}'.format(header=header, abundance=ab)
            out_fasta_handler.write( '>{header}\n{seq}\n'.format(header=header, seq=format_seq(seq)) )

    in_fasta_handler.close()
    out_fasta_handler.close()


def get_abundance_from_fasta(fasta, regexp='count=(\d+\.\d+|\d+)'):
    abundance = {}
    in_fasta_handler = open(fasta, 'r')
    p = re.compile(regexp)
    for i, [header, _] in enumerate(read_fasta_file_handle(in_fasta_handler)):
        id = header.split()[0].strip()
        m = re.search(p, header)
        if not m:
            logger.fatal("Can't retrieve abundance information:\nfasta_path:%s\nfasta_header:%s" % (fasta, header))
            sys.exit("Cant't retrieve abundance")
        abundance[id] = float(m.group(1))
    in_fasta_handler.close()
    return abundance
