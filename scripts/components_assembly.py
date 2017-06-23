#!/usr/bin/env python3
import os
import sys
import subprocess
import logging
from collections import defaultdict
import multiprocessing
import argparse

from fasta_utils import read_fasta_file_handle, format_seq
from fastq_utils import read_fastq_file_handle

from assembler_factory import AssemblerFactory

logger = logging.getLogger(__name__)

def extract_reads_by_component(fastq, read_metanode_component_filepath):
    if not os.path.isfile(fastq):
        logger.fatal('The input reads file does not exists:%s' % fastq)
        sys.exit("An error occured. Can't extract component's reads")

    if not os.path.isfile(read_metanode_component_filepath):
        logger.fatal('The file storing the correspondance between reads and components does not exists:%s' % read_metanode_component_filepath)
        sys.exit("An error occured. Can't extract component's reads")

    # Reading read --> component file
    logger.debug('Reading read-->component from {}'.format(read_metanode_component_filepath))
    read_component_dict = dict()
    with open(read_metanode_component_filepath, 'r') as read_metanode_component_fh:
        read_component_dict = {t[0]:t[2] for t in (l.split() for l in read_metanode_component_fh) if t[2] != 'NULL'}

    # Storing reads for each component
    logger.debug('Storing reads by component from {}'.format(fastq))
    component_reads_dict = defaultdict(list)
    with open(fastq, 'r') as fastq_fh:
        for header, seq, qual in read_fastq_file_handle(fastq_fh):
            try:
                component_reads_dict[read_component_dict[header]].append(tuple((header, seq, qual)))
            except KeyError:
                pass
    return component_reads_dict


def extract_lca_by_component(components_lca_filepath):
    if not os.path.isfile(components_lca_filepath):
        logger.fatal('The file storing the correspondance between lca info and components does not exists:%s' % components_lca_filepath)
        sys.exit("An error occured. Can't extract component's lca")

    # Reading components LCA and storing them in a dict
    logger.debug('Reading components LCA assignment from {0}'.format(components_lca_filepath))
    component_lca_dict = dict()
    with open(components_lca_filepath, 'r') as component_lca_fh:
        component_lca_dict = {t[0]:t[1] for t in (l.split() for l in component_lca_fh) if len(t) == 2}
    return component_lca_dict


def save_components(component_reads_dict, directory):
    """
    Take a dict (key=component_id, value=[header,seq,qual])
    and save each component to a file into the given directory:
    directory/
        component_%s.fq % component_id

    Return a dict (key=component_id, value=fastq_path)
    """

    try:
        os.mkdir(directory)
    except FileExistsError as fee:
        pass

    components_fq = {}
    for component_id, reads_list in component_reads_dict.items():
        fq_name = "component%s_reads.fq" % component_id
        fq_path = os.path.join(directory, fq_name)
        #logger.debug("Save component %s into %s" % (component_id, fq_path))
        components_fq[component_id] = fq_path

        with open(fq_path, 'w') as component_fh:
            for (header, seq, qual) in reads_list:
                component_fh.write('@{}\n{}\n+\n{}\n'.format(header, seq, qual))
    return components_fq


def isfastq(filepath):
    """
    Determine if filepath is a fastq file based on the extension
    """

    suffix = ('.fq', '.fastq')
    return any( [ filepath.endswith(s) for s in suffix] )


def isfasta(filepath):
    """
    Determine if filepath is a fasta file based on the extension
    """

    suffix = ('.fa', '.fasta', '.fna')
    return any( [filepath.endswith(s) for s in suffix] )


def nucleotidic_number(fastx):
    """
    Compute the total number of nucleotides in the fastx file
    """

    parser = None
    if isfasta(fastx):
        parser = read_fasta_file_handle
    elif isfastq(fastx):
        parser = read_fastq_file_handle

    if parser is None:
        logger.fatal("Can't dertermine whether this file is a fasta or fastq")
        sys.exit('Aborting assembly step')


    count = 0
    with open(fastx, 'r') as fastx_handle:
        for rec in parser(fastx_handle):
            seq = rec[1]
            count += len(seq)
    return count

def estimate_coverage(reads_fq, contigs_fa):
    """
    estimated_cov = reads_nt/contigs_nt
    """

    reads_nt = nucleotidic_number(reads_fq)
    contigs_nt = nucleotidic_number(contigs_fa)
    estimated_cov = None

    try:
        estimated_cov = reads_nt/contigs_nt
        #logger.debug("Estimated coverage:%s" % estimated_cov)

    except ZeroDivisionError:
        #logger.warning("Can't estimate the coverage for 0 length contigs:%s" % contigs_fa)
        pass

    return estimated_cov


def assemble_component(assembler_name,
                       in_fastq, workdir,
                       read_correction, cpu, coverage_threshold):
    logger.debug('Assembling: %s' % in_fastq)
    assembler_factory = AssemblerFactory()
    assembler = assembler_factory.get(assembler_name)
    assembler.build_command_line(in_fastq, workdir, read_correction, cpu)

    fasta_file = assembler.run()
    estimated_cov = estimate_coverage(in_fastq, fasta_file)
    logger.debug("Estimated coverage:%s, %s" % (estimated_cov, in_fastq))

    # Re-run the assembly with error correction activated
    if estimated_cov is not None and estimated_cov > coverage_threshold:
        assembler.build_command_line(in_fastq, workdir, 'yes', cpu)
        fasta_file = assembler.run()
        estimated_cov2 = estimate_coverage(in_fastq, fasta_file)
        logger.debug("Estimated coverage, before:%s, after:%s, %s" % (estimated_cov, estimated_cov2, in_fastq))

    return fasta_file


def concat_components_fasta_with_lca(assembled_components_fasta, contigs_fasta, component_lca_dict):
    if os.path.isfile(contigs_fasta):
        #logger.debug("Remove old contig fasta file:%s" % contigs_fasta)
        os.unlink(contigs_fasta)

    component_lca = 'NULL'
    contigs_fh = open(contigs_fasta, 'w')
    contig_count = 0
    for component_id, component_fasta in assembled_components_fasta.items():
        if component_id in component_lca_dict:
            component_lca = component_lca_dict[component_id]
        with open(component_fasta, 'r') as assembler_contigs_fh:
            for header, seq in read_fasta_file_handle(assembler_contigs_fh):
                if len(seq):
                    contig_count += 1
                    contigs_fh.write('>{0} component={1} '.format(contig_count, component_id))
                    contigs_fh.write('lca={0}\n{1}\n'.format(component_lca, format_seq(seq)))

    contigs_fh.close()


def _get_workdir(fq):
    """
    Convenient function to build a workdir from the name of the fastqfile
    """
    assembler_wkdir_basename, _ = os.path.splitext(fq)
    return '%s_assembly_wkdir' % assembler_wkdir_basename


def assemble_all_components(assembler_name,
                            fastq, read_metanode_component_filepath, components_lca_filepath,
                            out_contigs_fasta, workdir,
                            cpu, read_correction, coverage_threshold=20):


    logger.info("Save components to fastq files")
    components_dict = extract_reads_by_component(fastq, read_metanode_component_filepath)
    components_reads_fq = save_components(components_dict, workdir).items()
    assembled_components_fasta = {}

    logger.info("Assemble components")

    # Foreach component, build the parameters used by assemble_component and save
    # them into a list to be able to apply a map function on it
    params = []
    component_id_list = []
    for component_id, fq in components_reads_fq:
        params.append((assembler_name, fq, _get_workdir(fq), read_correction, 1, coverage_threshold))
        component_id_list.append(component_id)

    with multiprocessing.Pool(processes=cpu) as pool:
        try:
            fasta_list = pool.starmap(assemble_component, params)
        except subprocess.CalledProcessError as cpe:
            logger.fatal('Command %s returned non-zero exit status %s' % (cpe.cmd, cpe.returncode))
            sys.exit('Components assembly step failed')

    # Make the correspondance between the component_id and the fasta file
    assembled_components_fasta = dict(zip(component_id_list, fasta_list))

    lca_dict = extract_lca_by_component(components_lca_filepath)
    logger.info("Pool components contigs into: %s" % out_contigs_fasta)
    concat_components_fasta_with_lca(assembled_components_fasta,
                                     out_contigs_fasta, lca_dict)


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description='Assemble components')
    parser.add_argument('-a', '--assembler',
                        choices=[a.name() for a in AssemblerFactory.ASSEMBLER_ENGINES],
                        help="Select the assembler to be used. Default is %(default)s",
                        default="SGA")
    parser.add_argument('-i', '--input_fastq',
                        type=argparse.FileType('r'),
                        help='input fastq file',
                        required=True)
    parser.add_argument('-m', '--reads_metanode',
                        type=argparse.FileType('r'),
                        help='This file makes the correspondance between reads and the components',
                        required=True)
    parser.add_argument('-l', '--components_lca',
                        type=argparse.FileType('r'),
                        help='This file make the correspondance between lca and the components',
                        required=True)
    parser.add_argument('-w', '--workdir',
                        action = 'store',
                        type = str,
                        help = 'Working  directory')
    parser.add_argument('-o', '--output_fasta',
                        help='output fasta file',
                        required=True)
    parser.add_argument('--cpu',
                        action = 'store',
                        metavar = 'CPU',
                        type = int,
                        default = 1,
                        help = 'Max number of CPU to use. '
                        'Default is %(default)s cpu')
    parser.add_argument('--read_correction',
                        action = 'store',
                        type = str,
                        choices = ['no', 'yes', 'auto'],
                        default = 'no',
                        help = 'Set the assembler read correction step. '
                        'Default is %(default)s')

    args = parser.parse_args()
    assemble_all_components(args.assembler,
                            args.input_fastq.name, args.reads_metanode.name, args.components_lca.name,
                            args.output_fasta, args.workdir,
                            args.cpu, args.read_correction)
