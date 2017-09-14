#!/usr/bin/env python3

import sys
import re
import random
import argparse



def read_tab_file_handle_sorted(tab_file_handle, factor_index=0):
    """
    Parse a tab file (sorted by a column) and return a generator
    """
    previous_factor_id = ''
    factor_tab_list = list()
    # Reading tab file
    for line in tab_file_handle:
        #ignore header lines
        if line.startswith('@'): continue
        l = line.strip()
        #ignore blanck lines
        if not l: continue
        tab = l.split()
        current_factor = tab[factor_index]
        # Yield the previous factor tab list
        if current_factor != previous_factor_id:
            if previous_factor_id:
                yield factor_tab_list
                factor_tab_list = list()
        factor_tab_list.append(tab)
        previous_factor_id = current_factor
    # Yield the last tab list
    yield factor_tab_list
    # Close tab file
    tab_file_handle.close()


def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in fasta_file_handle:
        if line[0] == '>':
            # Yield the last read header and sequence
            if sequence_nb:
                yield (header, ''.join(seqlines))
                del seqlines[:]
            # Get header
            header = line[1:].rstrip()
            sequence_nb += 1
        else:
            # Concatenate sequence
            seqlines.append(line.strip())
    # Yield the input file last sequence
    yield (header, ''.join(seqlines))
    # Close input file
    fasta_file_handle.close()


def decode_cigar(cigar_values, current_pos=0):
    """
    Return a list of positions (on the ref) covered by this CIGAR (read)
    """
    covered_pos = set()
    for c_val in cigar_values:
        c_len = int(c_val[:-1])
        c_letter = c_val[-1].upper()
        if c_letter in ('M','X', '='):
            start = current_pos
            end = start + c_len - 1
            covered_pos |= set(range(start, end + 1))
            current_pos += c_len
        elif c_letter in ('N', 'D'):
            current_pos += c_len
    return covered_pos


def split_cigar(cigar_string):
    """
    Return all cigar's components
    One cigar component is a number of any digit, followed by a letter or =
    """
    cigar_pattern = re.compile("[\\d]+[a-zA-Z|=]")
    cigar_elems = re.findall(cigar_pattern, cigar_string)
    return cigar_elems


def compute_depth(reads_by_pos, discarded_reads=None):
    if discarded_reads is None:
        discarded_reads=[]
    return [len([r for r in reads if r not in discarded_reads]) for reads in reads_by_pos]


def restrict_depth(reads_by_pos, threshold):
    """
    Return the name of reads to be discarded to fulfill
    the depth threshold
    """
    discarded_reads = set()
    for pos, reads in enumerate(reads_by_pos):
        depth = len([r for r in reads if r not in discarded_reads])
        while depth > threshold:
            #randomly pick a read
            read_id = random.choice(reads)
            #pick a reads not previously picked
            if read_id in discarded_reads: continue
            discarded_reads.add(read_id)
            depth = len([r for r in reads if r not in discarded_reads])
    return discarded_reads


def get_reads_by_pos(alignment_tabs_list, ref_len):
    """
    Take a list of alignments (SAM alignment lines)
    and return foreach pos on the reference, the list
    of reads constributing to this pos
    """

    #reads_by_pos = [ set() for p in range(ref_len) ]
    #do not use a set because the same exact read can be mapped at the same pos more than once
    #see https://github.com/biocore/sortmerna/issues/137
    #Init an empty list foreach pos of the reference
    cigar_cache = {}
    reads_by_pos = [ list() for p in range(ref_len) ]
    for alignment_tab in alignment_tabs_list:
        read_id = alignment_tab[0]
        cigar = alignment_tab[5]
        #get zero-based position
        leftmost_mapping_pos = int(alignment_tab[3]) - 1
        #only mapped reads are considered
        if leftmost_mapping_pos < 0: continue
        if cigar in cigar_cache:
            mapping_positions = cigar_cache[cigar]
        else:
            mapping_positions = decode_cigar(split_cigar(cigar))
            cigar_cache[cigar] = mapping_positions

        for pos in mapping_positions:
            reads_by_pos[leftmost_mapping_pos + pos].append(read_id)
    return reads_by_pos


def sample_by_depth(sam_handler, fasta_ref_handler, threshold, sampled_out_sam_handler):
    """
    sample the SAM file based on the depth at each pos
    """
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(fasta_ref_handler):
            seqid = header.split()[0]
            ref_seq_dict[seqid] = seq.upper()

    # Reading sam file reference by reference
    for alignment_tabs_list in read_tab_file_handle_sorted(sam_handler, 2):
        ref_id = alignment_tabs_list[0][2]
        ref_len = len(ref_seq_dict[ref_id])
        #print("@SQ	SN:%s	LN:%s" % (ref_id, ref_len), file=sys.stderr)
        reads_by_pos = get_reads_by_pos(alignment_tabs_list, ref_len)
        reads_to_discard = restrict_depth(reads_by_pos, threshold)
        for alignment_tab in alignment_tabs_list:
            read_id = alignment_tab[0]
            cigar = alignment_tab[5]
            leftmost_mapping_pos = int(alignment_tab[3]) - 1
            if read_id not in reads_to_discard:
                print('{}\n'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler, end='')


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_sam
    parser.add_argument('-i', '--input_sam',
                        metavar='INSAM',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input sam file, sorted by reference and position')
    # -o / --output_sam
    parser.add_argument('-o', '--output_sam',
                        metavar='OUTSAM',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output filtered sam file')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar='REF',
                        type=argparse.FileType('r'),
                        required=True,
                        help='References fasta file')
    # -c / --cov_threshold
    parser.add_argument('-c', '--cov_threshold',
                        metavar='COV',
                        type=int,
                        default=500,
                        help='Identity threshold. '
                             'Default is %(default)s')

    args = parser.parse_args()
    sample_by_depth(args.input_sam, args.references, args.cov_threshold, args.output_sam)
