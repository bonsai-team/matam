import os
import sys
import tempfile
import subprocess
from collections import defaultdict
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
DB_DIR = os.path.join(CURRENT_DIR, '..', 'db')
sys.path.append(SCRIPTS_DIR)


SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from sample_sam_by_coverage import read_fasta_file_handle
from sample_sam_by_coverage import read_tab_file_handle_sorted
from sample_sam_by_coverage import get_reads_by_pos, compute_depth
from sample_sam_by_coverage import sample_by_depth
import pytest


def _get_theo_cov_by_ref(samfile):
    """
    this the depth given by samtools depth -a
    """
    cmd = 'samtools depth -a %s' % samfile
    str_depth = subprocess.check_output(cmd, shell=True).decode('utf8')
    cov_by_ref = defaultdict(list)
    for line in str_depth.split('\n'):
        line = line.strip()
        if not line: continue
        ref_id, pos, depth = line.split('\t')
        depth = int(depth)
        cov_by_ref[ref_id].append(depth)
    return cov_by_ref


def _get_cov_by_ref(samfile, fasta_file):
    """
    this the depth given by our method
    """
    ref_seq_dict = {}
    with open(fasta_file, 'r') as fasta:
        for header, seq in read_fasta_file_handle(fasta):
            seqid = header.split()[0]
            ref_seq_dict[seqid] = seq.upper()

    cov_by_ref = {}
    with open(samfile, 'r') as sam:
        for alignment_tabs_list in read_tab_file_handle_sorted(sam, 2):
            ref_id = alignment_tabs_list[0][2]
            ref_len = len(ref_seq_dict[ref_id])
            reads_by_pos = get_reads_by_pos(alignment_tabs_list, ref_len)
            ref_cov = compute_depth(reads_by_pos)
            cov_by_ref[ref_id] = ref_cov
    return cov_by_ref

@pytest.mark.skipif(not os.path.isdir(DB_DIR),
                    reason= 'DB_DIR is missing:%s' % DB_DIR)
@pytest.mark.parametrize('samfile , fasta_ref_file', [
    [ os.path.join(SAMPLE_DIR, 'simple.sam'), os.path.join(SAMPLE_DIR, 'simple_sam_ref.fa')],
    [ os.path.join(SAMPLE_DIR, '16sp.sam'), os.path.join(DB_DIR, 'SILVA_128_SSURef_NR95.clustered.fasta')]
], indirect=False)
def test_coverage_calculation(samfile, fasta_ref_file):
    """
    assert that our coverage calculation method is the same as
    samtools depth -a
    """
    theoretical_cov = _get_theo_cov_by_ref(samfile)
    computed_cov = _get_cov_by_ref(samfile, fasta_ref_file)

    assert theoretical_cov == computed_cov


@pytest.mark.skipif(not os.path.isdir(DB_DIR),
                    reason= 'DB_DIR is missing:%s' % DB_DIR)
def test_max_cov():
    """
    assert that the sampled result respect the given threshold
    """
    samfile = os.path.join(SAMPLE_DIR, '16sp.sam')
    fasta_file = os.path.join(DB_DIR, 'SILVA_128_SSURef_NR95.clustered.fasta')

    threshold = 40
    with open(samfile,'r') as sam_handler,\
         open(fasta_file, 'r') as  fasta_ref_handler, \
         tempfile.NamedTemporaryFile(mode='w') as out_sam_handler:

        sample_by_depth(sam_handler, fasta_ref_handler, threshold, out_sam_handler)
        cov_by_ref = _get_cov_by_ref(out_sam_handler.name, fasta_file)
        for ref, cov_list in cov_by_ref.items():
            assert all( [ cov <= threshold for cov in cov_list] )
