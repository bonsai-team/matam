import os
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)

SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from rdp import run_rdp_classifier, read_rpd_file, get_lineage, filter_rdp_file
from binary_utils import Binary


def test_run_rdp_classifier_ok():
    fasta = os.path.join(SAMPLE_DIR, 'scaffolds.fa')
    jar = Binary.which('classifier.jar')

    if jar is not None:
        java = Binary.assert_which('java')
        rdp_exe = '{java} -Xmx1g -jar {jar}'.format(java=java, jar=jar)
    else:
        rdp_exe = Binary.assert_which('classifier')

    result_file = tempfile.NamedTemporaryFile()
    run_rdp_classifier(rdp_exe, fasta, result_file.name)


def test_read_rdp_results():
    rdp_file = os.path.join(SAMPLE_DIR, 'rdp.txt')
    lines = list(read_rpd_file(rdp_file))
    assert len(lines) == 23  # 23 scaffolds


def test_filter_rdp_file():
    rdp_file = os.path.join(SAMPLE_DIR, 'rdp.txt')
    result_file = tempfile.NamedTemporaryFile()
    filter_rdp_file(rdp_file, result_file.name)
    for line in read_rpd_file(result_file.name):
        assert len(line) == 19  # seqid + 6 taxonomic levels * 3
        if line[0] == "87":
            assert get_lineage(line) == ['unclassified'] * 6
