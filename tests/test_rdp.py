import os
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)

SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from rdp import run_rdp_classifier
from binary_utils import Binary

def test_run_rdp_classifier_ok():
    fasta = os.path.join(SAMPLE_DIR, 'scaffolds.fa')
    bin = Binary.which('java')
    jar = Binary.which('classifier.jar')
    result_file = tempfile.NamedTemporaryFile()
    assert run_rdp_classifier(bin, jar, fasta, result_file.name) == 0
    with open(result_file.name, 'r') as h:
        lines = h.readlines()
        lines = [ l for l in lines if not l.startswith('#') ] #ignore comment lines
        lines = [ l for l in lines if l ] #ignore empty lines
        assert len(lines) == 22 # 22 scaffolds



