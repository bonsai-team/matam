import os
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)


SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from krona import rdp_file_to_krona_text_file, make_krona_plot
from binary_utils import Binary

import pytest

@pytest.fixture()
def abundance():
    return {'32': 226,
            '45': 230,
            '3': 192,
            '33': 294,
            '11': 147,
            '49': 411,
            '31': 383,
            '51': 512,
            '38': 700,
            '15': 475,
            '48': 567,
            '47': 700,
            '12': 695,
            '35': 477,
            '43': 505,
            '34': 722,
            '40': 750,
            '39': 737,
            '14': 746,
            '16': 682,
            '8': 750,
            '4': 749,
            '87': 860
    }

def test_rdp_file_to_krona_text_file_without_abundance():
    rdp_file = os.path.join(SAMPLE_DIR, 'rdp.txt')
    krona_text_file = tempfile.NamedTemporaryFile()
    rdp_file_to_krona_text_file(rdp_file, krona_text_file.name)

    assert os.path.getsize(krona_text_file.name) > 0
    with open(krona_text_file.name,'r') as h:
        lines = h.readlines()
        assert set([int(l.split('\t')[0]) for l in lines]) == set([1])
        assert len(lines) == 23

def test_rdp_file_to_krona_text_file_with_abundance(abundance):
    rdp_file = os.path.join(SAMPLE_DIR, 'rdp.txt')
    krona_text_file = tempfile.NamedTemporaryFile()
    rdp_file_to_krona_text_file(rdp_file, krona_text_file.name, abundance=abundance)

    assert os.path.getsize(krona_text_file.name) > 0
    with open(krona_text_file.name,'r') as h:
        lines = h.readlines()
        assert len(lines) == 23
        assert set([int(l.split('\t')[0]) for l in lines]) == set(abundance.values())

def test_make_krona_plot():
    krona_bin = Binary.which('ktImportText')
    krona_file = os.path.join(SAMPLE_DIR, 'krona.txt')
    krona_html_file = tempfile.NamedTemporaryFile()

    make_krona_plot(krona_bin, krona_file, krona_html_file.name)
    assert os.path.getsize(krona_html_file.name) > 0
