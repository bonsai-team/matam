import os
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)


SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from compute_abundance import abundance_calculation
import pytest

@pytest.fixture()
def expected_abundance():
    return {'scaff1': 1.75,
            'scaff2': 3.75,
            'scaff3': 1.25,
            'scaff4': 1.25
    }


def test_abundance_calculation(expected_abundance):
    blast_file = os.path.join(SAMPLE_DIR, 'scaffolds.blast')
    abundance = abundance_calculation(blast_file)

    assert set(abundance.keys()) == set(expected_abundance.keys())
    assert set(abundance.values()) == pytest.approx(set(expected_abundance.values()))
